# WGCNA
library(WGCNA)

WGCNA_matrix <- t(expr[order(apply(expr, 1, mad), decreasing = T)[1:2000], ])
sel <- colnames(WGCNA_matrix)
datExpr0 <- expr[sel, ] %>%
  t() %>%
  as.data.frame()
datExpr <- datExpr0

sampleNames <- rownames(datExpr)
datTraits <- datTraits[, c("sample", names(datTraits)[-1] %>% str_sort())]
names(datTraits)[1] <- "sample"
traitRows <- match(sampleNames, datTraits$sample)
rownames(datTraits) <- datTraits[traitRows, 1]
datTraits <- datTraits %>% dplyr::select(-1)
identical(rownames(datTraits), rownames(datExpr))

powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

net <- blockwiseModules(
  datExpr,
  power = sft$powerEstimate,
  maxBlockSize = 6000,
  TOMType = "unsigned", minModuleSize = 30,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = F,
  verbose = 3
)
table(net$colors)

mergedColors <- labels2colors(net$colors)
table(mergedColors)

nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
moduleColors <- labels2colors(net$colors)
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10, 6)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")",
                    sep = ""
)
dim(textMatrix) <- dim(moduleTraitCor)

modNames <- substring(names(MEs), 3)
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")

TOM <- TOMsimilarityFromExpr(datExpr, power = sft$powerEstimate)
probes <- colnames(datExpr)
inModule <- (moduleColors == module)
modProbes <- probes[inModule]
modTOM <- TOM[inModule, inModule]
dimnames(modTOM) <- list(modProbes, modProbes)

cyt <- exportNetworkToCytoscape(
  modTOM,
  edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse = "-"), ".txt", sep = ""),
  nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse = "-"), ".txt", sep = ""),
  weighted = TRUE,
  threshold = 0.02,
  nodeNames = modProbes,
  nodeAttr = moduleColors[inModule]
)

# RWR
library(Matrix)
(load("PPI_Network_STRING_HINT.RData"))
genes <- fread("pathway.txt") %>% dplyr::filter(Category == "TCRsignalingPathway") %>% .$Symbol
x <- RandomWalkRestart(Network, intersect(genes, rownames(Network)))
res <- x %>% as.data.frame() %>% set_colnames("rank") %>% rownames_to_column("gene") %>% arrange(desc(rank))
write.csv(res, "rank_TCR.csv", row.names = F)

genes <- fread("pathway.txt") %>% dplyr::filter(Category == "BCRSignalingPathway") %>% .$Symbol
x <- RandomWalkRestart(Network, intersect(genes, rownames(Network)))
res <- x %>% as.data.frame() %>% set_colnames("rank") %>% rownames_to_column("gene") %>% arrange(desc(rank))
write.csv(res, "rank_BCR.csv", row.names = F)

genes <- fread("pathway.txt") %>% dplyr::filter(Category == "Antigen_Processing_and_Presentation") %>% .$Symbol
x <- RandomWalkRestart(Network, intersect(genes, rownames(Network)))
res <- x %>% as.data.frame() %>% set_colnames("rank") %>% rownames_to_column("gene") %>% arrange(desc(rank))
write.csv(res, "rank_APP.csv", row.names = F)

library(clusterProfiler)
library(org.Hs.eg.db)
library(GseaVis)
l <- list.files(pattern = "^rank_.*.csv")

FUN <- function(file) {
  x <- read.csv(file)
  df <- x %>% as.data.frame()
  mRNA_si <- df
  head(mRNA_si)
  names(mRNA_si)[1] <- "HUGO"
  names(mRNA_si)[2] <- "Weight"
  mRNA_si$Weight <- as.numeric(as.character(mRNA_si$Weight))
  mRNA_si <- mRNA_si[order(mRNA_si$Weight, decreasing = T), ]
  mRNA_si <- mRNA_si %>% dplyr::filter(!is.na(Weight))
  si.id <- mRNA_si$Weight
  si.id <- length(si.id):1
  names(si.id) <- mRNA_si$HUGO
  return(si.id)
}
all_glist <- lapply(l, FUN)

hallmark$ENTREZID <- as.character(hallmark$ENTREZID)
ID <- bitr(hallmark$ENTREZID, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
hallmark <- hallmark %>%
  inner_join(ID) %>%
  dplyr::select(-ENTREZID)

TERM2GENE <- hallmark[, c("KEGGID", "SYMBOL")]
TERM2NAME <- hallmark[, c("KEGGID", "DESCRIPTION")]

lapply(1:3, function(x){
  set.seed(2020)
  clustergsea <- GSEA(
    geneList = all_glist[[x]], TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, verbose = F,
    minGSSize = 0, maxGSSize = 500, nPerm = 1000, pvalueCutoff = 1
  )
  clustergsea <- clustergsea[clustergsea$pvalue < 0.05, asis = T]
  return(clustergsea)
}) -> m_gsea_list

for (proj in c("_APP.csv", "_TCR.csv", "_BCR.csv")) {
  df <- read.csv(str_c("rank", proj))
  p_sel <- df %>%
    dplyr::slice(1:(0.1 * (nrow(df)))) %>%
    .$rank
  top1 <- df %>% dplyr::filter(rank >= p_sel)
  genes <- top1$gene %>% unique()
  expr <- expr[genes, ] %>% na.omit()
  df <- read.csv("easier.csv", row.names = 1)
  expr <- expr[, rownames(df)]
  FUN <- function(i) {
    exprSet <- rbind(df[, i], expr)
    rownames(exprSet)[1] <- i
    l <- lapply(2:nrow(exprSet), FUN = FUN_cor, exprSet = exprSet)
    tmp_res <- do.call(rbind, l)
    rownames(tmp_res) <- rownames(expr)
    tmp_res <- data.frame(term = i, tmp_res %>% rownames_to_column("gene"))
    return(tmp_res)
  }
  l <- lapply(names(df), FUN = FUN)
  l2 <- lapply(l, function(i) {
    i %>%
      dplyr::filter(cor > 0.6 & p < 0.05) %>%
      .$gene
  })
  genes <- purrr::reduce(l2, intersect)
  df2 <- data.frame(gene = genes)
  write.csv(df2, str_c("top_immune_genes", proj), row.names = F)
}

# ICBnetIS
library(GSVA)
exprSet <- surv_expr %>% dplyr::select(-c("time", "status")) %>% column_to_rownames("sample") %>%
  t() %>% as.data.frame()
(load("ICBnetIS.RData"))
tmp <- as.data.frame(t(gsva(as.matrix(exprSet), ICBnetIS, method = "ssgsea")))
surv_df <- cbind(surv_expr %>% dplyr::select(1:3), tmp)

library(survival)
library(survminer)
cutoff <- survminer::surv_cutpoint(surv_df, time = "time", event = "status", variables = "ICBnetIS", minprop = 0.1) %>% .$cutpoint
cutoff <- cutoff$cutpoint
surv_df$group <- ifelse(surv_df$ICBnetIS >= cutoff, "High", "Low")
surv_df$time <- surv_df$time/ 30
sfit <- survfit(Surv(time, status) ~ group, data = surv_df)
diff <- survdiff(formula = Surv(time, status) ~ group, data = surv_df, rho = 0)
pval <- pchisq(diff$chisq, length(diff$n) - 1, lower.tail = FALSE)

# Function
library(clusterProfiler)
library(org.Hs.eg.db)
deg <- fread("DEG_ICBnetIS.txt")

KEGG <- read.csv("../HSA_KEGG.csv")
TERM2GENE <- KEGG[, c("KEGGID", "ENTREZID")]
TERM2NAME <- KEGG[, c("KEGGID", "DESCRIPTION")]

ID <- bitr(deg$symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
yy <- enricher(ID$ENTREZID, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 1)
head(yy)
yy1 <- setReadable(yy, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
head(yy1)
yy2 <- yy1[yy1$pvalue < 0.05, asis = T]

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
df <- fread("nrDEG_ICBnetIS.txt")
df <- df %>% dplyr::select(symbol, logFC)
mRNA_si <- df
head(mRNA_si)
names(mRNA_si)[1] <- "HUGO"
names(mRNA_si)[2] <- "Weight"
mRNA_si$Weight <- as.numeric(as.character(mRNA_si$Weight))
mRNA_si <- mRNA_si[order(mRNA_si$Weight, decreasing = T), ]
mRNA_si <- mRNA_si %>% dplyr::filter(!is.na(Weight))
si.id <- mRNA_si$Weight
names(si.id) <- mRNA_si$HUGO
head(si.id)

hallmark <- toTable(org.Hs.egGO) %>% dplyr::filter(Ontology == "BP")
names(hallmark)[1] <- "ENTREZID"
ID <- bitr(hallmark$ENTREZID, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
hallmark <- hallmark %>%
  inner_join(ID) %>%
  dplyr::select(-ENTREZID)

TERM2GENE <- hallmark[, c("go_id", "SYMBOL")]
TERM2NAME <- go2term(hallmark$go_id)
set.seed(2020)
clustergsea <- GSEA(
  geneList = si.id, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, verbose = F,
  minGSSize = 10, maxGSSize = 500, nPerm = 10000, pvalueCutoff = 1
)
clustergsea <- clustergsea[clustergsea$pvalue < 0.05, asis = T]

library(progeny)
pathways <- progeny(expr %>% as.matrix(), scale=TRUE,
                    organism="Human",
                    top = 100, perm = 1)

load("immunogram.RData")
library(GSVA)
gsym.expr <- exprSet
immunogram_gsva <- gsva(as.matrix(gsym.expr), Immunogram, method = "ssgsea")

load("cancer.cell.RData")
library(GSVA)
gsym.expr <- exprSet
cancer.cell_gsva <- gsva(as.matrix(gsym.expr), gs, method = "ssgsea")

# Immune
setwd("./IOBR-cell/")
l <- list.files()
l2 <- lapply(l, FUN_ann)
anno <- do.call(rbind, l2)
anno$Category[which(anno$Category == "estimate")] <- "ESTIMATE"
anno$Category[which(anno$Category == "ssGSEA")] <- "Pornpimol"
anno$Category[which(anno$Category == "timer")] <- "TIMER"
anno$Category[which(anno$Category == "mcp")] <- "MCPcounter"

FUN_cor <- function(i, temp) {
  res <- cor.test(temp[1, ] %>% as.numeric(), temp[i, ] %>% as.numeric())
  p <- res$p.value
  cor <- res$estimate
  df <- data.frame(p = p, cor = cor)
  rownames(df) <- rownames(temp)[i]
  return(df)
}

FUN_cor_df <- function(i) {
  temp <- read.csv(i) %>% column_to_rownames(names(.)[1])
  names(temp) <- sapply(names(temp), FUN_name) %>% as.character()
  temp <- temp %>% t() %>% as.data.frame()
  temp <- temp[, ICBnetIS$sample]
  temp <- rbind(ICBnetIS$ICBnetIS, temp)
  temp_l <- lapply(2:nrow(temp), FUN = FUN_cor, temp = temp)
  res <- do.call(rbind, temp_l) %>% rownames_to_column("cell")
  p.value <- res$p
  res$sig.label <- ifelse(p.value < 0.001,"****",
                          ifelse(p.value < 0.005,"***",
                                 ifelse(p.value < 0.01,"**",
                                        ifelse(p.value < 0.05,"*",""))))
  res$lab <- str_c(res$cell, res$sig.label)
  expr_df <- temp[-1, ]
  exprSet <- expr_df
  n <- t(scale(t(exprSet)))
  n[n > 1] <- 1
  n[n < -1] <- -1
  exprSet <- n %>% as.data.frame()
  rownames(exprSet) <- res$lab
  return(exprSet)
}
l3 <- lapply(l, FUN_cor_df)

x <- fread("immunomodulator.txt") %>% as.data.frame()
names(x)[1] <- "gene"
names(x)[2] <- "Super category"
x <- x %>% arrange(`Super category`)

inter <- intersect(rownames(expr), x$gene)
x <- x %>% dplyr::filter(gene %in% inter)
inter <- x$gene
expr <- expr[inter, ] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample")
anno <- anno[inter, ]
anno <- anno %>%
  dplyr::select(-2) %>%
  rownames_to_column("gene")
FUN <- function(gene) {
  temp <- df[, c(gene, "cluster")]
  names(temp)[1] <- "gene"
  pv <- wilcox.test(gene ~ cluster, data = temp)$p.value
  return(data.frame(gene = gene, pv = pv))
}
l <- lapply(names(df)[-1], FUN)
res <- do.call(rbind, l)
p.value <- res$pv
res$sig.label <- ifelse(p.value < 0.001, "****",
                        ifelse(p.value < 0.005, "***",
                               ifelse(p.value < 0.01, "**",
                                      ifelse(p.value < 0.05, "*", "")
                               )
                        )
)

# Immunotherapy
library(GSVA)
load("ICBnetIS.RData")
tmp <- gsva(as.matrix(exprSet), ICBnetIS, method = "ssgsea")
df <- tmp %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample")
df <- df %>% inner_join(pd)
df$Response <- factor(df$Response, levels = c("NR", "R"))
roc1 <- roc(df$Response, df$ICBnetIS,
            data = df, auc = TRUE,
            levels = c("NR", "R")
)

