# WGCNA
library(WGCNA)

sel <- c("TCRsignalingPathway", "BCRSignalingPathway", "Antigen_Processing_and_Presentation")
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
genes <- fread("pathway.txt") %>%
  dplyr::filter(Category == "TCRsignalingPathway") %>%
  .$Symbol
x <- RandomWalkRestart(Network, intersect(genes, rownames(Network)))
res <- x %>%
  as.data.frame() %>%
  set_colnames("rank") %>%
  rownames_to_column("gene") %>%
  arrange(desc(rank))
write.csv(res, "rank_TCR.csv", row.names = F)

genes <- fread("pathway.txt") %>%
  dplyr::filter(Category == "BCRSignalingPathway") %>%
  .$Symbol
x <- RandomWalkRestart(Network, intersect(genes, rownames(Network)))
res <- x %>%
  as.data.frame() %>%
  set_colnames("rank") %>%
  rownames_to_column("gene") %>%
  arrange(desc(rank))
write.csv(res, "rank_BCR.csv", row.names = F)

genes <- fread("pathway.txt") %>%
  dplyr::filter(Category == "Antigen_Processing_and_Presentation") %>%
  .$Symbol
x <- RandomWalkRestart(Network, intersect(genes, rownames(Network)))
res <- x %>%
  as.data.frame() %>%
  set_colnames("rank") %>%
  rownames_to_column("gene") %>%
  arrange(desc(rank))
write.csv(res, "rank_APP.csv", row.names = F)

library(clusterProfiler)
library(org.Hs.eg.db)
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

lapply(1:3, function(x) {
  set.seed(2020)
  clustergsea <- GSEA(
    geneList = all_glist[[x]], TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, verbose = F,
    minGSSize = 0, maxGSSize = 500, nPerm = 1000, pvalueCutoff = 1
  )
  clustergsea <- clustergsea[clustergsea$pvalue < 0.05, asis = T]
  return(clustergsea)
}) -> m_gsea_list

# ICBnetIS
library(survival)

surv_expr %<>% dplyr::select(c("time", "status", genes))
## month
surv_expr$time <- surv_expr$time / 30
FUN_surv <- function(i, surv_expr) {
  Sur <- Surv(surv_expr$time, surv_expr$status)
  cox <- coxph(Sur ~ surv_expr[, i], data = surv_expr)
  coxSummary <- summary(cox)
  gene <- names(surv_expr)[i]
  HR <- coxSummary$coefficients[, "exp(coef)"]
  lower <- coxSummary$conf.int[, 3]
  upper <- coxSummary$conf.int[, 4]
  PValue <- round(coxSummary$coefficients[, 5], 6)
  res <- data.frame(gene = gene, HR = HR, lower.95 = lower, upper.95 = upper, pvalue = PValue)
  return(res)
}
l <- pbapply::pblapply(3:ncol(surv_expr), FUN = FUN_surv, surv_expr = surv_expr)
Univar <- do.call(rbind, l)
Univar$gene <- Univar$gene %>% str_replace_all("_", "-")

head(Univar)
df <- Univar
names(df)[2] <- "HR"
names(df)[3] <- "lower"
names(df)[4] <- "upper"
names(df)[5] <- "pvalue"

df <- df %>% inner_join(imp2)
df$pvalue <- format(df$pvalue, digits = 3) %>% as.numeric()
df$pvalue <- ifelse(df$pvalue < 0.0001, "<0.0001", df$pvalue)
names(df)[1] <- "id"

df <- df %>% arrange(id)

column_info <- tribble(
  ~id, ~group, ~name, ~geom, ~palette, ~options,
  "id", NA, "gene", "text", NA, list(width = 6),
  "HR", "Unicox", "HR", "rect", "palette1", list(scale = FALSE),
  # "HR",   "Unicox",        "",           "text",        "palette1",  list(label = "HR"),
  "lower", "Unicox", "lower", "rect", "palette1", list(),
  "upper", "Unicox", "upper", "rect", "palette1", list(),
  "pvalue", "Unicox", "pvalue", "text", "palette1", list(width = 6),
  "imp", "RandomSurvivalForest", "importance", "bar", "palette2", list(width = 10)
)
column_groups <- tribble( # tribble_start
  ~Category, ~group, ~palette,
  "Unicox", "Unicox", "palette1",
  "RandomSurvivalForest", "RandomSurvivalForest", "palette2"
) # tribble_end

palettes <- tribble(
  ~palette, ~colours,
  "palette1", grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu") %>% rev())(50),
  "palette2", grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu") %>% rev())(50),
)

g <- funky_heatmap(df,
  column_info = column_info,
  column_groups = column_groups,
  palettes = palettes,
  expand = list(xmax = 4)
)

tmp <- as.data.frame(t(gsva(as.matrix(exprSet), ICBnetIS, method = "ssgsea")))
surv_df <- cbind(surv_expr %>% dplyr::select(1:3), tmp)

library(survival)
library(survminer)
cutoff <- survminer::surv_cutpoint(surv_df, time = "time", event = "status", variables = "ICBnetIS", minprop = 0.1) %>% .$cutpoint
cutoff <- cutoff$cutpoint
surv_df$group <- ifelse(surv_df$ICBnetIS >= cutoff, "High", "Low")
surv_df$time <- surv_df$time / 30
sfit <- survfit(Surv(time, status) ~ group, data = surv_df)
diff <- survdiff(formula = Surv(time, status) ~ group, data = surv_df, rho = 0)
pval <- pchisq(diff$chisq, length(diff$n) - 1, lower.tail = FALSE)

# Function
library(clusterProfiler)
library(org.Hs.eg.db)

KEGG <- read.csv("../HSA_KEGG.csv")
TERM2GENE <- KEGG[, c("KEGGID", "ENTREZID")]
TERM2NAME <- KEGG[, c("KEGGID", "DESCRIPTION")]

ID <- bitr(deg$symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
yy <- enricher(ID$ENTREZID, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 1)
head(yy)
yy1 <- setReadable(yy, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
head(yy1)
yy2 <- yy1[yy1$pvalue < 0.05, asis = T]

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

load("immunogram.RData")
library(GSVA)
gsym.expr <- exprSet
immunogram_gsva <- gsva(as.matrix(gsym.expr), Immunogram, method = "ssgsea")

load("cancer.cell.RData")
library(GSVA)
gsym.expr <- exprSet
cancer.cell_gsva <- gsva(as.matrix(gsym.expr), gs, method = "ssgsea")

# Immunotherapy

generateInputFileForSubMap <- function(in_gct, gct_file, cls_file, sam_info, type_name = "type") {
  in_gct <- data.frame(
    GeneID = rownames(in_gct),
    description = "na",
    in_gct,
    stringsAsFactors = F,
    check.names = F
  )
  cat("#1.2\n", file = gct_file)
  cat(nrow(in_gct), "\t", ncol(in_gct) - 2, "\n", file = gct_file, append = T)
  cat(paste(colnames(in_gct), collapse = "\t"), "\n", file = gct_file, append = T)
  for (i in 1:nrow(in_gct)) cat(paste(in_gct[i, ], collapse = "\t"), "\n", file = gct_file, append = T)

  cat(nrow(sam_info), length(levels(factor(sam_info$rank))), 1, "\n", file = cls_file)
  cat("#", paste0(levels(factor(sam_info[, type_name])), collapse = " "), "\n", file = cls_file, sep = "", append = T)
  cat(as.numeric(factor(sam_info[, type_name])), file = cls_file, append = T)
}

skcm.immunotherapy.logNC <- read.table("skcm.immunotherapy.47samples.log2CountsNorm.txt", sep = "\t", row.names = 1, header = T, check.names = F, stringsAsFactors = F)
rownames(skcm.immunotherapy.logNC) <- toupper(rownames(skcm.immunotherapy.logNC))
skcm.immunotherapy.info <- read.table("skcm.immunotherapy.47sampleInfo.txt", sep = "\t", row.names = 1, header = T, check.names = F, stringsAsFactors = F)

skcm.immunotherapy.info <- skcm.immunotherapy.info[order(skcm.immunotherapy.info$label), ]
skcm.immunotherapy.info$rank <- rep(c(1, 2, 3, 4), times = as.character(table(skcm.immunotherapy.info$label)))

tmp <- TCGA %>%
  column_to_rownames("sample") %>%
  dplyr::select(-c(1:2)) %>%
  t() %>%
  as.data.frame()
tmp <- tmp[, rownames(ann)]
identical(names(tmp), rownames(ann))
GENELIST <- intersect(rownames(tmp), rownames(skcm.immunotherapy.logNC))

sam_info <- skcm.immunotherapy.info
in_gct <- skcm.immunotherapy.logNC[GENELIST, rownames(skcm.immunotherapy.info)]

generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

samples.C1 <- rownames(ann %>% dplyr::filter(ImmClust == "C1"))
samples.C2 <- rownames(ann %>% dplyr::filter(ImmClust == "C2"))

sam_info <- data.frame("ImmClust" = c(samples.C1, samples.C2), row.names = c(samples.C1, samples.C2))
sam_info$rank <- rep(c(1, 2), times = c(length(samples.C1), length(samples.C2))) #

in_gct <- tmp[GENELIST, rownames(sam_info)] 
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")


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
