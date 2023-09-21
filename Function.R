library(GSVA)
load("ICBnetIS.RData")
gs <- list(ICBnetIS = ICBnetIS)
tmp <- gsva(as.matrix(exprSet), gs, method = "ssgsea", ssgsea.norm = F)
df <- tmp %>%
  t() %>%
  as.data.frame()
