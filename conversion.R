library(mixOmics)
load("stemcells.rda")
X <- stemcells$gene
write.csv(X, file = "stemcells.csv")
