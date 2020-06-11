## install BiocManager if not installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
## install mixOmics
BiocManager::install('mixOmics')

## install devtools if not installed
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
## install mixOmics
devtools::install_github("mixOmicsTeam/mixOmics")

library(mixOmics) 
#data(nutrimouse)
#data("breast.tumors")
#X <- nutrimouse$gene
#MyResult.pca <- pca(X)
#plotIndiv(MyResult.pca)
#plotVar(MyResult.pca)


data(liver.toxicity)
X <- liver.toxicity$gene
MyResult.pca <- pca(X)     # 1 Run the method
plotIndiv(MyResult.pca)    # 2 Plot the samples
plotVar(MyResult.pca, cutoff = 0.8)
plotIndiv(MyResult.pca, group = liver.toxicity$treatment$Dose.Group, 
          legend = TRUE)
plotIndiv(MyResult.pca, ind.names = FALSE,
          group = liver.toxicity$treatment$Dose.Group,
          pch = as.factor(liver.toxicity$treatment$Time.Group),
          legend = TRUE, title = 'Liver toxicity: genes, PCA comp 1 - 2',
          legend.title = 'Dose', legend.title.pch = 'Exposure')
plotIndiv(MyResult.pca, group = liver.toxicity$treatment$Dose.Group, 
          legend = TRUE)
MyResult.pca2 <- pca(X, ncomp = 3)
plotIndiv(MyResult.pca2, comp = c(1,3), legend = TRUE,
          group = liver.toxicity$treatment$Time.Group,
          title = 'Multidrug transporter, PCA comp 1 - 3')
MyResult.pca2 <- pca(X, ncomp = 3)
plotIndiv(MyResult.pca2, comp = c(1,3), legend = TRUE,
          group = liver.toxicity$treatment$Time.Group,
          title = 'Multidrug transporter, PCA comp 1 - 3')
plot(MyResult.pca2)

