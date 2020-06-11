#mixOmics on SRBCT
library(mixOmics)
data("srbct")
X <- srbct$gene

MyResult.pca <- pca(X)     # 1 Run the method
plotIndiv(MyResult.pca)    # 2 Plot the samples
plotVar(MyResult.pca, cutoff = 0.8)      # 3 Plot the variables
plotIndiv(MyResult.pca, ind.names = TRUE,legend = FALSE, title = 'PCA on SRBCT')

#pls-da
MyResult.plsda <- plsda(X,Y) # 1 Run the method
plotIndiv(MyResult.plsda, title = 'PLSDA on SRBCT')    # 2 Plot the samples
plotVar(MyResult.plsda, cutoff = 0.7) 


#spls-da
X <- srbct$gene
Y <- srbct$class 

MyResult.splsda <- splsda(X, Y, keepX = c(50,50)) # 1 Run the method
plotIndiv(MyResult.splsda)                          # 2 Plot the samples (coloured by classes automatically)
plotIndiv(MyResult.splsda, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, star = TRUE, title = 'sPLS-DA on SRBCT',
          X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2')

plotVar(MyResult.splsda, var.names=FALSE)

background <- background.predict(MyResult.splsda, comp.predicted=2,
                                 dist = "max.dist") 
plotIndiv(MyResult.splsda, comp = 1:2, group = srbct$class,
          ind.names = FALSE, title = "Maximum distance",
          legend = TRUE,  background = background)

auc.plsda <- auroc(MyResult.splsda)
MyResult.splsda2 <- splsda(X,Y, ncomp=3, keepX=c(15,10,5))
selectVar(MyResult.splsda2, comp=1)$value
auc.plsda <- auroc(MyResult.splsda2)
plotLoadings(MyResult.splsda2, contrib = 'max', method = 'mean')
install.packages("rgl")
plotIndiv(MyResult.splsda2, style="3d")
MyResult.plsda2 <- plsda(X,Y, ncomp=10)
set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 3, 
                     progressBar = FALSE, nrepeat = 50) # we suggest nrepeat = 50

plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

