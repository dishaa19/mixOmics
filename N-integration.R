library(mixOmics)
data(breast.TCGA)
# extract training data and name each data frame
X <- list(mRNA = breast.TCGA$data.train$mrna, 
          miRNA = breast.TCGA$data.train$mirna, 
          protein = breast.TCGA$data.train$protein)
Y <- breast.TCGA$data.train$subtype
summary(Y)
list.keepX <- list(mRNA = c(16, 17), miRNA = c(18,5), protein = c(5, 5))
MyResult.diablo <- block.splsda(X, Y, keepX=list.keepX)
plotIndiv(MyResult.diablo) ## sample plot
plotVar(MyResult.diablo) ## variable plot
MyResult.diablo2 <- block.plsda(X, Y)
plotIndiv(MyResult.diablo, 
          ind.names = FALSE, 
          legend=TRUE, cex=c(1,2,3),
          title = 'BRCA with DIABLO')
plotVar(MyResult.diablo, var.names = c(FALSE, FALSE, TRUE),
        legend=TRUE, pch=c(16,16,1))
plotDiablo(MyResult.diablo, ncomp = 1)
circosPlot(MyResult.diablo, cutoff=0.7)
# minimal example with margins improved:
# cimDiablo(MyResult.diablo, margin=c(8,20))
# extended example:
cimDiablo(MyResult.diablo, color.blocks = c('darkorchid', 'brown1', 'lightgreen'), comp = 1, margin=c(8,20), legend.position = "right")
#plotLoadings(MyResult.diablo, contrib = "max")
plotLoadings(MyResult.diablo, comp = 2, contrib = "max")
network(MyResult.diablo, blocks = c(1,2,3),
        color.node = c('darkorchid', 'brown1', 'lightgreen'), 
        cutoff = 0.6, save = 'jpeg', name.save = 'DIABLOnetwork')
set.seed(123) # for reproducibility in this vignette
MyPerf.diablo <- perf(MyResult.diablo, validation = 'Mfold', folds = 5, 
                      nrepeat = 10, 
                      dist = 'centroids.dist')

#MyPerf.diablo  # lists the different outputs

# Performance with Majority vote
#MyPerf.diablo$MajorityVote.error.rate
Myauc.diablo <- auroc(MyResult.diablo, roc.block = "miRNA", roc.comp = 2)
# prepare test set data: here one block (proteins) is missing
X.test <- list(mRNA = breast.TCGA$data.test$mrna, 
               miRNA = breast.TCGA$data.test$mirna)

Mypredict.diablo <- predict(MyResult.diablo, newdata = X.test)
# the warning message will inform us that one block is missing
#Mypredict.diablo # list the different outputs
confusion.mat <- get.confusion_matrix(
  truth = breast.TCGA$data.test$subtype, 
  predicted = Mypredict.diablo$MajorityVote$centroids.dist[,2])
kable(confusion.mat)
get.BER(confusion.mat)
MyResult.diablo$design
MyDesign <- matrix(c(0, 0.1, 0.3,
                     0.1, 0, 0.9,
                     0.3, 0.9, 0),
                   byrow=TRUE,
                   ncol = length(X), nrow = length(X),
                   dimnames = list(names(X), names(X)))
MyDesign
MyResult.diablo.design <- block.splsda(X, Y, keepX=list.keepX, design=MyDesign)
