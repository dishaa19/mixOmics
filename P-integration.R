library(mixOmics)
data(stemcells)
data = stemcells$gene
type.id = stemcells$celltype
exp = stemcells$study
res = mint.splsda(X=data,Y=type.id,ncomp=3,keepX=c(10,5,15),study=exp)
out = tune.mint.splsda(X=data,Y=type.id,ncomp=2,near.zero.var=FALSE,
                       study=exp,test.keepX=seq(1,10,1))
plot(out)
out$choice.ncomp
out$choice.keepX

## Not run: 

out = tune.mint.splsda(X=data,Y=type.id,ncomp=2,near.zero.var=FALSE,
                       study=exp,test.keepX=seq(1,10,1))

out$choice.keepX


## only tune component 2 and keeping 10 genes on comp1
out = tune.mint.splsda(X=data,Y=type.id,ncomp=2, study=exp,
                       already.tested.X = c(10),
                       test.keepX=seq(1,10,1))
out$choice.keepX
plot(out)

## End(Not run)
