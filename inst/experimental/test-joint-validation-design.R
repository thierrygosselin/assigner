# Read-only checks of the actual selector expression and generated nested design.
args<-commandArgs(TRUE)
if(length(args)!=1L) stop("Supply the generated experiment directory.")
expressions<-parse("inst/experimental/run-joint-validation.R")
selector<-Filter(function(x) is.call(x)&&identical(x[[1]],as.name("<-"))&&
  identical(x[[2]],as.name("choose.model")),as.list(expressions))
stopifnot(length(selector)==1L)
environment<-new.env(); eval(selector[[1]],environment)
choose<-environment$choose.model
x<-data.frame(K=c(3,1,1),STRENGTH=c(2,100,2),COMPLETE=TRUE,SCORE=c(0,0,0))
s<-choose(x); stopifnot(s$K==1,s$STRENGTH==2)
x$SCORE[1]<-1; stopifnot(choose(x)$K==3)
x$COMPLETE[1]<-FALSE; stopifnot(is.null(choose(x)))
x$COMPLETE[1]<-TRUE; x$SCORE[1]<-NA_real_; stopifnot(is.null(choose(x)))

root<-args[1]; manifest<-readRDS(file.path(root,"manifest.rds"))
stopifnot(nrow(manifest)==4L)
for(condition in unique(manifest$CONDITION)) {
  small<-readRDS(file.path(root,"truth",paste0(condition,"_small.rds")))
  large<-readRDS(file.path(root,"truth",paste0(condition,"_joint_large.rds")))
  stopifnot(!length(intersect(rownames(large$train$genotype),rownames(large$test$genotype))),
    identical(small$train$genotype,large$train$genotype[rownames(small$train$genotype),1:100,drop=FALSE]),
    identical(small$test$genotype,large$test$genotype[,1:100,drop=FALSE]),
    identical(small$train$frequency,large$train$frequency[1:100,,drop=FALSE]))
}
cat("Selector ties/incomplete grids, nested panels, and disjoint evaluation-cohort checks passed\n")
