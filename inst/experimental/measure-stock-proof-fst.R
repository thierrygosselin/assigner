args <- commandArgs(TRUE); stopifnot(length(args)==1L)
root <- normalizePath(args[1],mustWork=TRUE)
design <- readRDS(file.path(root,"design-and-audit.rds"))$design
calculate <- function(j) {
  row <- design[j,]
  value <- if(row$TRUE_K==1) NA_real_ else {
    d <- readRDS(file.path(root,paste(row$CONDITION,row$REPLICATE,sep="_"),"truth.rds"))
    g <- d$train[,seq_len(row$MARKERS),drop=FALSE]
    codes <- matrix(c(11L,12L,22L)[g+1L],nrow(g),ncol(g))
    hierfstat::wc(data.frame(STRATA=d$group,codes,check.names=FALSE))$FST
  }
  data.frame(row,MEASURED_WC84=value)
}
result <- do.call(rbind,lapply(seq_len(nrow(design)),calculate))
saveRDS(result,file.path(root,"measured-WC84.rds"));print(result,row.names=FALSE)
