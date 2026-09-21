# Rscript .../measure-marker-benchmark-fst.R RUN
args<-commandArgs(TRUE); stopifnot(length(args)==1)
input<-normalizePath(args[1]); path<-file.path(input,"measured-WC84.rds")
if(file.exists(path)) stop("Measured results already exist.")
manifest<-readRDS(file.path(input,"manifest.rds"))
calculate<-function(j) {
  row<-manifest[j,]
  if(row$TRUE_K==1) return(data.frame(ID=row$ID,WC84=NA_real_))
  d<-readRDS(file.path(input,paste(row$CONDITION,row$REPLICATE,sep="_"),"truth.rds"))
  g<-d$train[,seq_len(row$MARKERS),drop=FALSE]
  codes<-matrix(c(11L,12L,22L)[g+1L],nrow(g),ncol(g))
  value<-hierfstat::wc(data.frame(STRATA=d$group,codes,check.names=FALSE))$FST
  data.frame(ID=row$ID,WC84=value)
}
results<-parallel::mclapply(seq_len(nrow(manifest)),calculate,mc.cores=2,mc.set.seed=FALSE)
if(any(vapply(results,inherits,logical(1),"try-error"))) stop("FST calculation failed.")
results<-do.call(rbind,results); saveRDS(results,path); print(results)
