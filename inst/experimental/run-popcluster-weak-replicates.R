# Rscript inst/experimental/run-popcluster-weak-replicates.R NEW_OUTPUT
# Run from assigner. Fixed design: eight new independent weak panels, K supplied.
args <- commandArgs(TRUE); stopifnot(length(args) %in% 1:2)
repo <- normalizePath(".")
sources <- c("../grur/R/simulation_benchmark.R",
 "inst/experimental/benchmark-stockr-paired.R",
 "inst/experimental/diagnose-truth-starts.R",
 "inst/experimental/prepare-ppr-popcluster.R")
for(p in sources) source(p)
output <- args[1]
if(file.exists(output)) stop("Use a new output directory; existing results will not be overwritten.")
dir.create(output,recursive=TRUE); output <- normalizePath(output)
binary <- "/Users/gos040/Downloads/PopCluster/Bin/PopClusterMac"
design <- list(replicates=3:10, N=300L, L=10000L, K=3L, F=.001,
 starts=3L, workers=2L, truth.seeds=911000+3:10,
 note="Eight fresh complete-genotype panels; no K selection or held-out prediction.")
panels <- data.frame(ID=paste0("replicate_",3:10),REPLICATE=3:10,
 N=300L,L=10000L,K=3L,F=.001,SEED=911000+3:10)
if(length(args)==2L) {
 panels <- readRDS(args[2])
 stopifnot(is.data.frame(panels),nrow(panels)>0,
   all(c("ID","REPLICATE","N","L","K","F","SEED") %in% names(panels)),
   !anyNA(panels),!anyDuplicated(panels$ID),!anyDuplicated(panels$SEED),
   all(grepl("^[A-Za-z0-9_-]+$",panels$ID)),
   all(panels$N>=6),all(panels$K>=2),all(panels$N %% panels$K==0),
   all(panels$L>=2),all(panels$F>=0 & panels$F<1))
 design <- list(panels=panels,starts=3L,workers=2L,
   note="Matched complete-genotype panels; supplied K only, not K selection or prediction.")
}
saveRDS(list(design=design,session=sessionInfo(),stockR=as.character(packageVersion("stockR")),
 binary.hash=tools::md5sum(binary),source.hash=tools::md5sum(c(sources,
 "inst/experimental/run-popcluster-weak-replicates.R")),
 source.text=lapply(c(sources,"inst/experimental/run-popcluster-weak-replicates.R"),readLines)),
 file.path(output,"design-locked.rds"))
run <- function(j) {
 panel <- panels[j,]; r <- panel$REPLICATE
 dest <- file.path(output,panel$ID); dir.create(dest)
 setwd(dest); on.exit(setwd(repo))
 cat("Starting panel",panel$ID,"\n")
 truth <- .grur_structure_truth(rep(as.integer(panel$N/panel$K),panel$K),
   panel$L,panel$F,panel$SEED)
 saveRDS(truth,"truth.rds")
 keep <- which(apply(truth$genotype,2,function(x)length(unique(x))>1L))
 g <- truth$genotype[,keep,drop=FALSE]
 ids <- sprintf("sample%04d",seq_len(nrow(g)))
 encoded <- ppr_encode(g)
 stopifnot(identical(unname(ppr_decode(encoded)),unname(g)))
 writeLines(as.vector(rbind(ids,encoded)),"train.dat")
 settings <- c(nrow(g),ncol(g),1,0,if(length(args)==1L) 940000+r else panel$SEED+100000L,"train.dat","mixture",
 0,panel$K,panel$K,3,0,1,0,0,1,0,0,0,1,0,0,0,2)
 stopifnot(length(settings)==24L)
 writeLines(as.character(settings),"mixture.PcPjt")
 fits <- list(); diagnostics <- list()
 for(s in 1:3) {
   seed <- if(length(args)==1L) 820000+r*10+s else panel$SEED*10+s
   set.seed(seed)
   warnings <- character(); before <- proc.time()[3]
   fit <- withCallingHandlers(suppressMessages(stockR::stockSTRUCTURE(t(g),K=panel$K,
     sample.grps=factor(seq_len(nrow(g)),levels=seq_len(nrow(g))),
     start.grps=NULL,control=list(quiet=TRUE))),warning=function(w) {
       warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning")
     })
   it <- length(fit$margLogLTrace)
   if(it>fit$control$EM.maxit) it <- it-1L
   fits[[s]] <- fit
   diagnostics[[s]] <- data.frame(START=s,SEED=seed,
     ITERATIONS=it,CONVERGED=it<=fit$control$EM.maxit,
     LOGLIK=fit$margLogL,SECONDS=unname(proc.time()[3]-before),WARNING=paste(warnings,collapse=" | "))
   saveRDS(list(fits=fits,diagnostics=diagnostics),"stock-checkpoint.rds")
 }
 tab <- do.call(rbind,diagnostics)
 eligible <- which(tab$CONVERGED & is.finite(tab$LOGLIK))
 if(!length(eligible)) stop("No eligible stockR fit")
 selected <- eligible[which.max(tab$LOGLIK[eligible])]
 saveRDS(list(fits=fits,diagnostics=tab,selected=selected),"stock.rds")
 started <- Sys.time()
 status <- system2("/usr/bin/arch",c("-x86_64",shQuote(binary),"INP:mixture.PcPjt"),
                   stdout="console.log",stderr="stderr.log")
 elapsed <- as.numeric(difftime(Sys.time(),started,units="secs"))
 saveRDS(list(status=status,seconds=elapsed),"execution.rds")
 if(status!=0L) stop("PopCluster failed")
 pop <- lapply(1:3,function(s) {
   lines <- readLines(paste0("mixture_K_",panel$K,"_R_",s))
   for(setting in c("Model = 1","Scaling = 0","AllSNP = 1","PopDataQ = 0","PopFlagQ = 0","DataForm = 1"))
     stopifnot(any(trimws(lines)==setting))
   h <- grep("^\\s*Label.*ML_C",lines); stopifnot(length(h)==1)
   z <- read.table(text=paste(lines[h+seq_len(nrow(g))],collapse="\n"),
                   col.names=c("ID","HOM","HET","MIS","F","CLUSTER"))
   stopifnot(identical(z$ID,ids),all(z$MIS==0),all(z$HOM+z$HET==ncol(g)),
             all(z$HET==rowSums(g==1L)))
   ll <- grep("LogL, NumIter",lines,value=TRUE); stopifnot(length(ll)==1)
   list(labels=z$CLUSTER,loglik=as.numeric(strsplit(trimws(sub(".*=","",ll)),",")[[1]][1]))
 })
 best <- which.max(vapply(pop,`[[`,numeric(1),"loglik"))
 # Selection is frozen before scoring against simulation truth.
 saveRDS(list(stock=selected,pop=best),"choices-locked.rds")
 pred <- .bp_stock_predict(fits[[selected]],g)
 stopifnot(isTRUE(all.equal(sum(pred$sample.loglik),fits[[selected]]$margLogL,tolerance=1e-8)))
 labels <- max.col(pred$membership,ties.method="first")
 stopifnot(all(rowSums(pred$membership==apply(pred$membership,1,max))==1L))
 codes <- matrix(c(11L,12L,22L)[g+1L],nrow(g),ncol(g))
 measured <- hierfstat::wc(data.frame(STRATA=truth$group,codes,check.names=FALSE))$FST
 result <- data.frame(ID=panel$ID,REPLICATE=r,K=panel$K,N=nrow(g),N_PER_GROUP=panel$N/panel$K,
   L_NOMINAL=panel$L,L=ncol(g),F_GENERATING=panel$F,
   F_MEASURED=measured,S_GENERATING=nrow(g)*ncol(g)*panel$F^2,
   S_MEASURED=nrow(g)*ncol(g)*measured^2,
   POP_ARI=.tb_ari(pop[[best]]$labels,truth$group),STOCK_ARI=.tb_ari(labels,truth$group),
   POP_START=best,STOCK_START=selected,POP_SECONDS=elapsed,
   STOCK_FAILED_STARTS=sum(!tab$CONVERGED),STOCK_WARNINGS=sum(nzchar(tab$WARNING)))
 saveRDS(list(result=result,pop=pop,keep=keep),"result.rds")
 cat("Completed replicate",r,"PopCluster",result$POP_ARI,"stockR",result$STOCK_ARI,"\n")
 result
}
results <- parallel::mclapply(seq_len(nrow(panels)),run,mc.cores=2L,mc.preschedule=FALSE,
                             mc.set.seed=FALSE)
saveRDS(results,file.path(output,"results.rds"))
if(any(vapply(results,inherits,logical(1),"try-error"))) stop("Some replicates failed; inspect checkpoints.")
print(do.call(rbind,results),row.names=FALSE)
