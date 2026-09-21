# Rscript inst/experimental/prepare-ppr-equal-information.R NEW_DESIGN_RDS
# Eight independent balanced panels, supplied K=3 and generating PPR S=3.
args <- commandArgs(TRUE); stopifnot(length(args)==1L,!file.exists(args[1]))
panels <- expand.grid(N=c(120L,360L),L=c(1000L,10000L),REPLICATE=1:2)
panels$K <- 3L
panels$F <- sqrt(3/(panels$N*panels$L))
panels$SEED <- 921100L+seq_len(nrow(panels))
panels$ID <- with(panels,paste0("N",N,"_L",L,"_rep",REPLICATE))
stopifnot(all(abs(panels$N*panels$L*panels$F^2-3)<1e-12),
 !anyDuplicated(panels$SEED),all(panels$N %% panels$K==0))
saveRDS(panels,args[1]);print(panels,row.names=FALSE)
