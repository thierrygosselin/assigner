# Rscript inst/experimental/prepare-fixed-f002.R NEW_DESIGN_RDS
# Fixed generating F=0.02; measured WC84 differentiation is an evaluation output.
args <- commandArgs(TRUE);stopifnot(length(args)==1L,!file.exists(args[1]))
panels <- expand.grid(N=c(120L,360L),L=c(1000L,10000L),REPLICATE=1:10)
panels$K <- 3L;panels$F <- .02
panels$SEED <- 924100L+seq_len(nrow(panels))
panels$ID <- with(panels,paste0("F002_N",N,"_L",L,"_rep",REPLICATE))
stopifnot(nrow(panels)==40L,all(table(panels$N,panels$L)==10L),
 !anyDuplicated(panels$ID),!anyDuplicated(panels$SEED),
 all(panels$N %% panels$K==0),all(panels$F==.02))
saveRDS(panels,args[1]);print(panels,row.names=FALSE)
