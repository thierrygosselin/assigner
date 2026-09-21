# Rscript inst/experimental/prepare-ppr-equal-information-extension.R PILOT_DESIGN NEW_DESIGN
# Preserve the eight-panel pilot and add eight independent replicates per cell.
args <- commandArgs(TRUE);stopifnot(length(args)==2L,!file.exists(args[2]))
pilot <- readRDS(args[1])
panels <- expand.grid(N=c(120L,360L),L=c(1000L,10000L),REPLICATE=3:10)
panels$K <- 3L
panels$F <- sqrt(3/(panels$N*panels$L))
panels$SEED <- 922100L+seq_len(nrow(panels))
panels$ID <- with(panels,paste0("N",N,"_L",L,"_rep",REPLICATE))
combined <- rbind(pilot,panels)
stopifnot(nrow(panels)==32L,all(table(combined$N,combined$L)==10L),
 !anyDuplicated(combined$ID),!anyDuplicated(combined$SEED),
 all(abs(combined$N*combined$L*combined$F^2-3)<1e-12),
 all(combined$N %% combined$K==0))
saveRDS(panels,args[2]);print(panels,row.names=FALSE)
