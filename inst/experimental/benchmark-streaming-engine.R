# Experimental bounded-block reader for the vectorized likelihood engine.
# Owns a read-only GDS handle; callers MUST arrange on.exit(context$close()).
# No full sample-by-marker likelihood or informative matrix is retained.
.bs_prepare <- function(x, encoding=c("GL","PL","log")) {
  encoding<-match.arg(encoding)
  if(!is.numeric(x)||length(dim(x))!=3||dim(x)[2]!=3||any(dim(x)==0)||any(is.nan(x)))
    stop("Invalid likelihood array.")
  if(encoding=="PL"&&any(x<0|x!=floor(x),na.rm=TRUE))
    stop("PL must contain non-negative whole numbers.")
  if(encoding=="PL") x<--x*log(10)/10
  if(encoding=="GL") x<-x*log(10)
  n<-dim(x)[1]; m<-dim(x)[3]
  a<-matrix(x[,1,],n,m); b<-matrix(x[,2,],n,m); c<-matrix(x[,3,],n,m)
  missing<-is.na(a)&is.na(b)&is.na(c)
  if(any((is.na(a)|is.na(b)|is.na(c))&!missing)||any(x==Inf,na.rm=TRUE)||
     any(a==-Inf&b==-Inf&c==-Inf,na.rm=TRUE)) stop("Malformed or partially missing likelihood vector.")
  a[missing]<-b[missing]<-c[missing]<-0
  peak<-pmax(a,b,c); a<-a-peak; b<-b-peak; c<-c-peak
  x[,1,]<-a; x[,2,]<-b; x[,3,]<-c
  list(log.likelihood=x,informative=(a!=0|b!=0|c!=0))
}

.bs_context <- function(path, block.snps=250L, max.cells=100000L,
                        sample.ids=NULL, variant.ids=NULL) {
  stopifnot(length(block.snps)==1L, is.finite(block.snps), block.snps>=1,
    block.snps==floor(block.snps), length(max.cells)==1L,
    is.finite(max.cells), max.cells>=1)
  g<-SeqArray::seqOpen(path,readonly=TRUE)
  success<-FALSE
  on.exit(if(!success) SeqArray::seqClose(g))
  if(!is.null(sample.ids)) SeqArray::seqSetFilter(g,sample.id=sample.ids,verbose=FALSE)
  if(!is.null(variant.ids)) SeqArray::seqSetFilter(g,variant.id=variant.ids,verbose=FALSE)
  ids<-as.character(SeqArray::seqGetData(g,"sample.id"))
  variants<-SeqArray::seqGetData(g,"variant.id")
  n<-length(ids); m<-length(variants)
  if(!n||!m||anyDuplicated(ids)||anyDuplicated(variants)) stop("Invalid active IDs.")
  alleles<-SeqArray::seqGetData(g,"allele")
  if(anyNA(alleles)||any(lengths(strsplit(alleles,",",fixed=TRUE))!=2))
    stop("Only biallelic variants are supported.")
  size<-min(block.snps,floor(max.cells/n))
  if(size<1) stop("One marker exceeds the block memory guard.")
  blocks<-split(seq_len(m),ceiling(seq_len(m)/size))
  state<-new.env(parent=emptyenv()); state$closed<-FALSE
  state$max.read.cells<-0; state$reads<-0L
  close<-function() {
    if(!state$closed) {SeqArray::seqClose(g); state$closed<-TRUE}
    invisible(NULL)
  }
  reader<-function(index) {
    if(state$closed) stop("GDS context is closed.")
    if(!length(index)||anyNA(index)||any(index!=floor(index))||
       any(index<1|index>m)||anyDuplicated(index)||is.unsorted(index,strictly=TRUE))
      stop("Use increasing unique marker indices.")
    if(length(index)*n>max.cells) stop("Read exceeds the block memory guard.")
    SeqArray::seqSetFilter(g,variant.id=variants[index],verbose=FALSE)
    raw<-SeqArray::seqGetData(g,"annotation/format/GL")
    if(!is.list(raw)||length(raw$length)!=length(index)||anyNA(raw$length)||
       any(raw$length!=3L)||!identical(dim(raw$data),c(as.integer(n),as.integer(3*length(index)))))
      stop("Expected three diploid GL values per sample and marker.")
    x<-array(as.numeric(raw$data),c(n,3,length(index))); x[is.nan(x)]<-NA_real_
    state$max.read.cells<-max(state$max.read.cells,n*length(index))
    state$reads<-state$reads+1L
    .bs_prepare(x,"GL")
  }
  information<-numeric(n); marker.information<-numeric(m)
  for(index in blocks) {
    b<-reader(index)
    information<-information+rowSums(b$informative)
    marker.information[index]<-colSums(b$informative)
  }
  if(any(marker.information==0)) stop("All-uninformative marker.")
  # Do not retain the last preflight block in the reader's closure environment.
  rm(b,alleles)
  success<-TRUE
  list(ids=ids,variants=as.character(variants),n=n,m=m,blocks=blocks,
    reader=reader,information=information,close=close,audit=state)
}
