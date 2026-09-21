# Experimental accelerator ONLY for complete near-certain GL simulations.
# Same penalized EM objective/updates as .sm_block_fit in the hard-genotype limit.
# Never use for uncertain GLs, missing calls or real-data likelihood inference.
.hl_decode <- function(b) {
  x <- b$log.likelihood; n <- dim(x)[1]; m <- dim(x)[3]
  a <- matrix(x[,1,],n,m); c <- matrix(x[,2,],n,m); e <- matrix(x[,3,],n,m)
  if(anyNA(x)||any(!b$informative)||any((a==0)+(c==0)+(e==0)!=1)||
     any(abs(x[x!=0]+80*log(10))>1e-5)) stop("Not a complete near-certain benchmark GDS.")
  (c==0)+2*(e==0)
}
.hl_read <- function(path) {
  ctx <- .bs_context(path); on.exit(ctx$close())
  if(ctx$n*ctx$m>5e6) stop("Hard-limit benchmark memory guard exceeded.")
  g <- matrix(NA_real_,ctx$n,ctx$m,dimnames=list(ctx$ids,ctx$variants))
  for(index in ctx$blocks) {
    g[,index] <- .hl_decode(ctx$reader(index))
  }
  g
}
.hl_predict <- function(g,p,pi) {
  # Binomial combinatorial term must be retained for comparable likelihoods.
  joint <- g %*% (log(p)-log1p(-p))
  joint <- sweep(joint,2,2*colSums(log1p(-p))+log(pi),"+")
  joint <- joint+rowSums(g==1)*log(2)
  peak <- apply(joint,1,max); score <- peak+log(rowSums(exp(joint-peak)))
  list(membership=exp(joint-score),sample.loglik=score,loglik=sum(score))
}
.hl_fit <- function(g,p,pi,strength,max.iter=300L,tolerance=1e-8,epsilon=1e-8) {
  stopifnot(!anyNA(g),all(g%in%0:2),strength>=0)
  center <- colMeans(g)/2
  penalty <- function(p) strength*sum(center*log(p)+(1-center)*log1p(-p))
  state <- .hl_predict(g,p,pi); objective <- state$loglik+penalty(p)
  history <- objective; converged <- FALSE
  for(iter in seq_len(max.iter)) {
    tau <- state$membership; counts <- colSums(tau); pi <- counts/nrow(g)
    if(any(pi<=0)) stop("Empty component.")
    p <- sweep(crossprod(g,tau)+strength*center,2,2*counts+strength,"/")
    p[] <- pmin(1-epsilon,pmax(epsilon,p))
    next.state <- .hl_predict(g,p,pi); next.objective <- next.state$loglik+penalty(p)
    change <- next.objective-objective
    if(change < -1e-7*(1+abs(objective))) stop("EM objective decreased.")
    state <- next.state; objective <- next.objective; history <- c(history,objective)
    if(abs(change)<=tolerance*(1+abs(objective))) {converged<-TRUE;break}
  }
  c(list(frequency=p,proportion=pi,objective=objective,objective.history=history,
    converged=converged,iterations=iter,prior.strength=strength),state)
}
.hl_models <- function(g,k,starts,strength) {
  best <- NULL; diagnostics <- list()
  for(s in seq_along(starts)) {
    tau <- if(k==1) matrix(1,nrow(g),1) else stockR:::get.init.taus(starts[[s]],2/(k+1),k)
    p <- sweep(crossprod(g,tau),2,2*colSums(tau),"/")
    p[] <- pmin(1-1e-8,pmax(1e-8,p))
    before <- proc.time()[3]
    fit <- tryCatch(.hl_fit(g,p,colMeans(tau),strength),error=function(e)e)
    ok <- !inherits(fit,"error") && isTRUE(fit$converged)
    diagnostics[[s]] <- data.frame(START=s,CONVERGED=ok,
      SECONDS=unname(proc.time()[3]-before),
      ITERATIONS=if(inherits(fit,"error")) NA else fit$iterations,
      OBJECTIVE=if(inherits(fit,"error")) NA else fit$objective,
      ERROR=if(inherits(fit,"error")) conditionMessage(fit) else "")
    if(ok&&(is.null(best)||fit$objective>best$objective)) best <- fit
  }
  list(fit=best,diagnostics=do.call(rbind,diagnostics))
}
