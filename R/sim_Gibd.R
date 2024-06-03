#' Genomic IBD Matrix in AlphaSimR
#'
#' Genomic IBD Matrix in AlphaSimR
#'
#' Details
#' 
#' @param pop AlphaSimR Pop-class
#' @param SP simulation parameters for AlphaSimR
#' @param n.mark number of markers to sample (evenly spaced in Morgans)
#' @param n.core number of cores
#' 
#' @return genomic IBD matrix
#'
#' @import AlphaSimR
#' @importFrom utils write.csv
#' @import parallel
#' @export

sim_Gibd <- function(pop, SP, n.mark, n.core) {

  Gibd <- function(X,haps,ploidy) {
    m <- ncol(X)
    n <- nrow(X)/ploidy
    
    Wlist <- lapply(as.list(1:m),function(i){
      z <- split(X[,i],f=rep(1:n,each=ploidy))
      sapply(z,FUN=function(z){table(factor(z,levels=haps))})
    }) 
    K <- matrix(0,nrow=n,ncol=n)
    for (i in 1:m) {
      K <- K + crossprod(Wlist[[i]])
    }
    return(K/m/ploidy)
  }  
  
  stopifnot(inherits(pop,"Pop"))
  ploidy <- pop@ploidy
  
  x <- lapply(SP$genMap,function(map){
    u <- split(names(map),cut(map,breaks=n.mark))
    sapply(u,"[[",1)
  })

  chrom <- as.list(1:length(x))
  y <- mapply(FUN=function(chr,x){
    list(pullIbdHaplo(pop,chr=chr,simParam = SP)[,x])
  },chr=chrom,x=x)
  
  haps <- unique(unlist(lapply(y,function(v){unique(as.integer(v))})))
  
  if (n.core > 1) {
    cl <- makeCluster(n.core)
    clusterExport(cl=cl,varlist="Gibd",envir=environment())
    Glist <- parLapply(cl,X=y,fun=Gibd,haps=haps,ploidy=ploidy)
    stopCluster(cl)
  } else {
    Glist <- lapply(y,FUN=Gibd,haps=haps,ploidy=ploidy)
  }
  
  nchr <- length(chrom)
  G <- Glist[[1]]
  for (i in 2:nchr) {
    G <- G + Glist[[i]]
  }
  G <- G/nchr
  dimnames(G) <- list(pop@id,pop@id)
  return(G)
}
