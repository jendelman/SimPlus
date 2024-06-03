#' Truncation selection of parents
#'
#' Truncation selection of parents
#'
#' Details
#' 
#' @param pop parental candidates, AlphaSimR Pop-class
#' @param SP simulation parameters for AlphaSimR
#' @param n.progeny number of total progeny
#' @param n.parent number of parents
#' @param geno.file geno filename with marker effects
#' @param weighted TRUE/FALSE for weighting marker effects
#' 
#' @return list containing
#' \describe{
#' \item{parents}{Pop-class of selected parents}
#' \item{matings}{data frame for sim_mate}
#' }
#'
#' @import AlphaSimR
#' @importFrom utils read.csv
#' @export

sim_trunc <- function(pop, SP, n.progeny, n.parent, geno.file, weighted) {

  stopifnot(inherits(pop,"Pop"))
  stopifnot(requireNamespace("COMA"))
  ploidy <- pop@ploidy
  
  effects <- read.csv(geno.file)[,1:2]
  colnames(effects) <- c("marker","add")
  
  geno <- t(pullSegSiteGeno(pop, simParam=SP))[effects$marker,]
  n <- min(ncol(geno), n.parent)
    
  if (weighted) {
    af <- apply(geno,1,mean)/ploidy
    p <- ifelse(effects$add > 0, af,1-af)
    effects$effects <- ifelse(p > 1e-6, effects$add/sqrt(p), 0)
  } else {
    effects$effects <- effects$add
  }
  GEBV <- as.numeric(crossprod(geno,effects$effects))
  tmp <- sort(GEBV,decreasing=TRUE,index.return=TRUE)
  parents <- colnames(geno)[tmp$ix[1:n]]
  
  matings <- expand.grid(parent1=parents,parent2=parents,stringsAsFactors = F)
  om <- rep(1/n^2,n^2)
  m <- n^2
  om.cumulative <- c(0,apply(array(1:m),1,function(i){sum(om[1:i])}))
  matings$n.progeny <- as.integer(table(cut(runif(n.progeny),
                                            breaks=om.cumulative)))
  ix <- which(matings$n.progeny > 0)
  matings <- matings[ix,]
  selected <- union(matings$parent1,matings$parent2)
  parent.pop <- selectInd(setPheno(pop,varE=0,simParam = SP), 
                          nInd=length(selected), 
                          candidates=match(selected,pop@id),
                          simParam = SP)
  return(list(parents=parent.pop, matings=matings))
}
