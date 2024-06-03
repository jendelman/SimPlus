#' Parent selection and mating by OCS in AlphaSimR
#'
#' Parent selection and mating by OCS in AlphaSimR
#'
#' Details
#' 
#' @param pop parental candidates, AlphaSimR Pop-class
#' @param SP simulation parameters for AlphaSimR
#' @param n.progeny number of progeny to simulate
#' @param dF target inbreeding rate
#' @param dF.max maximum inbreeding rate
#' @param geno.file geno filename for COMA
#' @param K.file kinship filename for COMA
#' @param min.c minimum contribution
#' @param solver CVXR solver
#' 
#' @return list containing
#' \describe{
#' \item{parents}{Pop-class of selected parents}
#' \item{ocs}{solution from ocs}
#' \item{matings}{data frame for sim_mate}
#' }
#'
#' @import AlphaSimR
#' @importFrom stats runif
#' @export

sim_OCS <- function(pop, SP, n.progeny, dF, dF.max, geno.file, K.file, 
                    min.c=1e-3, solver="ECOS") {

  stopifnot(inherits(pop,"Pop"))
  stopifnot(requireNamespace("COMA"))
  ploidy <- pop@ploidy

  ans1 <- COMA::read_data(geno.file=geno.file,
                    kinship.file=K.file,
                    ploidy=ploidy,matings="none",standardize=T)
  ans2 <- COMA::ocs(parents=data.frame(ans1$parents,min=0,max=1),
              ploidy=ploidy,K=ans1$K,dF=dF,min.c = min.c,
              dF.adapt=list(step=0.005,max=dF.max),
              solver=solver)
  if (nrow(ans2$oc)==0) {
    stop("No solution possible.")
  }
  parents <- ans2$oc$id
  
  oc <- ans2$oc$value
  matings <- expand.grid(parent1=parents,parent2=parents,stringsAsFactors = F)
  om <- oc[match(matings$parent1,parents)]*oc[match(matings$parent2,parents)]
  m <- length(om)
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
  return(list(parents=parent.pop, ocs=ans2, matings=matings))
}
