#' Parent selection and mating by OMA in AlphaSimR
#'
#' Parent selection and mating by OMA in AlphaSimR
#' 
#' Argument \code{dF} is a numeric vector of length 2, for lower and upper bounds on inbreeding rate.
#' 
#' Parent pre-selection is done with \code{ocs} to satisfy the \code{max.parent} limit for \code{oma}.
#' 
#' @param pop parental candidates, AlphaSimR Pop-class
#' @param SP simulation parameters for AlphaSimR
#' @param n.progeny number of progeny to simulate
#' @param dF target inbreeding rate
#' @param COMA.file geno filename for COMA
#' @param K.file kinship filename for COMA
#' @param max.parent maximum number of candidate parents for oma
#' @param solver CVXR solver
#' 
#' @return list containing
#' \describe{
#' \item{response}{named vector with realized dF, merit, Shannon diversity}
#' \item{oc}{data frame of optimal contributions for each individual}
#' \item{om}{data frame of optimal allocations for each mating}
#' }
#'
#' @import AlphaSimR
#' @import COMA
#' @export

sim_OMA <- function(pop, SP, n.progeny, dF, COMA.file, K.file, 
                    max.parent, solver="ECOS") {

  stopifnot(inherits(pop,"Pop"))
  stopifnot(requireNamespace("COMA"))
  ploidy <- pop@ploidy
  
  stopifnot(length(dF)==2L)
  stopifnot(dF[1]<=dF[2])

  ans1 <- COMA::read_data(geno.file=COMA.file,
                    kinship.file=K.file,
                    ploidy=ploidy,matings="none",standardize=T)
  ans2 <- COMA::ocs(parents=data.frame(ans1$parents,min=0,max=1/max.parent),
              ploidy=ploidy,K=ans1$K,dF=dF[2],
              dF.adapt=list(step=0.005,max=0.1),
              solver=solver)
  if (nrow(ans2$oc)==0) {
    stop("No solution possible.")
  }
  sel1 <- ans2$oc$id[order(ans2$oc$value,decreasing=T)]
  if (length(sel1) > max.parent)
    sel1 <- sel1[1:max.parent]
  
  ans1 <- COMA::read_data(geno.file=COMA.file,
                    kinship.file=K.file,
                    ploidy=ploidy,matings=sel1,standardize=T)
  ans3 <- COMA::oma(parents=data.frame(id=sel1,min=0,max=1),
              matings=data.frame(ans1$matings,min=0,max=1),
              ploidy=ploidy,K=ans1$K,dF=dF,
              dF.adapt=list(step=0.005,max=0.1),
              solver=solver)
  if (nrow(ans3$om)==0) {
    stop("No solution possible.")
  }
  
  return(ans3)
}
