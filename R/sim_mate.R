#' Perform mating in AlphaSimR
#'
#' Perform mating in AlphaSimR
#'
#' Details
#' 
#' @param pop parental candidates, AlphaSimR Pop-class
#' @param SP simulation parameters for AlphaSimR
#' @param matings data frame with columns parent1, parent2, n.progeny
#' 
#' @return list containing
#' \describe{
#' \item{progeny}{Pop-class of progeny}
#' \item{matings}{matings input with additional column for mu}
#' }
#'
#' @import AlphaSimR
#' @export

sim_mate <- function(pop, SP, matings) {

  stopifnot(inherits(pop,"Pop"))
  crossPlan <- cbind(match(matings$parent1,pop@id),
                     match(matings$parent2,pop@id))
  n.mate <- nrow(crossPlan)
  progeny.pop <- makeCross(pop=pop,
                           crossPlan=crossPlan[1,,drop=FALSE],
                           nProgeny=matings$n.progeny[1],simParam = SP)
  matings$mu <- numeric(n.mate)
  matings$mu[1] <- genParam(progeny.pop, simParam = SP)$mu
  
  if (n.mate > 1) {
    for (k in 2:n.mate) {
      tmp <- makeCross(pop=pop,
                crossPlan=crossPlan[k,,drop=FALSE],
                nProgeny=matings$n.progeny[k],simParam = SP)
      progeny.pop <- c(progeny.pop,tmp)
      matings$mu[k] <- genParam(tmp, simParam = SP)$mu
    }
  }
  return(list(progeny=progeny.pop, matings=matings))
}
