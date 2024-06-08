#' Perform mating in AlphaSimR
#'
#' Perform mating in AlphaSimR
#'
#' Details
#' 
#' @param pop parental candidates, AlphaSimR Pop-class
#' @param SP simulation parameters for AlphaSimR
#' @param matings data frame of mate allocations, with columns parent1, parent2, value
#' @param total.progeny total progeny to generate
#' @param min.progeny minimum progeny per mating
#' 
#' @return list containing
#' \describe{
#' \item{progeny}{Pop-class of progeny}
#' \item{matings}{matings with additional column for mu}
#' }
#'
#' @import AlphaSimR
#' @importFrom stats runif
#' @export

sim_mate <- function(pop, SP, matings, total.progeny, min.progeny) {

  stopifnot(inherits(pop,"Pop"))
  f1 <- function(x,n) {
    m <- length(x)
    cdf.x <- c(0,apply(array(1:m),1,function(i){sum(x[1:i])}))
    as.integer(table(cut(runif(n),breaks=cdf.x)))
  }
  
  x <- matings$value <- matings$value/sum(matings$value)
  
  n <- floor(total.progeny/min.progeny)
  matings$n.progeny <- f1(x,n)*min.progeny + f1(x,1)*(total.progeny%%min.progeny)
  matings <- matings[matings$n.progeny > 0,]
  m <- nrow(matings)
  selected <- union(matings$parent1,matings$parent2)
  n.parent <- length(selected)
#  parent.pop <- selectInd(setPheno(pop,varE=0,simParam = SP), 
#                          nInd=length(selected), 
#                          candidates=match(selected,pop@id),
#                          simParam = SP)
  
  crossPlan <- cbind(match(matings$parent1, pop@id),
                     match(matings$parent2, pop@id))
  progeny.pop <- makeCross(pop=pop,
                           crossPlan=crossPlan[1,,drop=FALSE],
                           nProgeny=matings$n.progeny[1],simParam = SP)
  matings$mu <- numeric(m)
  matings$mu[1] <- genParam(progeny.pop, simParam = SP)$mu
  
  if (m > 1) {
    for (k in 2:m) {
      tmp <- makeCross(pop=pop,
                crossPlan=crossPlan[k,,drop=FALSE],
                nProgeny=matings$n.progeny[k],simParam = SP)
      progeny.pop <- c(progeny.pop,tmp)
      matings$mu[k] <- genParam(tmp, simParam = SP)$mu
    }
  }
  return(list(n.parent=n.parent, progeny=progeny.pop, matings=matings))
}
