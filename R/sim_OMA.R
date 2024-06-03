#' Parent selection and mating by OMA in AlphaSimR
#'
#' Parent selection and mating by OMA in AlphaSimR
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
#' @param max.parent maximum number of candidate parents for oma
#' @param min.a minimum mate allocation
#' @param solver CVXR solver
#' 
#' @return list containing
#' \describe{
#' \item{parents}{Pop-class of selected parents}
#' \item{oma}{solution from oma}
#' \item{matings}{data frame for sim_mate}
#' \item{dev}{mean absolute deviation from random mating}
#' }
#'
#' @import AlphaSimR
#' @export

sim_OMA <- function(pop, SP, n.progeny, dF, dF.max, geno.file, K.file, 
                    max.parent, min.a=1e-3, solver="ECOS") {

  stopifnot(inherits(pop,"Pop"))
  stopifnot(requireNamespace("COMA"))
  ploidy <- pop@ploidy

  ans1 <- COMA::read_data(geno.file=geno.file,
                    kinship.file=K.file,
                    ploidy=ploidy,matings="none",standardize=T)
  ans2 <- COMA::ocs(parents=data.frame(ans1$parents,min=0,max=1/max.parent),
              ploidy=ploidy,K=ans1$K,dF=dF,
              dF.adapt=list(step=0.005,max=dF.max),
              solver=solver)
  if (nrow(ans2$oc)==0) {
    stop("No solution possible.")
  }
  sel1 <- ans2$oc$id[order(ans2$oc$value,decreasing=T)]
  if (length(sel1) > max.parent)
    sel1 <- sel1[1:max.parent]
  
  ans1 <- COMA::read_data(geno.file=geno.file,
                    kinship.file=K.file,
                    ploidy=ploidy,matings=sel1,standardize=T)
  ans3 <- COMA::oma(parents=data.frame(id=sel1,min=0,max=1),
              matings=data.frame(ans1$matings,min=0,max=1),
              ploidy=ploidy,K=ans1$K,dF=dF,min.a=min.a,
              dF.adapt=list(step=0.005,max=dF.max),
              solver=solver)
  if (nrow(ans3$om)==0) {
    stop("No solution possible.")
  }
  
  sel2 <- ans3$oc$id
  parent.pop <- selectInd(setPheno(pop,varE=0,simParam = SP), 
                          nInd=length(sel2), 
                          candidates=match(sel2,pop@id),simParam = SP)
  matings <- data.frame(ans3$om[,c("parent1","parent2")],
                        n.progeny=round(ans3$om$value*n.progeny,0))
  
  #calculate deviation from random mating
  np <- nrow(ans3$oc)
  y <- x <- matrix(0,nrow=np,ncol=np,dimnames=list(ans3$oc$id,ans3$oc$id))
  x[cbind(ans3$om$parent1,ans3$om$parent2)] <- ans3$om$value/2
  x <- x + t(x)
  tmp <- expand.grid(p1=ans3$oc$id,p2=ans3$oc$id,stringsAsFactors = F)
  ix1 <- match(tmp$p1,ans3$oc$id)
  ix2 <- match(tmp$p2,ans3$oc$id)
  tmp$value <- ans3$oc$value[ix1]*ans3$oc$value[ix2]/2
  y[cbind(tmp$p1,tmp$p2)] <- tmp$value 
  y <- y + t(y)
  
  return(list(parents=parent.pop, oma=ans3, 
              matings=matings, dev=mean(abs(x-y))))
}
