#' Parent selection and mating by OCS in AlphaSimR
#'
#' Parent selection and mating by OCS in AlphaSimR
#'
#' Auto detects whether COMA.file has marker effects or phenotype data based on 
#' the first column name.
#' 
#' @param pop parental candidates, AlphaSimR Pop-class
#' @param SP simulation parameters for AlphaSimR
#' @param n.progeny number of progeny to simulate
#' @param dF target inbreeding rate
#' @param dF.max maximum inbreeding rate
#' @param COMA.file marker effect or pheno filename for COMA
#' @param K.file kinship filename for COMA
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

sim_OCS <- function(pop, SP, n.progeny, dF, dF.max, COMA.file, K.file, 
                    solver="ECOS") {

  stopifnot(inherits(pop,"Pop"))
  stopifnot(requireNamespace("COMA"))
  ploidy <- pop@ploidy

  data <- fread(COMA.file, header=T, sep=",")
  if (colnames(data)[1]=="id") {
    data2 <- fread(K.file, header=T, sep=",")
    K <- as.matrix(data2[,-1])
    rownames(K) <- as.character(data2$V1)
    data <- data.frame(data[,1:2])
    colnames(data) <- c("id","merit")
    ans1 <- list(parents=data[match(rownames(K),data$id),], K=K)
  } else {
    ans1 <- COMA::read_data(geno.file=COMA.file,
                    kinship.file=K.file,
                    ploidy=ploidy,matings="none",standardize=T)
  } 
  
  ans2 <- COMA::ocs(parents=data.frame(ans1$parents,min=0,max=1),
              ploidy=ploidy,K=ans1$K,dF=dF,
              dF.adapt=list(step=0.005,max=dF.max),
              solver=solver)
  if (nrow(ans2$oc)==0) {
    stop("No solution possible.")
  }
  
  parents <- ans2$oc$id
  oc <- ans2$oc$value
  ans2$om <- expand.grid(parent1=parents,parent2=parents,stringsAsFactors = F)
  ans2$om$value <- oc[match(ans2$om$parent1,parents)]*oc[match(ans2$om$parent2,parents)]
  
  return(ans2)
}
