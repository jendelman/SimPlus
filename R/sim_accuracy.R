#' Estimate selection accuracy
#' 
#' Estimate selection accuracy
#' 
#' Details
#' 
#' @param pop parental candidates, AlphaSimR Pop-class
#' @param SP simulation parameters for AlphaSimR
#' @param geno.file geno filename for COMA
#' @param K.file kinship filename for COMA
#' @param n.mate number of matings to sample
#' @param n.progeny number of progeny per mating
#' 
#' @return list of accuracies for OCS and OMA
#' @importFrom stats cor
#' @export

sim_accuracy <- function(pop, SP, COMA.file, K.file, n.mate=100, n.progeny=50) {

  stopifnot(inherits(pop,"Pop"))
  stopifnot(requireNamespace("COMA"))
  
  matings <- expand.grid(parent1=pop@id, parent2=pop@id,stringsAsFactors = F)
  ix <- sample(nrow(matings), n.mate)
  matings <- matings[ix,]
  data <- COMA::read_data(geno.file=COMA.file,
                          kinship.file=K.file,
                          ploidy=pop@ploidy,
                          matings=matings)
  
  matings$value <- 1/n.mate
  data2 <- sim_mate(pop, SP, matings, total.progeny=n.mate*n.progeny, 
                    min.progeny=1)
  matings <- merge(data$matings,data2$matings)
  acc.oma <- cor(matings$merit, matings$mu)
  ocs.pred <- (data$parents$merit[match(matings$parent1,data$parents$id)] + 
    data$parents$merit[match(matings$parent1,data$parents$id)])/2
  acc.ocs <- cor(ocs.pred,matings$mu)
  return(list(ocs=acc.ocs, oma=acc.oma))
}
    