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
#' @return numeric
#' @importFrom stats cor
#' @export


sim_acc <- function(pop, SP, geno.file, K.file, n.mate=200, n.progeny=50) {

  stopifnot(inherits(pop,"Pop"))
  stopifnot(requireNamespace("COMA"))
  
  matings <- expand.grid(parent1=pop@id, parent2=pop@id,stringsAsFactors = F)
  ix <- sample(nrow(matings), n.mate)
  matings <- matings[ix,]
  ploidy <- pop@ploidy
      
  data <- COMA::read_data(geno.file=geno.file,
                          kinship.file=K.file,
                          ploidy=ploidy,
                          matings=matings)
  matings$n.progeny <- n.progeny
  data2 <- sim_mate(pop, SP, matings)
  cor(data$matings$merit, data2$matings$mu)
}
    