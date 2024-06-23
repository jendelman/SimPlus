#' Estimate selection accuracy
#' 
#' Estimate selection accuracy
#' 
#' @param pop parental candidates, AlphaSimR Pop-class
#' @param SP simulation parameters for AlphaSimR
#' @param COMA.file COMA filename
#' @param K.file kinship filename for COMA
#' @param n.mate number of matings to sample
#' @param n.progeny number of progeny per mating
#' 
#' @return list of accuracies for OCS and OMA
#' @importFrom stats cor
#' @export

sim_accuracy <- function(pop, SP, COMA.file, K.file=NULL, 
                         n.mate=100, n.progeny=50) {

  stopifnot(inherits(pop,"Pop"))
  stopifnot(requireNamespace("COMA"))
  
  matings <- expand.grid(parent1=pop@id, parent2=pop@id,stringsAsFactors = F)
  ix <- sample(nrow(matings), n.mate)
  matings <- matings[ix,]
  data <- fread(COMA.file, header=T, sep=",")
  if (colnames(data)[1]=="id") {
    pheno <- TRUE
  } else {
    stopifnot(!is.null(K.file))
    pheno <- FALSE
  }
  
  if (pheno) {
    data <- data.frame(data[,1:2])
    pred <- matings
    pred$GCA <- (data$value[match(pred$parent1,data$id)]+
                     data$value[match(pred$parent2,data$id)])/2
  } else {
    data <- COMA::read_data(geno.file=COMA.file,
                          kinship.file=K.file,
                          ploidy=pop@ploidy,
                          matings=matings)
    pred <- data$matings
    pred$GCA <- (data$parents$merit[match(pred$parent1,data$parents$id)] + 
                   data$parents$merit[match(pred$parent1,data$parents$id)])/2
  }
  
  matings$value <- 1/n.mate
  data2 <- sim_mate(pop, SP, matings, total.progeny=n.mate*n.progeny, 
                    min.progeny=1)
  matings <- merge(pred,data2$matings)
  if (pheno) {
    return(list(ocs=cor(matings$GCA,matings$mu), 
                oma=as.numeric(NA)))
  } else {
    return(list(ocs=cor(matings$GCA,matings$mu),
                oma=cor(matings$merit,matings$mu)))
  }
}
    