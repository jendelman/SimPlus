#' Simulate one cycle of mass selection
#'
#' Simulate one cycle of mass selection
#'
#' Options for selection method are "PS", "OCS", "OMA", "WGEBV". Results appended to \code{results.file}. Cycle progression: genomic prediction, mate selection, mating, phenotyping, updating TP files. The \code{params} variable is a list with the following: geno.file, pheno.file, K.file, results.file, K.method, gen.TP, gen.record, p.ref, ibd.loci, max.parent, solver, n.core. See \code{\link{sim_TP}}, \code{\link{sim_StageWise}}, \code{\link{sim_OCS}}, \code{\link{sim_OMA}} for details. 
#'
#' @param dF inbreeding rate
#' @param sel.method selection method
#' @param total.progeny number of progeny for next generation
#' @param min.progeny minimum progeny per mating
#' @param pop starting population
#' @param SP AlphaSimR parameters
#' @param params other parameters (Details)
#' 
#' @return list containing
#' \describe{
#' \item{pop}{Pop-class of progeny}
#' \item{SP}{simulation parameters}
#' \item{stats}{list with n.parent, dF, F.A, F.G, accuracy}
#' }
#'
#' @import AlphaSimR
#' @export

sim_mass <- function(dF, sel.method, total.progeny, min.progeny, 
                     pop, SP, params){

  stopifnot(inherits(pop,"Pop"))
  
  if (sel.method=="PS") {
    ans2 <- sim_OCS(pop, SP, n.progeny=total.progeny, dF=dF, dF.max=0.1, 
                     COMA.file=params$pheno.file, 
                     K.file=params$K.file,
                     solver=params$solver)
    acc <- as.numeric(NA)
  } else {
    ans1 <- sim_StageWise(params$geno.file,
                          params$pheno.file,
                          pop@ploidy,
                          params$COMA.file,
                          params$gen.record,
                          params$asreml.workspace)
    tmp <- sim_accuracy(pop, SP, params$COMA.file, params$K.file)
    if (K.method=="OCS") {
      acc <- tmp$acc.ocs
    } else {
      acc <- tmp$acc.oma
    }
  
    if (sel.method=="OCS") {
      ans2 <- sim_OCS(pop, SP, n.progeny=total.progeny, dF=dF, dF.max=0.1, 
                       COMA.file=params$COMA.file, 
                       K.file=params$K.file,
                       solver=params$solver)
      COMA.dF <- ans2$response$dF
    }
    if (sel.method=="OMA") {
      ans2 <- sim_OMA(pop, SP, n.progeny=total.progeny, dF=dF, dF.max=0.1, 
                      COMA.file=params$COMA.file, 
                      K.file=params$K.file,
                      max.parent=params$max.parent,
                      solver=params$solver)
      COMA.dF <- ans2$response$dF2
    }
  }

  ans3 <- sim_mate(pop, SP, matings=ans2$om, total.progeny, min.progeny)
  
  pop <- setPheno(ans3$progeny,simParam=SP)
  ans4 <- sim_TP(pop, SP,
                   params$geno.file,
                   params$pheno.file,
                   params$ped.file,
                   params$K.file,
                   params$K.method,
                   params$gen.record,
                   params$p.ref,
                   params$ibd.loci,
                   params$n.core)
  
  return(list(pop=pop,SP=SP,
          stats=list(n.parent=ans3$n.parent, n.mate=nrow(ans3$matings),
                     F.A=ans4$F.A, F.G=ans4$F.G, accuracy=acc, dF=COMA.dF)))
}
