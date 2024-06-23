#' Load simulation files
#'
#' Load simulation files
#' 
#' @param rda.file R data filename
#' @param geno.file genotype CSV file
#' @param pheno.file phenotype CSV filename
#' @param ped.file pedigree CSV filename
#'
#' @importFrom data.table fwrite
#' @importFrom utils write.csv
#' @export
#' 
#' @return list containing pop,SP,p.ref

sim_load <- function(rda.file, geno.file, pheno.file, ped.file) {

  load(rda.file)
  fwrite(geno, file=geno.file)
  write.csv(ped, file=ped.file, row.names=F)
  write.csv(pheno, file=pheno.file, row.names=F)
  return(list(pop=pop,SP=SP,p.ref=p.ref))
}
