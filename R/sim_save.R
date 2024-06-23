#' Save simulation files
#'
#' Save simulation files
#' 
#' @param rda.file R data filename
#' @param pop variable of AlphaSimR Pop-class
#' @param SP simulation parameters for AlphaSimR
#' @param geno.file genotype CSV file
#' @param pheno.file phenotype CSV filename
#' @param ped.file pedigree CSV filename
#' @param p.ref reference allele frequencies (for G.VR1)
#'
#' @importFrom data.table fread
#' @importFrom utils read.csv
#' @export

sim_save <- function(rda.file, pop, SP, geno.file, pheno.file, ped.file, 
                     p.ref) {

  geno <- fread(geno.file, sep=",", header=T)
  pheno <- read.csv(pheno.file,check.names=F)
  ped <- read.csv(ped.file,check.names=F)
  save(pop,SP,p.ref,geno,pheno,ped,file=rda.file)
}
