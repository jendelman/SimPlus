#' Update the Training Population (TP)
#'
#' Update genotype, genotyped, and pedigree files for the TP
#'
#' Pedigree file has three columns: id, mother, father. Phenotype file has three columns: id, value, gen. To start a new TP, use \code{gen.TP} = 0. 
#' 
#' @param pop variable of AlphaSimR Pop-class
#' @param SP simulation parameters for AlphaSimR
#' @param geno.file genotype CSV file
#' @param pheno.file phenotype CSV filename
#' @param ped.file pedigree CSV filename
#' @param gen.TP # generations recorded for the TP (including current)
#' 
#' @return vector of allele frequencies for pop if gen.TP==0, otherwise NULL
#'
#' @import AlphaSimR
#' @importFrom dplyr filter
#' @importFrom data.table fwrite fread
#' @importFrom pedigree trimPed
#' @export

update_TP <- function(pop, SP, geno.file, pheno.file, ped.file, gen.TP) {

  stopifnot(inherits(pop,"Pop"))
  
  if (gen.TP > 0) {
    pheno <- read.csv(pheno.file)
    cyc <- max(pheno$gen) + 1
    new.pheno <- data.frame(id=pop@id, value=pheno(pop)[,1], gen=cyc)
    pheno <- filter(rbind(pheno,new.pheno), gen >= cyc - (gen.TP-1))
  } else {
    pheno <- data.frame(id=pop@id, value=pheno(pop)[,1], gen=0)
  }
  write.csv(pheno, file=pheno.file, row.names=F)
  
  if (gen.TP > 0) {
    ped <- read.csv(ped.file)
  } else {
    ped <- NULL
  }
  ped <- rbind(ped, getPed(pop))
  keep <- pedigree::trimPed(ped, ped$id %in% pheno$id)
  ped <- ped[keep,]
  write.csv(ped, file=ped.file, row.names=F)
  
  new.geno <- t(pullSegSiteGeno(pop, simParam=SP))
  if (gen.TP > 0) {
    geno <- as.matrix(fread(geno.file, sep=",", header = T, drop=1))
    geno2 <- cbind(geno,new.geno)[,pheno$id]
    rownames(geno2) <- rownames(new.geno)
  } else {
    p.ref <- apply(new.geno,1,mean)/pop@ploidy
    geno2 <- new.geno
  }
  suppressMessages(fwrite(geno2, file=geno.file, row.names=T))
  
  if (gen.TP==0) {
    return(p.ref)
  } else {
    return(NULL)
  }
}
