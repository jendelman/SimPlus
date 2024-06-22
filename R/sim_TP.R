#' Add new data to the TP
#'
#' Add new data to the TP
#'
#' Pedigree file has three columns: id, mother, father. Phenotype file has three columns: id, value, gen. To create a new TP, use \code{gen.record} = 0. Argument \code{p.ref} is only used when \code{K.method}="G". 
#' 
#' @param pop variable of AlphaSimR Pop-class
#' @param SP simulation parameters for AlphaSimR
#' @param geno.file genotype CSV file
#' @param pheno.file phenotype CSV filename
#' @param ped.file pedigree CSV filename
#' @param K.file kinship CSV filename
#' @param K.method kinship method
#' @param gen.record number of generations to record (including current)
#' @param p.ref reference pop allele frequencies for "VR1" method
#' @param ibd.loci number of loci per chromosome for "AM" method 
#' @param n.core number cores for "AM" method
#' 
#' @return list containing
#' \describe{
#' \item{F.A}{pedigree inbreeding coefficient}
#' \item{F.G}{genomic inbreeding coefficient}
#' \item{p.ref}{allele freq in reference population}
#' }
#'
#' @import AlphaSimR
#' @importFrom dplyr filter
#' @importFrom data.table fwrite fread
#' @importFrom pedigree trimPed
#' @importFrom polyBreedR A_mat G_mat
#' @export

sim_TP <- function(pop, SP, geno.file, pheno.file, ped.file, K.file, K.method, 
                   gen.record, p.ref=NULL, ibd.loci=NULL, n.core=1) {

  stopifnot(inherits(pop,"Pop"))
  
  if (gen.record > 0) {
    pheno <- read.csv(pheno.file)
    cyc <- max(pheno$gen) + 1
    new.pheno <- data.frame(id=pop@id, value=pheno(pop)[,1], gen=cyc)
    pheno <- filter(rbind(pheno,new.pheno), gen >= cyc - (gen.record-1))
  } else {
    pheno <- data.frame(id=pop@id, value=pheno(pop)[,1], gen=0)
  }
  write.csv(pheno, file=pheno.file, row.names=F)
  
  if (gen.record > 0) {
    ped <- read.csv(ped.file)
  } else {
    ped <- NULL
  }
  ped <- rbind(ped, getPed(pop))
  keep <- pedigree::trimPed(ped, ped$id %in% pheno$id)
  ped <- ped[keep,]
  write.csv(ped, file=ped.file, row.names=F)
  
  new.geno <- t(pullSegSiteGeno(pop, simParam=SP))
  if (gen.record > 0) {
    geno <- as.matrix(fread(geno.file, sep=",", header = T, drop=1))
    geno2 <- cbind(geno,new.geno)[,pheno$id]
    rownames(geno2) <- rownames(new.geno)
  } else {
    p.ref <- apply(new.geno,1,mean)/pop@ploidy
    geno2 <- new.geno
  }
  suppressMessages(fwrite(geno2, file=geno.file, row.names=T))
  
  A <- A_mat(ped=ped, ploidy=pop@ploidy, order.ped=F)[pop@id,pop@id]
  F.A <- (mean(diag(A))-1)/(pop@ploidy-1)
  F.G <- as.numeric(NA)
  if (K.method=="G.VR1") {
    A <- G_mat(new.geno, pop@ploidy, p.ref=p.ref, method="VR1")
    F.G <- (mean(diag(A))-1)/(pop@ploidy-1)
  } 
  if (K.method=="G.IBD") {
    marks <- unlist(lapply(SP$genMap,function(map){
      u <- split(names(map),cut(map,breaks=ibd.loci))
      sapply(u,"[[",1)
    }))
    geno <- t(pullIbdHaplo(pop, simParam=SP))[marks,]
    A <- G_mat(geno, ploidy=pop@ploidy, method="AM",n.core=n.core)
    F.G <- (mean(diag(A))-1)/(pop@ploidy-1)
  }
  write.csv(A/pop@ploidy,K.file,row.names=T)

  return(list(F.A=F.A, F.G=F.G, p.ref=p.ref))
}
