#' Update kinship file
#'
#' Update kinship file
#'
#' Kinship is computed for selection candidates in \code{pop}. For \code{K.method}="A", argument \code{ped.file} is the three-column pedigree file. \code{K.method}="G" specifies VanRaden Method 1, and argument \code{p.ref} is the vector of reference allele frequencies. To reduce computing time for \code{K.method}="G.IBD", only \code{ibd.loci} are sampled per chromosome. 
#' 
#' @param pop variable of AlphaSimR Pop-class
#' @param SP simulation parameters for AlphaSimR
#' @param K.method kinship method: A, G, G.IBD
#' @param K.file kinship CSV filename
#' @param ped.file pedigree CSV filename
#' @param p.ref reference pop allele frequencies for "G" method
#' @param ibd.loci number of loci per chromosome for "G.IBD" method 
#' @param n.core number cores for "G.IBD" method
#' 
#' @return average inbreeding coefficient of the population
#'
#' @import AlphaSimR
#' @importFrom dplyr filter
#' @importFrom utils read.csv
#' @importFrom pedigree trimPed
#' @importFrom polyBreedR A_mat G_mat
#' @export

update_K <- function(pop, SP, K.method, K.file, ped.file=NULL, 
                     p.ref=NULL, ibd.loci=100, n.core=1) {

  stopifnot(inherits(pop,"Pop"))
  stopifnot(K.method %in% c("A","G","G.IBD"))
  
  if (K.method=="A") {
    ped <- read.csv(ped.file,colClasses = rep("character",3))
    A <- A_mat(ped=ped, ploidy=pop@ploidy, order.ped=F)[pop@id,pop@id]
  }
  
  if (K.method=="G") {
    geno <- t(pullSegSiteGeno(pop, simParam=SP))
    A <- G_mat(geno, pop@ploidy, p.ref=p.ref, method="VR1")
  } 
  
  if (K.method=="G.IBD") {
    marks <- unlist(lapply(SP$genMap,function(map){
      u <- split(names(map),cut(map,breaks=ibd.loci))
      sapply(u,"[[",1)
    }))
    geno <- t(pullIbdHaplo(pop, simParam=SP))[marks,]
    A <- G_mat(geno, ploidy=pop@ploidy, method="AM", n.core=n.core)
  }
  write.csv(A/pop@ploidy, K.file, row.names=T)

  return((mean(diag(A))-1)/(pop@ploidy-1))
}
