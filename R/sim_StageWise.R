#' Estimate marker effects by StageWise within AlphaSimR
#'
#' Estimate marker effects by StageWise within AlphaSimR
#'
#' Details
#' 
#' @param pop training population, AlphaSimR Pop-class
#' @param SP simulation parameters for AlphaSimR
#' @param geno.file geno filename for COMA
#' 
#' @return nothing
#'
#' @import AlphaSimR
#' @export

sim_StageWise <- function(pop, SP, geno.file) {

  stopifnot(inherits(pop,"Pop"))
  ploidy <- pop@ploidy
  geno <- t(pullSegSiteGeno(pop, simParam=SP))
  label <- sample(1e8:9e8,1)
  tmp.file <- sub("X",label,"sim_geno_X.csv")
  write.csv(geno,tmp.file)
  
  g1 <- StageWise::read_geno(tmp.file,ploidy=ploidy,map=F,min.minor.allele = 1,
                  dominance=T)
  data <- data.frame(env="none",id=colnames(geno),BLUE=pheno(pop)[,1])
  a2 <- StageWise::Stage2(data=data,geno=g1,non.add="dom",silent=F)
  prep <- StageWise::blup_prep(data=data,geno=g1,vars=a2$vars)
  AM <- StageWise::blup(prep,what="AM",geno=g1)
  DM <- StageWise::blup(prep,what="DM",geno=g1)
  marker.effects <- data.frame(marker=AM$marker,
                               add=AM$effect,
                               dom=DM$effect)
  write.csv(cbind(marker.effects,geno[marker.effects$marker,]),
            file=geno.file,row.names=F)
  file.remove(tmp.file)
  return()
}
