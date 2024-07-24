#' Estimate marker effects by StageWise 
#'
#' Estimate marker effects by StageWise 
#'
#' Details
#' 
#' @param geno.file TP genotype CSV filename
#' @param pheno.file TP phenotype CSV filename
#' @param ploidy ploidy
#' @param COMA.file input filename for COMA
#' @param gen.TP number of generations for TP (including current)
#' @param asreml.workspace memory workspace
#' 
#' @return list of dominance parameters
#'
#' @export
#' @importFrom data.table fread
#' @importFrom dplyr filter
#' @importFrom utils write.csv

sim_StageWise <- function(geno.file, pheno.file, ploidy, COMA.file, 
                          gen.TP, asreml.workspace="500mb") {

  stopifnot(gen.TP >= 1)
  g1 <- StageWise::read_geno(geno.file,ploidy=ploidy,map=F,
                             min.minor.allele = 1,dominance=T)
  pheno <- read.csv(pheno.file,colClasses = c("character","numeric","integer"))
  cyc <- max(pheno$gen)
  pheno <- filter(pheno, gen >= cyc - (gen.TP-1))
  data <- data.frame(env="none",id=pheno$id,
                     BLUE=pheno$value,
                     stringsAsFactors = F)
  a2 <- StageWise::Stage2(data=data,geno=g1,non.add="dom",silent=F,
                          workspace = asreml.workspace)
  prep <- StageWise::blup_prep(data=data,geno=g1,vars=a2$vars)
  AM <- StageWise::blup(prep,what="AM",geno=g1)
  DM <- StageWise::blup(prep,what="DM",geno=g1)
  marker.effects <- data.frame(marker=AM$marker,
                               add=AM$effect,
                               dom=DM$effect)
  tmp <- fread(geno.file,sep=",",header=T)
  geno <- as.matrix(tmp[,-1])
  rownames(geno) <- as.character(tmp$V1)
  
  pheno <- filter(pheno, gen == cyc)
  write.csv(cbind(marker.effects,
                  geno[marker.effects$marker,pheno$id]),
            file=COMA.file,row.names=F)

  x <- StageWise::dominance(a2$params)
  return(list(Vd=x$estimate[1],b=x$estimate[2]))
}
