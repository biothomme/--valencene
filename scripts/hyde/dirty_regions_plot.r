wordir = '/Users/Thomsn/Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/phylogenetto/data/real_data/vitis/hyde/plot/bonf_.05/'
setwd(wordir)
library(stringr)
df1 = read.csv('2020_04_vitis_imputed_hyde_hybr_ind_trp.csv', header=T)
ndf = read.csv('2020_04_vitis_imputed_hyde_hybr_not_sign_in_trp.csv', header=T)

oo = tolower(as.character(ndf$Pop.celine))
uu = apply(as.matrix(oo), 1,function(x) str_replace_all(x, "\\.", ""))
ndf$Pop.celine = apply(as.matrix(uu), 1,function(x) str_replace_all(x, "\\ ", ""))

signi = c(rep(TRUE, nrow(df1)), rep(FALSE, nrow(ndf)))
taxon = c(as.character(df1$hybrid), ndf$Pop.celine)
regions = c(as.character(df1$hyb_region), as.character(ndf$region.geo))
species = c(as.character(df1$hyb_taxon), as.character(ndf$taxon))

rm(hybrid_significance)
for (tax in unique(taxon)){
  loci = regions[tax == taxon]
  silo = signi[tax == taxon]
  for (reg in unique(loci)){
    subsilo = silo[reg == loci]
    perca = as.numeric(sum(subsilo) / length(subsilo))
    rowns = c('hybrid', 'region', 'sign_ind', 'total_ind', 'ratio_sign')
    dafo = data.frame('hybrid' = tax, 
                      'region' = reg, 
                      'sign_ind' = sum(subsilo), 
                      'total_ind' = length(subsilo), 
                      'ratio_sign' = perca)
    if (exists('hybrid_significance')){
      hybrid_significance = rbind(hybrid_significance, dafo)
    } else {
      hybrid_significance = dafo
    }
  }
}
percentage = c(hybrid_significance$ratio_sign, c(1 - hybrid_significance$ratio_sign))
region = c(as.character(hybrid_significance$region),as.character(hybrid_significance$region))
hyde_result = c(rep('sign. hybridization', nrow(hybrid_significance)), rep('no hybridization', nrow(hybrid_significance)))

library(ggplot2)

# specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
# condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
# value <- abs(rnorm(12 , 0 , 15))
data <- data.frame(region,hyde_result,percentage)

# Stacked + percent
ggplot(data, aes(fill=hyde_result, y=percentage, x=region)) + 
  geom_bar(position="fill", stat="identity")


