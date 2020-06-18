wordir = '/Users/Thomsn/Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/phylogenetto/data/real_data/vitis/hyde/admixed_asia/'
setwd(wordir)
bonfs = read.csv2('2020_04_vitisasia-ind-bonfiltered.txt', header = T, sep = '\t')
data = read.csv2('suspiciousindividuals.csv', header = T, sep = ';', dec = ',')
data$ADN.ID = tolower(data$ADN.ID)
bonfs[,5] = unname(sapply(sapply(bonfs[,5], as.character), as.numeric))
bonfs[,4] = unname(sapply(sapply(bonfs[,4], as.character), as.numeric))
bonfs[,6] = unname(sapply(sapply(bonfs[,6], as.character), as.numeric))
sorted = bonfs[order(bonfs[,5]),]

write("The following individuals show admixture between european and asian/american grapewines. \
The abbreviated p1-taxa are: acerif(olia) [USA], aesti(valis) [USA], pias(zekii) [Asia],\
coig(netiae) [Asia], thun(bergii) [Asia] and amur(ensis) [Asia].", file = "admixed_ind_hyde.txt")

for (ind in unique(sorted$Hybrid)){
  subsort = sorted[sorted$Hybrid == ind,]
  write('\n', file = "admixed_ind_hyde.txt", append = T)
  write(as.character(data[data$ADN.ID == ind,]$Name), file = "admixed_ind_hyde.txt", append = T)
  write.table(format(subsort[,c(1:6)], digits=5), sep="\t\t", file = "admixed_ind_hyde.txt", append = T, quote = F, row.names = F)
}

