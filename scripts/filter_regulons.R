library(tidyverse)
library(optparse)

option_list = list(
  make_option(c("-r", "--regfile"), type="character", default=NULL, 
              help="Path to file containing regulons", metavar="character"))
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
regfile <- opt$regfile
outfile <- gsub("\\.csv","_filtered\\.csv",regfile)
motifsDf <- read.delim(regfile, header = F, sep=",")
topframe <- motifsDf[1:3,]
motifsDfsub <- motifsDf %>% dplyr::filter(V4 >= 3,
                            V7 == "gene is directly annotated" | grepl("gene is orthologous to.*which is directly annotated for motif",V7))
qcols <- grep(",",motifsDfsub[1,])
for (coli in qcols){
  motifsDfsub[,coli] <- paste0("\"",motifsDfsub[,coli],"\"")
}
motifsDfsub <- rbind(topframe, motifsDfsub)
write.table(motifsDfsub, file = outfile, sep = ",", quote = F, row.names = F, col.names = F)
