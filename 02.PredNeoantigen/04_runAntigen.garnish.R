library(magrittr)
library(data.table)
library(antigen.garnish)

pep <- read.table("aml_0201.txt", header=FALSE)
pep_out <- foreignness_score(pep, db = "human")
write.table(pep_out, "pep_out.txt", sep = "\t")