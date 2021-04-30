rm(list = objects())
library("INLA")
library("tidyverse")
options(warn=-2)

REGIONE<-"lombardia"

MESE<-3
readRDS(glue::glue("result{MESE}_{REGIONE}.RDS"))->result
readRDS(glue::glue("dati{MESE}_{REGIONE}.RDS"))->dati
readRDS(glue::glue("mesh{MESE}_{REGIONE}.RDS"))->mesh


inla.posterior.sample(n=1000,result,num.threads = 3)->inlaSampleOut

saveRDS(inlaSampleOut,glue::glue("inlaSampleOut_mese{MESE}_{REGIONE}.RDS"))
