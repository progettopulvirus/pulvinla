rm(list=objects())
library("INLA")
library("tidyverse")
library("tictoc")

MESE<-4
REGIONE<-"lombardia"

inla.setOption(pardiso.license="~/pardiso/licenza.txt")


readRDS(glue::glue("result{MESE}_{REGIONE}.RDS"))->inla.out
tic()
inla.posterior.sample(n=1000,result=inla.out,num.threads = 3)->inlaSampleOut
toc()
saveRDS(inlaSampleOut,glue::glue("inlaSampleOut_mese{MESE}_{REGIONE}.RDS"))