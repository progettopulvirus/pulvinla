rm(list = objects())
library("INLA")
library("tidyverse")
options(warn=-2)

MESE<-4
readRDS(glue::glue("result{MESE}_lombardia.RDS"))->result
readRDS(glue::glue("dati{MESE}_lombardia.RDS"))->dati
readRDS(glue::glue("mesh{MESE}_lombardia.RDS"))->mesh


inla.posterior.sample(n=1000,result,num.threads = 3)->listaSamples

Aloc<-inla.spde.make.A(mesh=mesh,loc=cbind(dati$coordx/1000,dati$coordy/1000),group = dati$banda,n.group = length(unique(dati$date)))

inla.posterior.sample.eval(function(...,A=Aloc,df=dati){
  
  Intercept+as.numeric(A %*% i)->eta
  
  eta+weekend*df$weekend+d_a1*df$d_a1+d_a2*df$d_a2+
    i_surface*df$i_surface+day*df$day+daylockdown*df$daylockdown+
    altitudedem*df$altitudedem+cams*df$cams->eta
  
  # rnorm(n=nrow(dati),mean = 0,sd=sqrt(1/theta[1]))->noise
  # rnorm(n=nrow(dati),mean = 0,sd=sqrt(1/theta[2]))->sec
  
  eta

  
},samples = listaSamples)->dfOut

tibble(as.data.frame(dfOut))->dfOut
names(dfOut)<-paste0("sample",1:1000)
apply(dfOut,MARGIN = 1,FUN=mean,na.rm=TRUE)->media

bind_cols(dati[,c("value","date","banda","yy","mm","dd","station_eu_code","weekend")],dfOut)->dati2
dati2$media<-media

dati2 %>%
  mutate(banda=ifelse(banda>30,banda-30,banda))->dati2

write_delim(dati2,glue::glue("samples{MESE}.csv"),delim=";",col_names=TRUE)
