rm(list=objects())
library("tidyverse")
library("stazioniMonitoraggio")


list.files(pattern="^dati.+RDS$")->ffile


purrr::map_dfr(ffile,.f=function(nomeFile){
  
  readRDS(nomeFile)->dati
  ifelse(grepl("Train",nomeFile),"Train","Val")->dati$stage
  str_extract(str_extract(nomeFile,"_indice[1234]_"),"[[:digit:]]")->dati$trial
  
  str_extract(str_extract(nomeFile,"^dati[[:alpha:]]+[0-9]_"),"[0-9]")->mese
  
  dati  %>%
    mutate(mm=mese,yy=2020,dd=day) %>%
    mutate(value=exp(value)-1,fitted=exp(fitted)-1)
  
})->mydf


mydf %>%
  mutate(station_eu_code2=paste0(station_eu_code,stage,trial,mm))->mydf

left_join(mydf,stazioni %>% dplyr::select(-i_surface,-altitudedem,-d_a1,-d_a2,-clc_agricultural))->mydf

ggplot(data=mydf %>% filter(stage=="Val"))+
  geom_point(aes(x=value,y=fitted))+
  geom_abline(a=0,b=1,color="red")+
  facet_wrap(mm~tipo_zona)

mydf%>%
  mutate(error=(value-fitted)^2) %>%
  group_by(station_eu_code,mm,stage,tipo_zona) %>%
  summarise(rmse=sqrt(mean(error,na.rm=TRUE)))->valoriRmse


ggplot(data=valoriRmse,aes(x=mm,y=rmse))+
  geom_boxplot()+
  facet_wrap(~stage)


mydf%>%
  group_by(station_eu_code,mm,stage) %>%
  summarise(corr=cor(value,fitted,use="pairwise.complete.obs"))->valoriCor


ggplot(data=valoriCor,aes(x=mm,y=corr))+
  geom_boxplot()+
  facet_wrap(~stage)


cor(mydf[mydf$stage=="Val",c("fitted","value")],use = "pairwise.complete.obs")
cor(mydf[mydf$stage=="Train",c("fitted","value")],use = "pairwise.complete.obs")


cor(mydf[mydf$tipo_zona=="U" & mydf$stage=="Val",c("fitted","value")],use = "pairwise.complete.obs")
cor(mydf[mydf$tipo_zona=="R" & mydf$stage=="Val",c("fitted","value")],use = "pairwise.complete.obs")
cor(mydf[mydf$tipo_zona=="S" & mydf$stage=="Val",c("fitted","value")],use = "pairwise.complete.obs")

mydf %>%
  dplyr::select(station_eu_code,stage,trial,value,fitted,mm,yy,tipo_zona) %>%
  rename(observed=value) %>%
  tidyr::gather(key="tipo_valori",value="value",-tipo_zona,-yy,-mm,-station_eu_code,-trial,-stage)->gather_mydf

ggplot(data=gather_mydf )+
  geom_boxplot(aes(x=tipo_zona,y=value,fill=tipo_valori))+
  facet_wrap(~stage)

mydf %>%
  filter(stage=="Val") %>%
  mutate(station_eu_code=paste0(station_eu_code,stage,trial,mm))->validazione


pdf("validazione.pdf",width=12,height=8,onefile = TRUE)
purrr::walk(unique(validazione$station_eu_code),.f=function(codex){

  ggplot(data=validazione %>% filter(station_eu_code==codex) )+
    geom_ribbon(aes(x=dd,ymin=fitted025,ymax=fitted975),fill="lightgray")+
    geom_line(aes(x=dd,y=value,group=station_eu_code))+
    geom_line(aes(x=dd,y=fitted,group=station_eu_code),color="red")+
    facet_wrap(~yy)+
    theme_bw()->grafico
  
  print(grafico)

})
dev.off()


sqrt(mean((validazione$fitted-validazione$value)^2,na.rm=TRUE))
mean((validazione$fitted-validazione$value),na.rm=TRUE)

mydf %>%
  filter(stage=="Train") %>%
  mutate(station_eu_code=paste0(station_eu_code,stage,trial,mm))->training


pdf("training.pdf",width=12,height=8,onefile = TRUE)
purrr::walk(unique(training$station_eu_code),.f=function(codex){
  ggplot(data=training %>% filter(station_eu_code==codex) )+
    geom_line(aes(x=dd,y=value,group=station_eu_code))+
    geom_line(aes(x=dd,y=fitted,group=station_eu_code),color="red")+
    facet_wrap(~yy)+
    theme_bw()->grafico
  
  print(grafico)

})
dev.off()


sqrt(mean((training$fitted-training$value)^2,na.rm=TRUE))
mean((training$fitted-training$value),na.rm=TRUE)

