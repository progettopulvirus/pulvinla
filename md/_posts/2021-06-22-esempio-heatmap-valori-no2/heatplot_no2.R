#Possibile grafico per rappresentazione delle concentrazioni di NO2
#14 giugno 2021
#Il grafico si riferisce a tutti i dati? Si. 
#Se si vogliono invalidare valori con variazioni percentuali troppo elevate modificare il codice
rm(list=objects())
library("tidyverse")
library("sf")
library("sp")
library("rpulvinla")
library("stazioniMonitoraggio")
library("seplyr")
library("scico")
library("fontawesome")
library("patchwork")
library("latex2exp")
options(warn=0)

#coppia di anni da confrontare
ANNI<-c(2019,2020)

#su quali mesi lavorare
LISTA_MESI<-1:4
#giorni di festa variabili (pasqua, carnevale)
as.Date(c("2018-02-12","2018-02-13","2018-04-01","2018-04-02",
          "2019-03-04","2019-03-05","2019-04-21","2019-04-22",
          "2020-02-24","2020-02-25","2020-04-12","2020-04-13"))->giorniFestivi

purrr::map(ANNI,.f=function(qualeAnno){
  
  paste0(qualeAnno,c("-01-01","-01-06","-04-25","-05-01","-11-01"))
  
}) %>% unlist()->giorniFestiviFissi

c(as.Date(giorniFestiviFissi,format="%Y-%m-%d"),giorniFestivi)->giorniFestivi



#############################
###Regione e inquinante: fissare REGIONE e INQUINANTE
#############################
REGIONE<-"lombardia" #viene usato come suffisso per file output e per disegnare la regione nella mesh (solo disegnare, non per costruire la mesh)
INQUINANTE<-"no2" #inquinante su cui lavorare 

#gli output riportano il suffisso lombardia ma l'analisi riguarda pianura
c("lombardia","veneto","piemonte","emiliaromagna","toscana","valleaosta","patrento","pabolzano","friuliveneziagiulia")->pianura

#metto lockdown su tutto 2020 perche faccio il confronto fra 2019 vs 2020, quindi non il confronto fra gennaio febraio 2020 vs marzo aprile 2020
INIZIO_LOCKDOWN<-"2020-01-01"
FINE_LOCKDOWN<-"2020-12-31"

###############################
###Quali variabili meteo? 
###############################
METEO<-c("t2m","tp","ptp","dtr","wspeed","pblmax","pblmin","pwspeed","rh","nirradiance","sp")

###############################
###Quali variabili spaziali? 
###############################
SPATIAL<-c("i_surface","altitudedem","d_a2","d_a1","clc_arable","tipo_zona") #"clc_agricultural","clc_deciduous","clc_evergreen","clc_pasture","clc_shrub","clc_crop"

###############################
###Trasformazione logaritmica della variabile target?
###############################
LOGARITMO<-c(TRUE,FALSE)[1]

#####
# iid su centraline (station_eu_code)?
#####
IID_STAZIONI<-c(TRUE,FALSE)[1] 


####### i_surface, d_a1,altitudedem le standardizzo su TUTTE le stazioni dell'anagrafica..sono questi i parametri che ho poi usato per standardizzare i corrispettivi rasters
scala<-function(x){
  
  mean(x,na.rm=TRUE)->media
  sd(x,na.rm=TRUE)->deviazione
  
  (x-media)/deviazione
  
}#fine scala
####### 



stazioni$clc_arable<-stazioni$clc_agricultural+stazioni$clc_arable
stazioni$clc_arable<-scala(stazioni$clc_arable)
stazioni$altitudedem<-scala(stazioni$altitudedem)
stazioni$d_a2<-scala(stazioni$d_a2)
stazioni$d_a1<-scala(stazioni$d_a1)
stazioni$i_surface<-scala(stazioni$i_surface)


purrr::map_dfr(LISTA_MESI,.f=function(MESE){

  purrr::map_dfr(pianura[!grepl("bolzano",pianura)],.f=~caricaDati(pacchetto=.,inquinante=INQUINANTE))->datiTemp
  

  bind_rows(datiTemp,pabolzano::no2)->datiTemp
  
  #prendo solo le stazioni complete 2016-2020
  stazioni %>%
    seplyr::rename_se("inquinante" :=INQUINANTE) %>%
    filter(inquinante %in% c("completa 2016-2020")) %>%
    filter(grepl("^I.+",station_eu_code)) %>%
    filter(tolower(regione) %in% tolower(unique(datiTemp$regione)))->stazioni

  datiTemp %>% 
    filter(station_eu_code %in% stazioni$station_eu_code) %>%
    filter(yy %in% ANNI & mm %in% c(MESE))->datiTemp  

  datiTemp %>%
    mutate(oldvalue=value)->datiTemp
  
  left_join(datiTemp,stazioni[,c("station_eu_code","regione",SPATIAL)])->datiTemp
  
 
  #crea la formula di base con effetti fissi

  #####
  # Inserire trend temporale (esterno a spde)?
  #####
  DAY_TREND<-c(TRUE,FALSE)[1] #rw1 su giorno? ..usare variabile day

  WEEK_TREND<-c(TRUE,FALSE)[1] #rw1 su settimana? ... variabile week

  WDAY_TREND<-c(TRUE,FALSE)[1] #rw1 su giorno della settimana, lun, mar..dom? ... variabile wday  

  WEEKEND_TREND<-c(TRUE,FALSE)[1]

  #dati meteo standardizzati
  datiMeteo::meteo_standardizzati[,c("station_eu_code","date","coordx","coordy",METEO)] ->meteo

  #dati ok? verifichiamo
  which(is.na(datiTemp$pollutant_fk))->righe
  if(length(righe)) stop("pollutant_fk NA???")

  #rpulvinla::prepara_dati, aggiunge la variabile banda per l'SPDE, fa il logaritmo della variabile value (l'inquinante) e aggiunge la variabile lockdown
  prepara_dati(.x = datiTemp,logaritmo = LOGARITMO,lockdown = TRUE,inizio_lockdown = INIZIO_LOCKDOWN,fine_lockdown = FINE_LOCKDOWN,wday=WDAY_TREND,day=DAY_TREND,week=WEEK_TREND,weekend = WEEKEND_TREND)->dati
  rm(datiTemp)
  
  #Sia che lavoriamo con il logaritmo dell'inquinante, sia che lavoriamo con la variabile originale, nella creazione dello stack chiamiamo la variabile
  #target "value", in modo di avere un'unica interfaccia e toccare il meno possibile il codice
  if(LOGARITMO){
    
    dati$value<-dati$lvalue
    dati$lvalue<-NULL
    
  }#fine LOGARITMO  

  #associo dati meteo
  
  left_join(dati,meteo)->dati
  BACKUP<-dati
  
  dati$rest<-0
  dati[dati$date %in% giorniFestivi,]$rest<-1

  dati %>%
    filter(yy==2019)->dati19
  
  dati %>%
    filter(yy==2020)->dati20

  right_join(dati19 %>% dplyr::select(-day,-dd,-yy,-pollutant_fk,-regione,-mm),
             dati20 %>% dplyr::select(-mm,-clc_arable,  #-clc_agricultural,-clc_deciduous,-clc_evergreen,-clc_pasture,-clc_shrub,-clc_crop,
                                      regione,-tipo_zona,-coordx,-coordy,-Intercept,-dd,-yy,-pollutant_fk,-i_surface,-altitudedem,-d_a1,-d_a2,-banda,-day),
            by=c("week"="week","wday"="wday","station_eu_code"="station_eu_code"))->datiJoin

  #cerchiamo di cogliere l'effetto di un giorno festivo nel 2019 ma non festivo nel 2020 o viceversa
  datiJoin %>%
    mutate(rest=ifelse((rest.x+rest.y)==1 ,1,0))->datiJoin

  
  #creo variabili differenze
  purrr::map_dfc(c("value",METEO),.f=function(nomeVar){
  
    tibble(datiJoin[[paste0(nomeVar,".y")]]-datiJoin[[paste0(nomeVar,".x")]])->temp
    
    names(temp)<-nomeVar
    
    temp
    
  })->differenze

  rm(dati)

  bind_cols(datiJoin[,c("regione","station_eu_code","banda","wday","week","coordx","coordy","Intercept","rest","oldvalue.x","oldvalue.y",SPATIAL)],differenze)->dati

  dati %>%
    filter(!is.na(coordx))->dati

  if(WEEKEND_TREND){
    
    dati %>%
      mutate(weekend=ifelse(wday==7,1,0))->dati
    
  }


  # min(dati$banda)->banda_first
  
  dati %>%
    mutate(day=banda)->dati

  #invalidiamo dati anomali

  dati %>%
    mutate(evalue=exp(value)-1)->dati
    #mutate(value=ifelse(abs(evalue)>1,NA,value))->dati


#spaghetti plot
min(dati$week)->min_week


offset<-30*(MESE-1) 


dati %>%
  mutate(mm=MESE) %>%
  rename(`2019`=oldvalue.x,`2020`=oldvalue.y) %>%
  dplyr::select(`2019`,`2020`,station_eu_code,banda,wday,mm,regione,coordx,coordy) %>%
  gather(key="yy",value="no2",-banda,-wday,-mm,-station_eu_code,-regione,-coordx,-coordy) %>%
  mutate(indice=offset+banda)


})->finale

finale %>%
  mutate(no2=ifelse(no2>=50,50,no2))->finale

ggplot(data=finale,aes(x=indice,y=station_eu_code,fill=no2))+
  geom_tile( )+
  scale_fill_scico(palette="lisbon",na.value=NA)+
  facet_wrap(~yy,ncol=2,scales="free_x")+
  theme_bw()

ggplot(data=finale,aes(x=indice,y=station_eu_code,fill=no2))+
  geom_tile( )+
  scale_fill_scico(palette="nuuk",na.value=NA)+
  facet_wrap(~yy,ncol=2,scales="free_x")+
  theme_bw()

ggplot(data=finale,aes(x=indice,y=station_eu_code,fill=no2))+
  geom_tile( )+
  scale_fill_scico(palette="berlin",na.value=NA)+
  facet_wrap(~yy,ncol=2,scales="free_x")+
  theme_bw()

ggplot(data=finale,aes(x=indice,y=station_eu_code,fill=no2))+
  geom_tile( )+
  scale_fill_scico(palette="batlow",na.value=NA)+
  facet_wrap(~yy,ncol=2,scales="free_x")+
  theme_bw()

ggplot(data=finale,aes(x=indice,y=station_eu_code,fill=no2))+
  geom_tile( )+
  scale_fill_scico(palette="lapaz",na.value=NA)+
  facet_wrap(~yy,ncol=2,scales="free_x")+
  theme_bw()

ggplot(data=finale,aes(x=indice,y=station_eu_code,fill=no2))+
  geom_tile( )+
  scale_fill_scico(palette="oleron",na.value=NA)+
  facet_wrap(~yy,ncol=2,scales="free_x")+
  theme_bw()

ggplot(data=finale,aes(x=indice,y=station_eu_code,fill=no2))+
  geom_tile( )+
  scale_fill_scico(palette="vik",na.value=NA)+
  facet_wrap(~yy,ncol=2,scales="free_x")+
  theme_bw()

ggplot(data=finale,aes(x=indice,y=station_eu_code,fill=no2))+
  geom_tile( )+
  scale_fill_scico(palette="lajolla",na.value=NA)+
  facet_wrap(~yy,ncol=2,scales="free_x")+
  theme_bw()

ggplot(data=finale,aes(x=indice,y=station_eu_code,fill=no2))+
  geom_tile( )+
  scale_fill_scico(palette="bamako",na.value=NA)+
  facet_wrap(~yy,ncol=2,scales="free_x")+
  theme_bw()

ggplot(data=finale,aes(x=indice,y=station_eu_code,fill=no2))+
  geom_tile( )+
  scale_fill_scico(palette="imola",na.value=NA)+
  facet_wrap(~yy,ncol=2,scales="free_x")+
  theme_bw()

ggplot(data=finale,aes(x=indice,y=station_eu_code,fill=no2))+
  geom_tile( )+
  scale_fill_scico(palette="cork",na.value=NA)+
  facet_wrap(~yy,ncol=2,scales="free_x")+
  theme_bw()

ggplot(data=finale,aes(x=indice,y=station_eu_code,fill=no2))+
  geom_tile( )+
  scale_fill_scico(palette="broc",na.value=NA)+
  facet_wrap(~yy,ncol=2,scales="free_x")+
  theme_bw()
