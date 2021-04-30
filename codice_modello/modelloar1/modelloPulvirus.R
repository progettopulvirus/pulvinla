#modello che usa le serie di CAMS standardizzato come predittori: periodo gennaio-maggio 2020
rm(list=objects())
library("tidyverse")
library("INLA")
library("sf")
library("sp")
library("rpulvinla")
library("stazioniMonitoraggio")
library("seplyr")
options(warn=0)

set.seed(1)
inla.setOption(pardiso.license="~/pardiso/licenza.txt")

system("rm -rf *.html")
#system("rm -rf *.RDS")
#############################
###Salvare l'output del modello?
#############################
SAVE_OUTPUT<-c(TRUE,FALSE)[1] 
readRDS("../cams_standardizzato_gennaio_aprile_2019_2020.RDS") %>%
  gather(key="station_eu_code",value="value",-yy,-mm,-dd)->std_cams

#############################
###Regione e inquinante: fissare REGIONE e INQUINANTE
#############################
REGIONE<-"lombardia" #viene usato come suffisso per file output e per disegnare la regione nella mesh (solo disegnare, non per costruire la mesh)
INQUINANTE<-"no2" #inquinante su cui lavorare 

#metto lockdown su tutto 2020 perche faccio il confronto fra 2019 vs 2020, quindi non il confronto fra gennaio febraio 2020 vs marzo aprile 2020
INIZIO_LOCKDOWN<-"2020-01-01"
FINE_LOCKDOWN<-"2020-04-30"


####### i_surface, d_a1,altitudedem le standardizzo su TUTTE le stazioni dell'anagrafica..sono questi i parametri che ho poi usato per standardizzare i corrispettivi rasters
scala<-function(x){
  
  mean(x,na.rm=TRUE)->media
  sd(x,na.rm=TRUE)->deviazione
  
  (x-media)/deviazione
  
}

stazioni$i_surface<-scala(stazioni$i_surface)
stazioni$altitudedem<-scala(stazioni$altitudedem)
stazioni$d_a1<-scala(stazioni$d_a1)


purrr::walk(1:4,.f=function(MESE){

caricaDati(pacchetto=REGIONE,inquinante=INQUINANTE)->datiTemp

#prendo solo le stazioni complete 2016-2020
stazioni %>%
  seplyr::rename_se("inquinante" :=INQUINANTE) %>%
  filter(inquinante %in% c("completa 2016-2020")) %>%
  filter(grepl("^I.+",station_eu_code)) %>%
  filter(grepl(toupper(REGIONE),regione))->stazioni

datiTemp %>% 
  filter(station_eu_code %in% stazioni$station_eu_code) %>%
  filter(yy %in% c(2019,2020) & mm %in% c(MESE))->datiTemp #gennaio..maggio 2020 
 
#print(skimr::skim(datiTemp[,c("yy","mm","dd","date")]))

std_cams %>%
  rename(cams=value)->std_cams

left_join(datiTemp,std_cams)->datiTemp
left_join(datiTemp,stazioni[,c("station_eu_code","d_a1","i_surface","popolazione","tipo_stazione","altitudedem")])->datiTemp

###############################
###Trasformazione logaritmica della variabile target?
###############################
LOGARITMO<-c(TRUE,FALSE)[1] 

###############################
###Quali variabili meteo? 
###############################
METEO<-c("t2m","tp","ptp","dtr","wspeed","pbl12","pbl00","pwspeed","rh")

#crea la formula di base con effetti fissi
reduce(METEO,.f=str_c,sep="+")->formula.regressori
paste0("value~",formula.regressori)->myformula
as.formula(myformula)->myformula
update(myformula,.~.+i_surface+Intercept+d_a1+altitudedem-1)->myformula
###############################

#####
# iid su centraline (station_eu_code)?
#####
IID_STAZIONI<-c(TRUE,FALSE)[1] 

#16 luglio
list(theta = list(prior="pc.prec", param=c(1,0.1)))->prec_hyper
if(IID_STAZIONI) update(myformula,.~.+f(station_eu_code,model="iid", hyper = prec_hyper))->myformula

#####
# Inserire trend temporale (esterno a spde)?
#####
DAY_TREND<-c(TRUE,FALSE)[1] #rw1 su giorno? ..usare variabile day
if(DAY_TREND) update(myformula,.~.+day)->myformula

WEEK_TREND<-c(TRUE,FALSE)[2] #rw1 su settimana? ... variabile week
if(WEEK_TREND) update(myformula,.~.+f(week,model="rw2"))->myformula

WDAY_TREND<-c(TRUE,FALSE)[1] #rw1 su giorno della settimana, lun, mar..dom? ... variabile wday  
#if(WDAY_TREND) update(myformula,.~.+wday)->myformula

MONTH_TREND<-c(TRUE,FALSE)[2] #trend (lineare) sul mese?...variabile mm  
if(MONTH_TREND) update(myformula,.~.+mm)->myformula

WEEKEND_TREND<-c(TRUE,FALSE)[1]
if(WEEKEND_TREND) update(myformula,.~.+weekend)->myformula

#############################
###Plot della mesh?
#############################
DISEGNA_MESH<-c(TRUE,FALSE)[1]

#dati meteo standardizzati
datiMeteo::meteo_standardizzati[,c("station_eu_code","date","coordx","coordy",METEO)] ->meteo

#dati ok? verifichiamo
which(is.na(datiTemp$pollutant_fk))->righe
if(length(righe)) stop("pollutant_fk NA???")

#rpulvinla::prepara_dati, aggiunge la variabile banda per l'SPDE, fa il logaritmo della variabile value (l'inquinante) e aggiunge la variabile lockdown
prepara_dati(.x = datiTemp,logaritmo = LOGARITMO,lockdown = TRUE,inizio_lockdown = INIZIO_LOCKDOWN,fine_lockdown = FINE_LOCKDOWN,wday=WDAY_TREND,day=DAY_TREND,week=WEEK_TREND,weekend = WEEKEND_TREND)->dati
rm(datiTemp)

dati %>%
  mutate(tipo_stazione=as.factor(tipo_stazione)) %>%
  mutate(daylockdown=day*lockdown) %>%
  mutate(nolockdown=ifelse(lockdown==0,1,0)) %>%
  mutate(day=day*nolockdown) %>%
  mutate(weekendlockdown=ifelse(yy==2019,0,weekend)) %>%
  mutate(weekend=ifelse(yy==2020,0,weekend))->dati

update(myformula,.~.+daylockdown+weekendlockdown)->myformula

read_delim("emissioni.csv",delim=";",col_names=TRUE)->emissioni
caret::preProcess(emissioni[,!names(emissioni) %in% c("yy","mm","dd")],method=c("center","scale"))->proEm
predict(proEm,emissioni[,!names(emissioni) %in% c("yy","mm","dd")])->std_em
bind_cols(emissioni[,c("yy","mm","dd")],std_em)->std_em

#update(myformula,.~.+Transport)->myformula

#aggiungo giorni di festa
as.Date(c("2019-05-01","2020-05-01","2019-01-01","2020-01-01","2019-01-06","2020-01-06","2019-04-21","2019-04-22","2020-01-12","2020-04-13"))->giorniFestivi

if(WEEKEND_TREND){
dati %>%
  mutate(weekend=case_when(date %in% giorniFestivi~1,
         TRUE~weekend))->dati
  
}

if(WDAY_TREND){
  dati %>%
    mutate(wday=case_when(date %in% giorniFestivi~7,
                          TRUE~wday))->dati
  
  dati %>%
    mutate(wday=as.factor(wday))->dati
}

#Sia che lavoriamo con il logaritmo dell'inquinante, sia che lavoriamo con la variabile originale, nella creazione dello stack chiamiamo la variabile
#target "value", in modo di avere un'unica interfaccia e toccare il meno possibile il codice
if(LOGARITMO){
  
  dati$value<-dati$lvalue
  dati$lvalue<-NULL
  
}#fine LOGARITMO  

#associo dati meteo
left_join(dati,meteo)->dati
left_join(dati,std_em)->dati


######################################
#Per lo stack ho bisogno non delle coordinate del singolo punto, ma delle coordinate per singolo punto e per osservazione
######################################
st_as_sf(dati,coords = c("coordx","coordy"),crs=32632)->sfDati

#rpulvinla::in_km
in_km(sfDati)->sfDati
as_tibble(st_coordinates(sfDati))->coordinateOsservazioni  #<- per lo stack 
names(coordinateOsservazioni)<-c("coordx_km","coordy_km")
bind_cols(sfDati[,c("station_eu_code","date")],coordinateOsservazioni)->coordinateOsservazioni
st_geometry(coordinateOsservazioni)<-NULL
rm(sfDati)


##Coordinate standardizzate, se usate come predittori
coordinateOsservazioni %>%
  mutate(scoordx_km=scale(coordx_km),
         scoordy_km=scale(coordy_km))->coordinateOsservazioni

left_join(dati, coordinateOsservazioni %>%dplyr::select(station_eu_code,date,matches("^scoord[xy]_km")))->dati

#rpulvinla::numero_giorni
dati %>%
  distinct(yy,mm,dd) %>%
  count() %>%
  .$n->n_giorni

####Priors: AR1
list(prior="pc.cor0",param=c(0.7,0.1),fixed=FALSE)->theta_hyper #<-------------- questa la prior usata ne modello per ar1
####


######################## FORMULA MODELLO: la parte random (spde, etc) va aggiunta prima del comando inla()
terms(myformula)->termini
attr(termini,which="term.labels")->VARIABILI #<-mi serve nello stack  
########################

#Utilizzo coordinateOsservazioni per costruire la mesh: in questo caso ho bisogno di un punto per stazione (non per stazione E per osservazione)
coordinateOsservazioni[!duplicated(coordinateOsservazioni$station_eu_code),c("coordx_km","coordy_km")]->stazioni
as.matrix(stazioni)->stazioni

### Quale mesh?
#mesh<-inla.mesh.2d(loc =stazioni, max.edge = c(20,25),cutoff=10,min.angle = 30,offset=c(20,25)) # max.edge = c(12.5,30)
mesh<-inla.mesh.2d(loc =stazioni, max.edge = c(10,25),cutoff=5,min.angle = 30,offset=c(5,10)) # max.edge = c(12.5,30)

saveRDS(mesh,glue::glue("mesh{MESE}_{REGIONE}.RDS"))


if(DISEGNA_MESH){
  
  
  eval(substitute(`::`(regioniItalia,regione),list(regione=REGIONE)))->shapeRegione #oggetto sf della regione (da sovrapporre su grafico mesh), non fondamentale
  #`::`("regioniItalia","toscana")->shapeRegione 
  in_km(shapeRegione)->shapeRegione
  
  png("mesh.png",width=1024,height=1024)
  plot(mesh)
  plot(st_as_sf(as.data.frame(stazioni),coords=c("coordx_km","coordy_km"),crs=st_crs(shapeRegione)),add=TRUE,bg="red",pch=21)
  plot(shapeRegione,add=TRUE,fill="transparent",col="transparent",lwd=2)  
  dev.off()
  
}#fine DISEGNA_MESH


######################## SPDE: Priors & more
inla.spde2.pcmatern(mesh=mesh,alpha=2,constr=FALSE,prior.range = c(150,0.9),prior.sigma = c(20,0.8))->spde
saveRDS(spde,glue::glue("spde{MESE}_{REGIONE}.RDS"))

inla.spde.make.index(name="i",n.spde=spde$n.spde,n.group = n_giorni)->iset
saveRDS(iset,glue::glue("iset{MESE}_{REGIONE}.RDS"))



#training
inla.spde.make.A(mesh=mesh,loc=as.matrix(coordinateOsservazioni[,c("coordx_km","coordy_km")]),group =dati$banda,n.spde=spde$n.spde,n.group =n_giorni )->A.training


#La variabile target si chiama value, che sia il logaritmo o no
str_remove(str_remove(attr(termini,"term.labels"),pattern=",.+$"),pattern="^f\\(")->EFFETTI

inla.stack(data=list(value=dati$value),A=list(A.training,1),effects=list(iset,dati[c(EFFETTI)]),tag="training")->stack.training
saveRDS(stack.training,glue::glue("stack.training{MESE}_{REGIONE}.RDS"))

stack.training->mystack
saveRDS(mystack,glue::glue("mystack{MESE}_{REGIONE}.RDS"))

########################
#Random effects
########################
update(myformula,.~.+f(i,model=spde,group = i.group,control.group = list(model="ar1")))->myformula

print("###############################")
print("Formula finale modello:")
print(myformula)
print("\n\n ATTENZIONE:\n")
if(LOGARITMO) print("Nella formula value rappresenta il logaritmo della variabile target\n")
print("###############################")

######################## INLA Vai!
inla(myformula,
     data=inla.stack.data(mystack,spde=spde),
     family ="gaussian",
     verbose=TRUE,
     control.inla=list(int.strategy="eb"),
     control.compute = list(openmp.strategy="pardiso.parallel",cpo=F,waic=F,dic=F,config=TRUE),
     control.fixed = list(prec.intercept = 1, prec=0.01,mean.intercept=0),
     control.predictor =list(A=inla.stack.A(mystack),compute=TRUE) )->inla.out


#calcolo residui su scala logaritmica
mystack$data$index$training->righe
inla.out$summary.fitted.values$mean[righe]->dati$fitted
inla.out$summary.fitted.values$`0.025quant`[righe]->dati$fitted025
inla.out$summary.fitted.values$`0.975quant`[righe]->dati$fitted975

dati %>% 
  mutate(residui=value-fitted)->dati
saveRDS(dati,glue::glue("dati{MESE}_{REGIONE}.RDS")) 

rm(righe)

#calcolo residui su scala logaritmica
# mystack$data$index$validation->righe
# inla.out$summary.fitted.values$mean[righe]->datiValidazione$fitted
# inla.out$summary.fitted.values$`0.025quant`[righe]->datiValidazione$fitted025
# inla.out$summary.fitted.values$`0.975quant`[righe]->datiValidazione$fitted975
# 
# datiValidazione %>% 
#   mutate(residui=value-fitted)->dati
# saveRDS(dati,glue::glue("datiValidazione_{REGIONE}.RDS")) 

rmarkdown::render("analisi-covariate.Rmd",output_format = "html_document",params = list(regione=REGIONE,inquinante=INQUINANTE))
# rmarkdown::render("disegnaVariogramma.Rmd",output_format = "html_document",params = list(regione=REGIONE,inquinante=INQUINANTE))
# rmarkdown::render("analisi-residui-modello.Rmd",output_format = "html_document",params = list(regione=REGIONE,inquinante=INQUINANTE))
# rmarkdown::render("risultati-spde.Rmd",output_format = "html_document",params = list(regione=REGIONE,inquinante=INQUINANTE))
# rmarkdown::render("codice-modello-prima-versione.Rmd",output_format = "html_document",params = list(regione=REGIONE,inquinante=INQUINANTE))

system(glue::glue("mv analisi-covariate.html analisi-covariate{MESE}.html"))

if(SAVE_OUTPUT){

  print("Salvataggio su disco di inla.out")  
  saveRDS(inla.out,glue::glue("result{MESE}_{REGIONE}.RDS"))
  print("Fine scrittura su disco")
  
}


}) #su MESE
