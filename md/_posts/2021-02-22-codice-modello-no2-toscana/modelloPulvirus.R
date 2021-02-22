rm(list=objects())
library("tidyverse")
library("INLA")
library("sf")
library("sp")
#devtools::install_github("progettopulvirus/rpulvinla")
library("rpulvinla")
options(warn=0)

set.seed(1)
inla.setOption(pardiso.license="~/pardiso/licenza.txt")


#############################
###Regione e inquinante: fissare REGIONE e INQUINANTE
#############################
REGIONE<-"toscana" #viene usato come suffisso per file output e per disegnare la regione nella mesh (solo disegnare, non per costruire la mesh)
INQUINANTE<-"pm10" #inquinante su cui lavorare 

caricaDati(pacchetto=REGIONE,inquinante=INQUINANTE)->datiTemp

datiTemp %>% 
  filter(mm %in% seq(1,5) & yy==2020)->datiTemp #gennaio..maggio 2020 

###############################
###Trasformazione logaritmica della variabile target?
###############################
LOGARITMO<-c(TRUE,FALSE)[1] 

###############################
###Quali variabili meteo? 
###############################
METEO<-c("t2m","tp","ptp","sp","wdir","wspeed","pblmax","pblmin","altitudedem")

#crea la formula di base con effetti fissi
reduce(METEO,.f=str_c,sep="+")->formula.regressori
paste0("value~",formula.regressori)->myformula
as.formula(myformula)->myformula
update(myformula,.~.+lockdown)->myformula


###############################
###Inserire le coordinate?
###############################
SCOORDX_KM<-c(TRUE,FALSE)[1] 
SCOORDY_KM<-c(TRUE,FALSE)[1] 

if(SCOORDX_KM) update(myformula,.~.+scoordx_km)->myformula
if(SCOORDY_KM) update(myformula,.~.+scoordy_km)->myformula


#####
# iid su centraline (station_eu_code)?
#####
IID_STAZIONI<-c(TRUE,FALSE)[1] 
if(IID_STAZIONI) update(myformula,.~.+f(station_eu_code,model="iid"))->myformula

#####
# Inserire trend temporale (esterno a spde)?
#####
DAY_TREND<-c(TRUE,FALSE)[2] #rw1 su giorno? ..usare variabile day
if(DAY_TREND) update(myformula,.~.+f(day,model="rw1"))->myformula

WEEK_TREND<-c(TRUE,FALSE)[1] #rw1 su settimana? ... variabile week
if(WEEK_TREND) update(myformula,.~.+f(week,model="rw1"))->myformula

WDAY_TREND<-c(TRUE,FALSE)[1] #rw1 su giorno della settimana, lun, mar..dom? ... variabile wday  
if(WDAY_TREND) update(myformula,.~.+f(wday,model="rw1"))->myformula

MONTH_TREND<-c(TRUE,FALSE)[2] #trend (lineare) sul mese?...variabile mm  
if(MONTH_TREND) update(myformula,.~.+mm)->myformula

#############################
###Salvare l'output del modello?
#############################
SAVE_OUTPUT<-c(TRUE,FALSE)[1] 

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
prepara_dati(.x = datiTemp,logaritmo = LOGARITMO,lockdown = TRUE,wday=WDAY_TREND,day=DAY_TREND,week=WEEK_TREND)->dati

#Sia che lavoriamo con il logaritmo dell'inquinante, sia che lavoriamo con la variabile originale, nella creazione dello stack chiamiamo la variabile
#target "value", in modo di avere un'unica interfaccia e toccare il meno possibile il codice

if(LOGARITMO){
  
  dati$value<-dati$lvalue
  dati$lvalue<-NULL
  
}#fine LOGARITMO  



#associo dati meteo
left_join(dati,meteo)->dati
rm(datiTemp)


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

saveRDS(dati,glue::glue("dati_{REGIONE}.RDS"))

#Altre info
unique(dati$station_eu_code)->CODICI

#rpulvinla::numero_giorni
numero_giorni(dati)->n_giorni

####Priors: AR1
list(prior="pc.cor0",param=c(0.7,0.1),fixed=FALSE)->theta_hyper #<-------------- questa la prior usata ne modello per ar1
####

#16 luglio
list(theta = list(prior="pc.prec", param=c(1,0.1)))->prec_hyper


######################## FORMULA MODELLO: la parte random (spde, etc) va aggiunta prima del comando inla()
terms(myformula)->termini
attr(termini,which="term.labels")->VARIABILI #<-mi serve nello stack  
########################

#Utilizzo coordinateOsservazioni per costruire la mesh: in questo caso ho bisogno di un punto per stazione (non per stazione E per osservazione)
coordinateOsservazioni[!duplicated(coordinateOsservazioni$station_eu_code),c("coordx_km","coordy_km")]->stazioni
as.matrix(stazioni)->stazioni

### Quale mesh?
mesh<-inla.mesh.2d(loc =stazioni, max.edge = c(12.5,30),cutoff=5,min.angle = 30,offset=c(50,25))
saveRDS(mesh,glue::glue("mesh_{REGIONE}.RDS"))


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
inla.spde2.pcmatern(mesh=mesh,alpha=2,constr=FALSE,prior.range = c(70,0.8),prior.sigma = c(0.7,0.2))->spde
saveRDS(spde,glue::glue("spde_{REGIONE}.RDS"))

inla.spde.make.index(name="i",n.spde=spde$n.spde,n.group = n_giorni)->iset
saveRDS(iset,glue::glue("iset_{REGIONE}.RDS"))


#training
inla.spde.make.A(mesh=mesh,loc=as.matrix(coordinateOsservazioni[,c("coordx_km","coordy_km")]),group =dati$banda,n.spde=spde$n.spde,n.group =n_giorni )->A.training
#La variabile target si chiama value, che sia il logaritmo o no
str_remove(str_remove(attr(termini,"term.labels"),pattern=",.+$"),pattern="^f\\(")->EFFETTI

inla.stack(data=list(value=dati$value),A=list(A.training,1),effects=list(iset,dati[EFFETTI]),tag="training")->stack.training
saveRDS(stack.training,glue::glue("stack.training_{REGIONE}.RDS"))

########################
#Random effects
########################
update(myformula,.~.+f(i,model=spde,group = i.group,control.group = list(model="ar1",hyper=list(theta1=list(prior="pc.prec", param=c(1,0.1)),rho=theta_hyper))))->myformula

print("###############################")
print("Formula finale modello:")
print(myformula)
print("\n\n ATTENZIONE:\n")
if(LOGARITMO) print("Nella formula value rappresenta il logaritmo della variabile target\n")
print("###############################")


######################## INLA Vai!
inla(myformula,
     data=inla.stack.data(stack.training,spde=spde),
     family ="gaussian",
     verbose=TRUE,
     control.compute = list(openmp.strategy="pardiso.parallel",cpo=TRUE,waic=TRUE,dic=TRUE,config=TRUE),
     control.fixed = list(prec.intercept = 0.001, prec=1,mean.intercept=0),
     control.predictor =list(A=inla.stack.A(stack.training),compute=TRUE) )->inla.out


rmarkdown::render("analisi-covariate.Rmd",output_format = "html_document",params = list(regione=REGIONE,inquinante=INQUINANTE))
rmarkdown::render("analisi-residui-modello.Rmd",output_format = "html_document",params = list(regione=REGIONE,inquinante=INQUINANTE))
rmarkdown::render("risultati-spde.Rmd",output_format = "html_document",params = list(regione=REGIONE,inquinante=INQUINANTE))
rmarkdown::render("codice-modello-prima-versione.Rmd",output_format = "html_document",params = list(regione=REGIONE,inquinante=INQUINANTE))

if(SAVE_OUTPUT){

  print("Salvataggio su disco di inla.out")  
  saveRDS(inla.out,glue::glue("result_{REGIONE}.RDS"))
  print("Fine scrittura su disco")
  
}
