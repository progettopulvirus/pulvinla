#modello che usa le serie di CAMS standardizzato come predittori: periodo gennaio-maggio 2020
#rh ad aprile correlato con dtr: -0.72, mentre negli altri mesi non correlati. Inn gennaio febbraio e marzo, rh non significativo mentre significativo dtr
#Togliamo rh
#Sp correlata con t2m. Togliamo sp.
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

system("rm -rf *.html")
system("rm -rf *.RDS")
#############################
###Salvare l'output del modello?
#############################
SAVE_OUTPUT<-c(TRUE,FALSE)[1] 

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
METEO<-c("t2m","tp","dtr","wspeed","pblmax","pblmin","pwspeed","rh","nirradiance","sp")

###############################
###Quali variabili spaziali? 
###############################
SPATIAL<-c("altitudedem","d_a2","clc_arable_agri","tipo_zona") #"clc_agricultural","clc_deciduous","clc_evergreen","clc_pasture","clc_shrub","clc_crop"

###############################
###Trasformazione logaritmica della variabile target?
###############################
LOGARITMO<-c(TRUE,FALSE)[1]

#####
# iid su centraline (station_eu_code)?
#####
IID_STAZIONI<-c(TRUE,FALSE)[1] 

#############################
###Plot della mesh?
#############################
DISEGNA_MESH<-c(TRUE,FALSE)[1]

####### i_surface, d_a1,altitudedem le standardizzo su TUTTE le stazioni dell'anagrafica..sono questi i parametri che ho poi usato per standardizzare i corrispettivi rasters
scala<-function(x){
  
  mean(x,na.rm=TRUE)->media
  sd(x,na.rm=TRUE)->deviazione
  
  (x-media)/deviazione
  
}#fine scala
####### 



stazioni$clc_arable_agri<-stazioni$clc_agricultural+stazioni$clc_arable
stazioni$clc_arable_agri<-scala(stazioni$clc_arable_agri)
stazioni$altitudedem<-scala(stazioni$altitudedem)
stazioni$d_a2<-scala(stazioni$d_a2)
#stazioni$d_a1<-scala(stazioni$d_a1)
stazioni$d_a2_altitudedem<-stazioni$d_a2*stazioni$altitudedem
#stazioni$d_a1_altitudedem<-stazioni$d_a1*stazioni$altitudedem
#stazioni$i_surface<-scala(stazioni$i_surface)


purrr::walk(LISTA_MESI,.f=function(MESE){

  purrr::map_dfr(pianura,.f=~caricaDati(pacchetto=.,inquinante=INQUINANTE))->datiTemp

  
  
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
  
  left_join(datiTemp,stazioni[,c("station_eu_code",SPATIAL)])->datiTemp
  
 
  #crea la formula di base con effetti fissi
  reduce(METEO,.f=str_c,sep="+")->formula.regressori
  reduce(SPATIAL,.f=str_c,sep="+")->formula.regressori.spaziali
  paste0("value~",formula.regressori,"+",formula.regressori.spaziali)->myformula
  as.formula(myformula)->myformula
  
  update(myformula,.~.+Intercept-tipo_zona-1)->myformula    

  ###############################
 
  #16 luglio
  list(theta = list(prior="pc.prec", param=c(1,0.1)))->prec_hyper
  if(IID_STAZIONI) update(myformula,.~.+f(station_eu_code,model="iid", hyper = prec_hyper))->myformula
  
  #####
  # Inserire trend temporale (esterno a spde)?
  #####
  DAY_TREND<-c(TRUE,FALSE)[1] #rw1 su giorno? ..usare variabile day
  if(DAY_TREND) update(myformula,.~.+day)->myformula
  
  WEEK_TREND<-c(TRUE,FALSE)[1] #rw1 su settimana? ... variabile week
  #if(WEEK_TREND) update(myformula,.~.+f(week,model="rw2"))->myformula
  
  WDAY_TREND<-c(TRUE,FALSE)[1] #rw1 su giorno della settimana, lun, mar..dom? ... variabile wday  
  #if(WDAY_TREND) update(myformula,.~.+wday)->myformula
  
  MONTH_TREND<-c(TRUE,FALSE)[2] #trend (lineare) sul mese?...variabile mm  
  if(MONTH_TREND) update(myformula,.~.+mm)->myformula
  
  WEEKEND_TREND<-c(TRUE,FALSE)[1]
  if(WEEKEND_TREND) update(myformula,.~.+weekend)->myformula

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

  
  dati$rest<-0
  dati[dati$date %in% giorniFestivi,]$rest<-1

  dati %>%
    filter(yy==2019)->dati19
  
  dati %>%
    filter(yy==2020)->dati20

  right_join(dati19 %>% dplyr::select(-day,-dd,-yy,-pollutant_fk,-regione,-mm),
             dati20 %>% dplyr::select(-mm,-clc_arable_agri,  #-clc_agricultural,-clc_deciduous,-clc_evergreen,-clc_pasture,-clc_shrub,-clc_crop,
                                      -regione,-tipo_zona,-coordx,-coordy,-Intercept,-dd,-yy,-pollutant_fk,-altitudedem,-d_a2,-banda,-day),
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

  bind_cols(datiJoin[,c("station_eu_code","banda","wday","week","coordx","coordy","Intercept","rest","oldvalue.x","oldvalue.y",SPATIAL)],differenze)->dati

  dati %>%
    filter(!is.na(coordx))->dati

  if(WEEKEND_TREND){
    
    dati %>%
      mutate(weekend=ifelse(wday==7,1,0))->dati
    
  }


  min(dati$banda)->banda_first
  
  dati %>%
    mutate(banda=banda-banda_first+1,day=banda)->dati

  #invalidiamo dati anomali
  dati %>%
    mutate(evalue=exp(value)-1) %>%
    mutate(value=ifelse(abs(evalue)>1,NA,value))->dati


  #rpulvinla::numero_giorni
  nrow(dati[!duplicated(dati[,c("week","wday")]),])->n_giorni

  ######################################
  #Per lo stack ho bisogno non delle coordinate del singolo punto, ma delle coordinate per singolo punto e per osservazione
  ######################################
  st_as_sf(dati,coords = c("coordx","coordy"),crs=32632)->sfDati
  
  #rpulvinla::in_km
  in_km(sfDati)->sfDati
  as_tibble(st_coordinates(sfDati))->coordinateOsservazioni  #<- per lo stack 
  names(coordinateOsservazioni)<-c("coordx_km","coordy_km")
  bind_cols(sfDati[,c("station_eu_code","wday","week")],coordinateOsservazioni)->coordinateOsservazioni
  st_geometry(coordinateOsservazioni)<-NULL
  rm(sfDati)

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
  mesh<-inla.mesh.2d(loc =stazioni, max.edge = c(20,25),cutoff=10,min.angle = 30,offset=c(55,30)) 
  #mesh<-inla.mesh.2d(loc =stazioni, max.edge = c(10,25),cutoff=5,min.angle = 30,offset=c(5,10)) # fina
  
  saveRDS(mesh,glue::glue("mesh{MESE}_{REGIONE}.RDS"))


  if(DISEGNA_MESH){
    library("regioniItalia")
    st_union(lombardia,piemonte)->shapeRegione
    st_union(shapeRegione,veneto)->shapeRegione
    st_union(shapeRegione,emiliaromagna)->shapeRegione
    st_union(shapeRegione,toscana)->shapeRegione
    st_union(shapeRegione,valleaosta)->shapeRegione
    st_union(shapeRegione,trentino)->shapeRegione
    st_union(shapeRegione,friuliveneziagiulia)->shapeRegione
    # shapeRegione->lombardia
    # eval(substitute(`::`(regioniItalia,regione),list(regione=REGIONE)))->shapeRegione #oggetto sf della regione (da sovrapporre su grafico mesh), non fondamentale
    #`::`("regioniItalia","toscana")->shapeRegione 
    in_km(shapeRegione)->shapeRegione
    
    png("mesh.png",width=1024,height=1024)
    plot(mesh)
    plot(st_as_sf(as.data.frame(stazioni),coords=c("coordx_km","coordy_km"),crs=st_crs(shapeRegione)),add=TRUE,bg="red",pch=21)
    plot(shapeRegione,add=TRUE,fill="transparent",col="transparent",lwd=2)  
    dev.off()
    
  }#fine DISEGNA_MESH


  ######################## SPDE: Priors & more
  inla.spde2.pcmatern(mesh=mesh,alpha=2,constr=FALSE,prior.range = c(150,0.8),prior.sigma = c(0.5,0.2))->spde
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
  list(prior="pc.cor1",param=c(0.8,0.318))->theta_hyper
  update(myformula,.~.+f(i,model=spde,group = i.group,control.group = list(model="ar1",hyper=list(theta=theta_hyper))))->myformula
  
  print("###############################")
  print("Formula finale modello:")
  print(myformula)
  print("\n\n ATTENZIONE:\n")
  if(LOGARITMO) print("Nella formula value rappresenta il logaritmo della variabile target\n")
  print("###############################")

if(1==0){
  X<-dati[,c("value",SPATIAL[!SPATIAL %in% c("tipo_zona")])]
  lm(value~.-clc_arable,data=X)->zz
  car::vif(zz)
  ppcor::pcor(X %>% dplyr::select(-value))
}
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


  rmarkdown::render("analisi-covariate.Rmd",output_format = "html_document",params = list(regione=REGIONE,inquinante=INQUINANTE))
  # rmarkdown::render("disegnaVariogramma.Rmd",output_format = "html_document",params = list(regione=REGIONE,inquinante=INQUINANTE))
  # rmarkdown::render("analisi-residui-modello.Rmd",output_format = "html_document",params = list(regione=REGIONE,inquinante=INQUINANTE))
  rmarkdown::render("risultati-spde.Rmd",output_format = "html_document",params = list(regione=REGIONE,inquinante=INQUINANTE))
  # rmarkdown::render("codice-modello-prima-versione.Rmd",output_format = "html_document",params = list(regione=REGIONE,inquinante=INQUINANTE))
  
  system(glue::glue("mv analisi-covariate.html analisi-covariate{MESE}.html"))
  system(glue::glue("mv risultati-spde.html risultati-spde{MESE}.html"))


  if(SAVE_OUTPUT){
  
    print("Salvataggio su disco di inla.out")  
    saveRDS(inla.out,glue::glue("result{MESE}_{REGIONE}.RDS"))
    print("Fine scrittura su disco")
    
  }


}) #su MESE
