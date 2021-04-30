#27 aprile
rm(list=objects())
library("tidyverse")
library("INLA")
library("inlabru")
library("sf")
library("sp")
library("furrr")
library("future")
library("raster")
library("stazioniMonitoraggio")
library("regioniItalia")
options(warn=-2,error=recover)
options(future.globals.maxSize= '+Inf')

numeroSimulazioni<-1000

REGIONE<-"lombardia"

as_Spatial(lombardia)->maschera

plan(multicore,workers=15)

#leggo il layer vuoto
raster("../empty_layer.tif")->emptyLayer
extent(emptyLayer)->estensione

c("weekend","weekendlockdown","day","daylockdown")->nomiTemporali


########################################################################################
#nomiSpatial: variabili solo spaziali da utilizzare per le mappe (non spazio temporali)
########################################################################################
nomiSpatial<-c("i_surface","altitudedem","d_a1")

purrr::map(nomiSpatial,.f=function(covariataSpaziale){
  
  nomeFileRaster<-glue::glue("../{covariataSpaziale}_std.tif")
  raster(nomeFileRaster)->myRaster
  crop(myRaster,estensione)->myRaster
  extend(myRaster,estensione)
  
})->listaSpatial


#listaSpatial contiene i raster delle covariate spaziali standardizzate
names(listaSpatial)<-nomiSpatial

purrr::map(3,.f=function(MESE){
  
  readRDS(glue::glue("dati{MESE}_{REGIONE}.RDS"))->dati
  length(unique(dati$date))->numero_giorni
  
  dati %>%
    distinct(date,.keep_all=TRUE) %>%
    dplyr::select(yy,mm,dd,date,day,weekend,daylockdown,weekendlockdown) %>%
    mutate(wday=lubridate::wday(date,week_start=1),week=lubridate::isoweek(date)) %>%
    arrange(date)->calendario
  
  calendario %>%
    filter(yy==2019)->cal2019
  
  calendario %>%
    filter(yy==2020)->cal2020
  
  #questo ci serve per sapere quale giorno di aprile 2019 va confrontato con april 2020. In particolare: la mappa in posizione day.x va confrontata conla mappa daylockdown.y
  full_join(cal2019,cal2020,by=c("wday","week")) %>%
    filter(!is.na(day.x) & !is.na(daylockdown.y)) %>%
    mutate(weekend=ifelse(wday==7,1,0))->calendarioAllineato 

  calendarioAllineato %>%
    distinct(weekend,week)->elementiDaEstrarre
  
  nrow(elementiDaEstrarre[elementiDaEstrarre$weekend==0,])->numeroFeriali #numero di settimane nel mese
  nrow(elementiDaEstrarre[elementiDaEstrarre$weekend==1,])->numeroFestivi #numero di settimane nel mese  

  #lettura parametri per standardizzare i raster
  #-->read_delim("../parametriPerStandardizzareRaster_5giugno2020/rasters_parametri.csv",delim=";",col_names=TRUE)->parametri

  #nome del file di output in cui salvare i risultati dell'interpolazione
  nomeBrick<-glue::glue("exp_pm10_mese{MESE}_simulazioneMedia.nc")


  #leggo risultati modello
  if(!file.exists(glue::glue("result{MESE}_{REGIONE}.RDS"))) stop(glue::glue("Non trovo il file result{MESE}_{REGIONE}.RDS"))

  readRDS(glue::glue("result{MESE}_{REGIONE}.RDS"))->inla.out

  inla.out$misc$configs$contents->contents
  rm(inla.out)

  effect<-"i" #nome del latent field
  which(contents$tag==effect)->id.effect #id.effect mi da la posizione di "i" in contents

  #dentro length ho la dimensione di ciascun effetto
  contents$length[id.effect]/numero_giorni->lunghezzaEffettoi #corrisponde a mesh$n

  readRDS(glue::glue("mesh{MESE}_{REGIONE}.RDS"))->mesh
  #riporto le coordinate in km
  mesh$loc[,1]<-mesh$loc[,1]*1000
  mesh$loc[,2]<-mesh$loc[,2]*1000 
  stopifnot(mesh$n==lunghezzaEffettoi)

  readRDS(glue::glue("inlaSampleOut_mese{MESE}_{REGIONE}.RDS"))->simulazione.out


  #prendiamo le componenti latent da ciascuno dei 1000 samples: quidentro trovo gli effetti fissi
  #e random per costruire il linear predictor delle simulazioni
  purrr::map(simulazione.out,"latent")->latentComponent
  purrr::map(simulazione.out,"hyperpar")->hyperparComponent
  rm(simulazione.out)


  #La mesh mi serve per proiettare il latent field
  readRDS(glue::glue("iset{MESE}_{REGIONE}.RDS"))->iset

  #simulation.out lista di 1000 elemnti, ogni elemento contiene una matrice latent
  #utilizzando la posizione "id.effect" individuada precedentemente mediante contents
  #so dove trovare gli elementi "i" (latent field) di ciascuno dei 31 giorni del mese
  #a cui sommare i fixed effects7
  furrr::future_map(1:numeroSimulazioni,.f=function(SIM){
  
    #qui dentro trovo i nomi che mi permettono di individuare i vari elementi
    #percostruire il linear predictor
    row.names(latentComponent[[SIM]])->nomiMatrice
    
    #intercetta: e' il layer su cui sommiamo poi tutte le altre componenti
    #Inizializzo spatial layer con l'intercetta
    grep("Intercept",nomiMatrice)->riga
    if(length(riga)!=1) browser()
  
    latentComponent[[SIM]][riga,]->INTERCETTA
    setValues(emptyLayer,INTERCETTA)->spatialIntercetta

    if(!is.null(nomiSpatial)){
    
        #calcolo gli efeftti fissi per le variabili spaziali  
        purrr::map(nomiSpatial,.f=function(covariataSpaziale){
    
          grep(covariataSpaziale,nomiMatrice)->riga
          if(length(riga)!=1) browser()
          latentComponent[[SIM]][riga,]->BETA
          BETA*listaSpatial[[covariataSpaziale]]
        
        })->listaOut
        
        
        purrr::reduce(listaOut,.f=`+`)->sommaCovariateSpaziali
        sommaCovariateSpaziali+spatialIntercetta->sommaCovariateSpaziali

    }else{
      
        spatialIntercetta->sommaCovariateSpaziali
    }    
    
    
    purrr::map(nomiTemporali,.f=function(covariataTemporale){
      
      grep(glue::glue("^{covariataTemporale}:"),nomiMatrice)->riga

      if(length(riga)!=1) browser()
      latentComponent[[SIM]][riga,]
      
    })->listaBetaTemporali
    
    names(listaBetaTemporali)<-nomiTemporali
    
    
    #questa e' la posizione iniziale dell'effetto "i" latent field dentro
    #la componente "latent" in ciascuno dei samples
    contents$start[id.effect]->indice0
  
    purrr::map(1:numero_giorni,.f=function(qualeGiorno){
      
      sommaCovariateSpaziali->griglia

      if((calendario[qualeGiorno,]$weekendlockdown==1) && length(grep("weekendlockdown",nomiTemporali))) listaBetaTemporali$weekendlockdown+griglia->griglia
      if((calendario[qualeGiorno,]$weekend==1) && length(grep("^weekend$",nomiTemporali))) listaBetaTemporali$weekend+griglia->griglia
      if((calendario[qualeGiorno,]$day!=0) && length(grep("^day$",nomiTemporali))) listaBetaTemporali$day*calendario[qualeGiorno,]$day+griglia->griglia
      if((calendario[qualeGiorno,]$daylockdown!=0) && length(grep("^daylockdown$",nomiTemporali))) listaBetaTemporali$daylockdown*(calendario[qualeGiorno,]$daylockdown)+griglia->griglia
      
      
      if(length(grep("Transport",nomiTemporali)) && (calendario[qualeGiorno,]$daylockdown!=0)){
        
        as.Date(glue::glue("2020-{calendario[qualeGiorno,]$mm}-{calendario[qualeGiorno,]$dd}"))->giorno
        #Transport e' unicaper tutta Italia
        dati[which(dati$date==giorno)[1],]$Transport->trasporto
        listaBetaTemporali$Transport*trasporto+griglia->griglia
        
      }
      
      #start mi dice la posizione di partenza dell'effetto in nell'elemento "latent" di ciascuno
      #dei 1000 samples prodotti da inla.posterior.sample
      indice0-1+(1:lunghezzaEffettoi)->indiciEffettoi
    
      #aggiorno indice di partenza
      (indiciEffettoi[lunghezzaEffettoi]+1)->>indice0
  
      #spde (media)
      #questo il latent field del sample che poi vado a riproiettare
      latentComponent[[SIM]][indiciEffettoi,]->campo
      
      inla.mesh.projector(mesh,xlim=c(estensione@xmin,estensione@xmax),ylim=c(estensione@ymin,estensione@ymax),dims = c(218,231))->myproj
      inla.mesh.project(myproj,campo)->campoProj
      raster(list(x=myproj$x,y=myproj$y,z=campoProj))->myraster
      crs(myraster)<-CRS("+init=epsg:32632")
  
      projectRaster(myraster,emptyLayer)->SPDE #spde medio per giorno yymmdd
      rm(myraster)  
      rm(campo)

      crop(extend(SPDE,emptyLayer),emptyLayer)->SPDE

      #questa e' la media della variabile su scala logaritmica (meteoLayer gia include spatialLayer)
      griglia+SPDE->finaleLog 


      #Variance of the Gaussian observations 
      (1/hyperparComponent[[SIM]]["Precision for the Gaussian observations"])->mean_var_GO

 
      #campo medio 
      #exp(finaleLog+0.5*mean_var_GO)+0.1->inquinante
      finaleLog->inquinante
      
      #inquinante[inquinante<1]<-1
      trim(inquinante)
    
    
    })->listaMappe

    purrr::map(.x=1:nrow(calendarioAllineato),.f=function(.riga){
      
      calendarioAllineato[.riga,]$day.x->day.x
      calendarioAllineato[.riga,]$daylockdown.y->daylockdown.y
      calendarioAllineato[.riga,]$wday->giornoSettimana
      
      if(is.na(day.x) | is.na(daylockdown.y)) return()
      
      exp((listaMappe[[daylockdown.y]])-(listaMappe[[day.x]]))-1->out
      mask(out,maschera)
      # stack(listaMappe[c(daylockdown.y,day.x)])->mystack
      # calc(mystack,fun=function(x){(x[1])/x[2]})
      
    })->listaConfronti
    
    brick(listaConfronti)->mybrick
    

    
    purrr::map(1:nrow(elementiDaEstrarre),.f=function(.riga){
      
      which((calendarioAllineato$weekend==elementiDaEstrarre[.riga,]$weekend) & (calendarioAllineato$week==elementiDaEstrarre[.riga,]$week))->qualiRighe
      
      if(length(qualiRighe)!=1){ ##
        
        calc(mybrick[[qualiRighe]],fun = mean,na.rm=TRUE)
        
      }else{ #festivi..domenica...ma se in settimana ho un festivo e una domenica....ricado nel caso sopra
        
        mybrick[[qualiRighe]]
        
      }
      
    })->listaRastersMedi
    

    which(elementiDaEstrarre$weekend==0)->feriali
    brickFeriali<-brick(listaRastersMedi[feriali])
    
    which(elementiDaEstrarre$weekend==1)->festivi
    brickFestivi<-brick(listaRastersMedi[festivi])
    
    list(feriali=brickFeriali,festivi=brickFestivi)
    
  })->listaSimulazioni
  

  purrr::map(listaSimulazioni,"feriali")->listaFeriali

  purrr::map(1:numeroFeriali,.f=function(settimana){
    
    calc(brick(map(listaFeriali,.f=`[[`,settimana)),fun = mean,na.rm=TRUE)

    
  }) %>% brick->feriali_mediati_su_simulazioni

  
  purrr::map(listaSimulazioni,"festivi")->listaFestivi
  
  purrr::map(1:numeroFestivi,.f=function(settimana){
    
    calc(brick(map(listaFestivi,.f=`[[`,settimana)),fun = mean,na.rm=TRUE)
    
    
  }) %>% brick->festivi_mediati_su_simulazioni
  
  writeRaster(feriali_mediati_su_simulazioni,glue::glue("feriali_{MESE}_{REGIONE}.tif"),overwrite=TRUE)
  writeRaster(festivi_mediati_su_simulazioni,glue::glue("festivi_{MESE}_{REGIONE}.tif"),overwrite=TRUE)
  
})#FINE PURRR WALK SU MESE
