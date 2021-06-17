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
options(future.globals.maxSize= 13600*2048^2)

furrr_options(scheduling=7)
numeroSimulazioni<-1000

ANNI<-c(2019,2020)
MESI<-2
REGIONE<-"lombardia"
MODELLO<-"METEO1920"

try(dir.create(glue::glue("mese{MESI}")))

st_union(lombardia,piemonte)->shapeRegione
st_union(shapeRegione,veneto)->shapeRegione
st_union(shapeRegione,emiliaromagna)->shapeRegione
st_union(shapeRegione,toscana)->shapeRegione
st_union(shapeRegione,valleaosta)->shapeRegione
st_union(shapeRegione,trentino)->shapeRegione
st_union(shapeRegione,friuliveneziagiulia)->shapeRegione

st_transform(italia_senza_isole,crs=32632)->italia_senza_isole
st_intersection(shapeRegione,italia_senza_isole)->shapeRegione
as_Spatial(shapeRegione)->maschera
rm(shapeRegione)

plan(multicore,workers=9)

#leggo il layer vuoto
raster("../empty_layer.tif")->emptyLayer
raster::crop(emptyLayer,maschera)->emptyLayer
extent(emptyLayer)->estensione

ncell(emptyLayer)->numero_celle
nrow(emptyLayer)->numero_righe
ncol(emptyLayer)->numero_colonne


c("day","weekend")->nomiTemporali


########################################################################################
#nomiSpatial: variabili solo spaziali da utilizzare per le mappe (non spazio temporali)
########################################################################################
nomiSpatial<-c("altitudedem","d_a2","clc_arable_agri","d_a2_altitudedem")

purrr::map(nomiSpatial,.f=function(covariataSpaziale){
  
  nomeFileRaster<-glue::glue("../{covariataSpaziale}_std.tif")
  raster(nomeFileRaster)->myRaster
  crop(myRaster,estensione)->myRaster
  extend(myRaster,estensione)
  
})->listaSpatial


#listaSpatial contiene i raster delle covariate spaziali standardizzate
names(listaSpatial)<-nomiSpatial

purrr::walk(MESI,.f=function(MESE){
  
  readRDS(glue::glue("result{MESE}_{REGIONE}.RDS"))->inla.out
  
  if(!file.exists(glue::glue("inlaSampleOut_mese{MESE}_{REGIONE}.RDS"))){
    inla.posterior.sample(n=1000,inla.out,num.threads = 3)->simulazione.out
    saveRDS(simulazione.out,glue::glue("inlaSampleOut_mese{MESE}_{REGIONE}.RDS"))
  }else{
    readRDS(glue::glue("inlaSampleOut_mese{MESE}_{REGIONE}.RDS"))->simulazione.out
  }
  
  
  readRDS(glue::glue("dati{MESE}_{REGIONE}.RDS"))->dati
  max(dati$banda)->numero_giorni
  
  dati %>%
    mutate(mm=MESE,dd=day,yy=2020) %>%
    mutate(date=as.Date(glue::glue("{yy}-{mm}-{dd}")))->dati
  
  
  dati %>%
    distinct(date,.keep_all=TRUE) %>%
    dplyr::select(yy,mm,dd,date,weekend,week,wday,rest) %>%
    arrange(date)->calendarioAllineato 
  
  nrow(calendarioAllineato)->nCA
  
  calendarioAllineato %>%
    distinct(weekend,week)->elementiDaEstrarre

  nrow(elementiDaEstrarre)->eDE

  nrow(elementiDaEstrarre[elementiDaEstrarre$weekend==0,])->numeroFeriali #numero di settimane nel mese
  nrow(elementiDaEstrarre[elementiDaEstrarre$weekend==1,])->numeroFestivi #numero di settimane nel mese  

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
    if(length(riga)!=1) stop()
  
    latentComponent[[SIM]][riga,]->INTERCETTA

    if(!is.null(nomiSpatial)){
    
        #calcolo gli efeftti fissi per le variabili spaziali  
        purrr::map(nomiSpatial,.f=function(covariataSpaziale){
    
          grep(glue::glue("^{covariataSpaziale}:.+"),nomiMatrice)->riga
          if(length(riga)!=1) stop()
          latentComponent[[SIM]][riga,]->BETA
          BETA*listaSpatial[[covariataSpaziale]]
        
        })->listaOut
        
        
        purrr::reduce(listaOut,.f=`+`)->sommaCovariateSpaziali
        sommaCovariateSpaziali+INTERCETTA->sommaCovariateSpaziali

    }else{
      
        setValues(emptyLayer,INTERCETTA)->sommaCovariateSpaziali
    }    
    
    
    purrr::map(nomiTemporali,.f=function(covariataTemporale){
      
      grep(glue::glue("^{covariataTemporale}:"),nomiMatrice)->riga

      if(length(riga)!=1) stop()
      latentComponent[[SIM]][riga,]
      
    })->listaBetaTemporali
    
    names(listaBetaTemporali)<-nomiTemporali
    
    
    #questa e' la posizione iniziale dell'effetto "i" latent field dentro
    #la componente "latent" in ciascuno dei samples
    contents$start[id.effect]->indice0

    purrr::map(1:numero_giorni,.f=function(qualeGiorno){
      
      sommaCovariateSpaziali->griglia
      
      if(length(grep("^day$",nomiTemporali))) (listaBetaTemporali$day*calendarioAllineato[qualeGiorno,]$dd)+griglia->griglia
      if(length(grep("^weekend$",nomiTemporali)) & (calendarioAllineato[qualeGiorno,]$weekend==1)) (listaBetaTemporali$weekend*calendarioAllineato[qualeGiorno,]$weekend)+griglia->griglia
      if(length(grep("^rest$",nomiTemporali)) & (calendarioAllineato[qualeGiorno,]$rest==1)) (listaBetaTemporali$rest*calendarioAllineato[qualeGiorno,]$rest)+griglia->griglia
      
      
      #start mi dice la posizione di partenza dell'effetto in nell'elemento "latent" di ciascuno
      #dei 1000 samples prodotti da inla.posterior.sample
      indice0-1+(1:lunghezzaEffettoi)->indiciEffettoi
    
      #aggiorno indice di partenza
      (indiciEffettoi[lunghezzaEffettoi]+1)->>indice0
  
      #spde (media)
      #questo il latent field del sample che poi vado a riproiettare
      inla.mesh.projector(mesh,xlim=c(estensione@xmin,estensione@xmax),ylim=c(estensione@ymin,estensione@ymax),dims = c(218,231))->myproj
      inla.mesh.project(myproj,latentComponent[[SIM]][indiciEffettoi,])->campoProj
      raster(list(x=myproj$x,y=myproj$y,z=campoProj),crs=CRS("+init=epsg:32632"))->myraster

      projectRaster(myraster,emptyLayer)->SPDE #spde medio per giorno yymmdd
      crop(extend(SPDE,emptyLayer),emptyLayer)->SPDE

      #questa e' la media della variabile su scala logaritmica (meteoLayer gia include spatialLayer)
      griglia+SPDE->finale
      
      
      #Precision for station_eu_code
      (1/hyperparComponent[[SIM]]["Precision for station_eu_code"])->mean_var_iid
      rnorm(n=numero_celle,mean=0,sd=mean_var_iid)->vettore_iid
      matrix(vettore_iid,nrow=numero_righe,ncol=numero_colonne,byrow = TRUE)->matrice
      raster(matrice)->raster_iid
      crs(raster_iid)<-crs(emptyLayer)
      extent(raster_iid)<-estensione
      
      finale+raster_iid->finale

      #Variance of the Gaussian observations 
      (1/hyperparComponent[[SIM]]["Precision for the Gaussian observations"])->mean_var_GO

 
      #campo medio 
      (exp(finale+0.5*mean_var_GO))-1
      
 
    })->listaMappe

    raster::mask(brick(listaMappe),maschera)->mybrick

    
    purrr::map(1:eDE,.f=function(.riga){
      
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

  purrr::walk(1:numeroFeriali,.f=function(settimana){
    
    brick(map(listaFeriali,.f=`[[`,settimana))->daScrivere
    writeRaster(daScrivere,glue::glue("./mese{MESI}/feriali_settimana{settimana}_modello{MODELLO}_mese{MESE}_regione{REGIONE}.nc"),format="CDF",overwrite=TRUE)

    
  }) 

  
  purrr::map(listaSimulazioni,"festivi")->listaFestivi
  
  purrr::walk(1:numeroFestivi,.f=function(settimana){
    
    brick(map(listaFestivi,.f=`[[`,settimana))->daScrivere
    writeRaster(daScrivere,glue::glue("./mese{MESI}/festivi_settimana{settimana}_modello{MODELLO}_mese{MESE}_regione{REGIONE}.nc"),format="CDF",overwrite=TRUE)
    
    
  }) 

  

  
})#FINE PURRR WALK SU MESE

setwd(glue::glue("./mese{MESI}"))
system("./cdoSTATSDaily.sh")
