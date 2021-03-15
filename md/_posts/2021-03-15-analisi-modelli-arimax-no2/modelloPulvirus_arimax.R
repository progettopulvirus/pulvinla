rm(list=objects())
library("tidyverse")
library("sf")
library("sp")
library("rpulvinla") 
library("datiMeteo") 
library("regioniItalia")
library("stazioniMonitoraggio") #anagrafica stazioni
library(forecast)
library(lubridate)
library(tmap)
options(warn=2,error=browser)

set.seed(1)

future::plan(strategy = "multicore",workers=3)

### rimozione file preesistenti .RDS e .html
system("rm -rf *.RDS")
system("rm -rf *.html")
system("rm -rf *.csv")

####
### Calcola dtr
######
meteo$dtr<-meteo$tmax2m-meteo$tmin2m
meteo_standardizzati$dtr<-scale(meteo$dtr)

#############################
###Salvare l'output del modello? Puo' richiedere molto molto tempo, non necessario se i file markdown vengono processati a fine di questo codice con l'oggetto inla.out ancora in memoria
#############################
SAVE_OUTPUT<-c(TRUE,FALSE)[2] 

#############################
### QUale anno analizzare ? (tra 2013 e 2020)
#############################
ANNO<-2020

#############################
### QUale mesi ? (tra gennaio e maggio)
#############################
MESI<-1:5


#############################
###Regione e inquinante: fissare REGIONE e INQUINANTE
#############################
REGIONE<-"lombardia" #viene usato come suffisso per file output e per disegnare la regione nella mesh (solo disegnare, non per costruire la mesh)
INQUINANTE<-"no2" #inquinante su cui lavorare 

caricaDati(pacchetto=REGIONE,inquinante=INQUINANTE)->datiTemp

datiTemp %>% 
  filter(mm %in% MESI & yy==ANNO)->datiTemp #gennaio..maggio 2020 

stopifnot(nrow(datiTemp)!=0)

###############################
###Trasformazione logaritmica della variabile target?
###############################
LOGARITMO<-c(TRUE,FALSE)[1] 

###############################
###Quali variabili meteo? 
###############################
METEO<-c("t2m","tp","ptp","sp","wdir","wspeed","pblmax","pblmin","dtr")

#crea la formula di base con effetti fissi
reduce(METEO,.f=str_c,sep="+")->formula.regressori
paste0("value~",formula.regressori)->myformula
as.formula(myformula)->myformula
#aggiungiamo l'effetto lockdown e inseriamo l'intercetta e togliamo l'intercetta di default
update(myformula,.~.+lockdown)->myformula

###############################
###Usare come predittore il valore dell'inquinante nel giorno precedente? Questa variabile  si chiamera' "pvalue" (previous value). Se LOGARITMO==TRUE pvalue sara' il logaritmo del giorno precedente.
###############################
PREVIOUS<-c(TRUE,FALSE)[2]
if(PREVIOUS) update(myformula,.~.+pvalue)->myformula


#####
# Inserire trend temporale (esterno a spde)?
#####
DAY_TREND<-c(TRUE,FALSE)[2] #rw1 su giorno? ..usare variabile day
if(DAY_TREND) update(myformula,.~.+day)->myformula

WEEK_TREND<-c(TRUE,FALSE)[2] #rw1 su settimana? ... variabile week
if(WEEK_TREND) update(myformula,.~.+f(week,model="rw2"))->myformula

WDAY_TREND<-c(TRUE,FALSE)[2] #rw1 su giorno della settimana, lun, mar..dom? ... variabile wday  
if(WDAY_TREND) update(myformula,.~.+f(wday,model="rw1",cyclic = FALSE))->myformula

MONTH_TREND<-c(TRUE,FALSE)[2] #trend (lineare) sul mese?...variabile mm  
if(MONTH_TREND) update(myformula,.~.+mm)->myformula

WEEKEND<-c(TRUE,FALSE)[1] #effetto weekend  
if(WEEKEND) update(myformula,.~.+weekend)->myformula

########################################
####################### Inizio programma
########################################
#dati meteo standardizzati
meteo_standardizzati[,c("station_eu_code","date","coordx","coordy",METEO)] ->meteo

#dati ok? verifichiamo
which(is.na(datiTemp$pollutant_fk))->righe
if(length(righe)) stop("pollutant_fk NA???")

#rpulvinla::prepara_dati, aggiunge la variabile banda per l'SPDE, fa il logaritmo della variabile value (l'inquinante) e aggiunge la variabile lockdown
prepara_dati(.x = datiTemp,
             previous=PREVIOUS,
             logaritmo = LOGARITMO,
             lockdown = TRUE,
             inizio_lockdown = "2020-03-09",
             fine_lockdown = "2020-05-03",
             wday=WDAY_TREND,
             day=DAY_TREND,
             week=WEEK_TREND,
             weekend = WEEKEND)->dati

#Sia che lavoriamo con il logaritmo dell'inquinante, sia che lavoriamo con la variabile originale, nella creazione dello stack chiamiamo la variabile
#target "value", in modo di avere un'unica interfaccia e toccare il meno possibile il codice

if(LOGARITMO){
  
  dati$value<-dati$lvalue
  dati$lvalue<-NULL
  
  if(PREVIOUS){
    dati$pvalue<-dati$lpvalue
    dati$lpvalue<-NULL
    #dati$pvalue<-inla.group(x=dati$pvalue,n=50,method="quantile")
  }
  
}#fine LOGARITMO  


#associo dati meteo
left_join(dati,meteo)->dati
rm(datiTemp)

######################################
#
######################################
st_as_sf(dati,coords = c("coordx","coordy"),crs=32632)->sfDati
#rpulvinla::in_km
in_km(sfDati)->sfDati
as_tibble(st_coordinates(sfDati))->coordinateOsservazioni  #<- per lo stack 
names(coordinateOsservazioni)<-c("coordx_km","coordy_km")
bind_cols(sfDati[,c("station_eu_code","date")],coordinateOsservazioni)->coordinateOsservazioni
st_geometry(coordinateOsservazioni)<-NULL
rm(sfDati)

#Altre info
unique(dati$station_eu_code)->CODICI
length(CODICI)->numero_stazioni

#rpulvinla::numero_giorni
numero_giorni(dati)->n_giorni

######################## FORMULA MODELLO: la parte random (spde, etc) va aggiunta prima del comando inla()
terms(myformula)->termini
attr(termini,which="term.labels")->VARIABILI #<-mi serve nello stack  


######################## ARIMAX MODELS
system("rm -rf *.csv")

purrr::partial(.f=round,digits=3)->myround

#prepare objects for output
out_pvalues<-tibble(station_eu_code=rep(NA_character_,numero_stazioni),lockdown=rep(NA_real_,numero_stazioni),bandalockdown=rep(NA_real_,numero_stazioni))
out = tibble(station_eu_code=NULL,variabili=NULL,coeff = NULL, se = NULL,pvalue=NULL)

purrr::imap(CODICI,.f=function(.codice,i){
  
  # Select data and create matrices for estimation/prediction
  dati_staz = dati %>% filter(station_eu_code==.codice)

  # save station code and coordiantes for the output
  out_pvalues$station_eu_code[i]<-.codice
  
  #out_pvalues[i,4:5] = dplyr::select(dati_staz, coordx,coordy)

  X = dati_staz %>%  
    dplyr::select(all_of(VARIABILI),banda)
  
  X_est = X %>%
    mutate(bandalockdown = lockdown*banda)

  X_est = as.matrix(X_est)
  
  X_pred = X %>%
    mutate(lockdown=0,bandalockdown=0)

  X_pred = as.matrix(X_pred)
  
  # Estimate the ARIMAX model
  tryCatch({
    
    auto.arima(dati_staz$value, xreg=X_est, d=0)
  
  },error=function(e){
    NULL
  })->modts
  
  #il codice lo assegnamo anche se arimax non funziona
  .codice->>out_pvalues[i,"station_eu_code"]
  
  if(is.null(modts)) return()
   

  # Retrieve outputs
  tstat = modts$coef/sqrt(diag(modts$var.coef))
  pvalue = 2*(1-pnorm(abs(tstat)))
    
  temp<-tibble(station_eu_code=.codice,
                   variabili=  attr(modts$coef,"names"),
                   coeff = myround(modts$coef),
                   se = myround(sqrt(diag(modts$var.coef))),
                   pvalue=myround(pvalue))

  bind_rows(temp,out)->>out #gl

  as.numeric(temp$pvalue[temp$variabili=="lockdown"])->>out_pvalues[i,"lockdown"] 
  as.numeric(temp$pvalue[temp$variabili=="bandalockdown"])->>out_pvalues[i,"bandalockdown"] 

  # compute residuals and fitted values
  regerr = residuals(modts, type="regression")
  armaerr = residuals(modts, type="innovation")
  
  if ("intercept" %in% names(modts$coef)){
      dati_staz$predfitted = regerr + (X_est %*% matrix(modts$coef[colnames(X_est)],ncol=1) + modts$coef["intercept"]) - armaerr 
      dati_staz$predcounterf = regerr + (X_pred %*% matrix(modts$coef[colnames(X_pred)],ncol=1) + modts$coef["intercept"]) - armaerr 
  } else { #no estimated intercept
      dati_staz$predfitted = regerr + (X_est %*% matrix(modts$coef[colnames(X_est)],ncol=1) ) - armaerr 
      dati_staz$predcounterf = regerr + (X_pred %*% matrix(modts$coef[colnames(X_pred)],ncol=1) ) - armaerr 
  }
    
    
    dati_staz
    
})->dati_staz_list

#scrivo out
write_delim(out,"_out.csv",col_names = TRUE,delim=",")

#eliminiamo eventualli NULL da dati_staz_list
purrr::compact(dati_staz_list)->dati_staz_list
if(!length(dati_staz_list)) stop("Nessun risultato ARIMAX disponibile")

saveRDS(dati_staz_list,glue::glue("_dati_staz_list_{INQUINANTE}_{REGIONE}.RDS")) 

out_pvalues %>%
  filter(!is.na(lockdown))->out_pvalues

# Create a table with the significative change in the level and/or trend slope
out_pvalues %>%
  mutate(lockdown_sign= ifelse(out_pvalues$lockdown<0.05, 1, 0)) %>%
  mutate(bandalockdown_sign= ifelse(out_pvalues$bandalockdown<0.05, 1, 0)) %>%
  mutate(sign=case_when(lockdown_sign==0 & bandalockdown_sign == 0 ~ 1,
                        lockdown_sign==1 & bandalockdown_sign == 0 ~ 2,
                        lockdown_sign==0 & bandalockdown_sign == 1 ~ 3,
                        lockdown_sign==1 & bandalockdown_sign == 1 ~ 4)) %>%
                        mutate(sign = factor(sign,levels = 1:4,labels=c("no effect", "level only","slope only", "level & slope")))->out_pvalues

#associo coordinate: 
# st_as_sf(stazioni,coords = c("st_x","st_y"),crs=4326)->sf_stazioni
# st_transform(sf_stazioni,crs=32632)->sf_stazioni
# tibble(as.data.frame(st_coordinates(sf_stazioni)))->sf_stazioni
# names(sf_stazioni)<-c("coordx","coordy")
# sf_stazioni$station_eu_code<-stazioni$station_eu_code
left_join(out_pvalues,stazioni[,c("station_eu_code","st_x","st_y","zona_tipo")])->out_pvalues

saveRDS(out_pvalues,glue::glue("_out_pvalues_{INQUINANTE}_{REGIONE}.RDS")) 


rmarkdown::render("analisi-modelli-arimax-no2.Rmd",
                  output_format = "html_document",
                  params = list(regione=REGIONE,
                                inquinante=INQUINANTE))

