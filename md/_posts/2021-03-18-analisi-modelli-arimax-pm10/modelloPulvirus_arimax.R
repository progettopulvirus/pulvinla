rm(list=objects())
library("tidyverse")
library("rpulvinla") 
library("datiMeteo") 
library("stazioniMonitoraggio") #anagrafica stazioni
library("forecast")
library("lubridate")
options(warn=1)

#future::plan(strategy = "multicore",workers=3)

### rimozione file preesistenti .RDS e .html
system("rm -rf *.RDS")
system("rm -rf *.html")
system("rm -rf *.csv")


#############################
### QUale anno analizzare ? (tra 2013 e 2020)
#############################
ANNO<-2020

#############################
### QUale mesi ? (tra gennaio e maggio)
#############################
MESI<-1:5

#############################
### Inizioe fine lockdown
#############################
inizioL<-"2020-03-09"
fineL<-"2020-05-03"

###############################
### Trasformazione logaritmica della variabile target (concentrazioni)?
###############################
LOGARITMO<-c(TRUE,FALSE)[1] 

###############################
###Quali variabili meteo? 
###############################
METEO<-c("t2m","tp","ptp","sp","wdir","wspeed","pblmax","pblmin","dtr")

###############################
###Usare come predittore il valore dell'inquinante nel giorno precedente? Questa variabile  si chiamera' "pvalue" (previous value). Se LOGARITMO==TRUE pvalue sara' il logaritmo del giorno precedente.
###############################
PREVIOUS<-c(TRUE,FALSE)[2]

#####
# Inserire trend temporale su giorno/settimana/wday (giorni da 1 a 7)/month/weekend?
#####
DAY_TREND<-c(TRUE,FALSE)[2] #rw1 su giorno? ..usare variabile day
WEEK_TREND<-c(TRUE,FALSE)[2] #rw1 su settimana? ... variabile week
WDAY_TREND<-c(TRUE,FALSE)[2] #rw1 su giorno della settimana, lun, mar..dom? ... variabile wday  
MONTH_TREND<-c(TRUE,FALSE)[2] #trend (lineare) sul mese?...variabile mm  
WEEKEND<-c(TRUE,FALSE)[1] #effetto weekend  

#############################
###Regione e inquinante: fissare REGIONE e INQUINANTE
#############################
REGIONI<-PULVIRUS_REGIONI[-2]
purrr::walk(REGIONI,.f=installa_pacchetto_pulvirus)
INQUINANTE<-"pm10" #inquinante su cui lavorare 


#############################
### Inizio programma
#############################

#############################
#crea la formula del modello
#############################
reduce(METEO,.f=str_c,sep="+")->formula.regressori
paste0("value~",formula.regressori)->myformula
as.formula(myformula)->myformula
update(myformula,.~.+lockdown)->myformula

if(PREVIOUS) update(myformula,.~.+pvalue)->myformula
if(DAY_TREND) update(myformula,.~.+day)->myformula
if(WEEK_TREND) update(myformula,.~.+f(week,model="rw2"))->myformula
if(WDAY_TREND) update(myformula,.~.+f(wday,model="rw1",cyclic = FALSE))->myformula
if(MONTH_TREND) update(myformula,.~.+mm)->myformula
if(WEEKEND) update(myformula,.~.+weekend)->myformula

print(myformula)

terms(myformula)->termini
attr(termini,which="term.labels")->VARIABILI #


purrr::map(REGIONI,.f=function(REGIONE){

tryCatch({
  
  caricaDati(pacchetto=REGIONE,inquinante=INQUINANTE)
  
},error=function(e){
  
  message(glue::glue("Verificare: pacchetto {REGIONE} mancante o l'inquinante {INQUINANTE} non disponibile nel pacchetto {REGIONE}"))
  
  NULL
  
})->datiTemp
  
if(is.null(datiTemp)) return()  
  

message(glue::glue("Elaboro regione {REGIONE}"))  
  
datiTemp %>% 
  filter(mm %in% MESI & yy==ANNO)->datiTemp #gennaio..maggio 2020 

stopifnot(nrow(datiTemp)!=0)

#dati ok? verifichiamo
which(is.na(datiTemp$pollutant_fk))->righe
if(length(righe)) stop("pollutant_fk NA???")

#dati meteo standardizzati
meteo_standardizzati[,c("station_eu_code","date","coordx","coordy",METEO)] ->meteo

#rpulvinla::prepara_dati
#aggiunge la variabile banda (per l'SPDE o trend temporale sui giorni), fa il logaritmo della variabile value (var target) e aggiunge la variabile lockdown
prepara_dati(.x = datiTemp,
             previous=PREVIOUS,
             logaritmo = LOGARITMO,
             lockdown = TRUE,
             inizio_lockdown = inizioL,
             fine_lockdown = fineL,
             wday=WDAY_TREND,
             day=DAY_TREND,
             week=WEEK_TREND,
             weekend = WEEKEND)->dati

rm(datiTemp)

#Sia che lavoriamo con il logaritmo dell'inquinante, sia che lavoriamo con la variabile originale, la variabile
#target si chiamera' "value"

if(LOGARITMO){
  
  dati$value<-dati$lvalue
  dati$lvalue<-NULL
  
  if(PREVIOUS){
    dati$pvalue<-dati$lpvalue
    dati$lpvalue<-NULL
  }
  
} #fine LOGARITMO  


#associo dati meteo
left_join(dati,meteo,by = c("station_eu_code", "date"))->dati

#Altre info
unique(dati$station_eu_code)->CODICI
length(CODICI)->numero_stazioni

#############################
#### ARIMAX MODELS
#############################

purrr::partial(.f=round,digits=3)->myround

purrr::imap(CODICI,.f=function(.codice,i){
  
  # Select data and create matrices for estimation/prediction
  dati %>% 
    filter(station_eu_code==.codice)->dati_staz
  
  if(!nrow(dati_staz)) stop(glue::glue("Nessun dato per la stazione {.codice}, impossibile!"))
  
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
  

  print(.codice)

  if(is.null(modts)){
    sink("_error_arimax.csv",append = TRUE)
    cat(glue::glue("{.codice};{REGIONE}\n"))
    sink()
    return()
  }
  # Retrieve outputs
  modts$coef/sqrt(diag(modts$var.coef))->tstat
  myround(2*(1-pnorm(abs(tstat))))->pvalue
    
  tibble(station_eu_code=.codice,
                   variabili=  attr(modts$coef,"names"),
                   coeff = myround(modts$coef),
                   se = myround(sqrt(diag(modts$var.coef))),
                   pvalue=pvalue)->out
  
  # compute residuals and fitted values
  residuals(modts, type="regression")->regerr
  residuals(modts, type="innovation")->armaerr

  if ("intercept" %in% names(modts$coef)){
      predfitted <- as.numeric(regerr + (X_est %*% matrix(modts$coef[colnames(X_est)],ncol=1) + modts$coef["intercept"]) - armaerr )
      predcounterf <- as.numeric(regerr + (X_pred %*% matrix(modts$coef[colnames(X_pred)],ncol=1) + modts$coef["intercept"]) - armaerr)
  } else { #no estimated intercept
      predfitted <- as.numeric(regerr + (X_est %*% matrix(modts$coef[colnames(X_est)],ncol=1) ) - armaerr )
      predcounterf <- as.numeric(regerr + (X_pred %*% matrix(modts$coef[colnames(X_pred)],ncol=1) ) - armaerr )
  }
    
  dati_staz %>% 
    mutate(predfitted=predfitted,predcounterf=predcounterf) %>%
    dplyr::select(station_eu_code,date,value,predfitted,predcounterf)->df_finale
    
  list("dati_staz"=df_finale,"out"=out)

    
})->dati_staz_list

#tolgo i NULL corrispondenti alle stazioni per cui auto.arima non ha funzionato
purrr::compact(dati_staz_list)->dati_staz_list

if(!length(dati_staz_list)){
  message(glue::glue("Nessuna stazione per la regione {REGIONE}"))
  return()
}

purrr::map_dfr(dati_staz_list,"out") %>%
  filter(grepl("lockdown",variabili)) %>%
  dplyr::select(station_eu_code,variabili,pvalue) %>%
  tidyr::spread(key=variabili,value=pvalue)->out_pvalues



purrr::map_dfr(dati_staz_list,"dati_staz")->dati_staz

list(out_pvalues=out_pvalues,dati_staz=dati_staz)

})->listaFinale

#tolgo i NULL corrispondenti alle stazioni per cui auto.arima non ha funzionato
purrr::compact(listaFinale)->listaFinale

if(!length(listaFinale)) stop(glue::glue("Nessuna stazione vaida per l''inquinante {INQUINANTE}"))

purrr::map_dfr(listaFinale,"dati_staz")->dati_staz
saveRDS(dati_staz,glue::glue("_dati_staz_{INQUINANTE}.RDS")) 

purrr::map_dfr(listaFinale,"out_pvalues")->out_pvalues
saveRDS(out_pvalues,glue::glue("_out_pvalues_{INQUINANTE}.RDS")) 

