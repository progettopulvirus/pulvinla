
rm(list=objects())
library("tidyverse")
library("lombardia")
library("datiMeteo")
library("stazioniMonitoraggio")
#library("regioniItalia") # shapefiles  da ISTAT
library("sf")
library("scico")
library("janitor")
library("GGally")
library("gridExtra")
library("patchwork")
#library("ggspatial")
#library("data.tree") da sviluppare
#library("raster")


# PM10
PARAMETRO<-c("pm10","no2")[1] 
assign("inquinante",eval(parse(text=PARAMETRO)))
#compare::compareIdentical(inquinante,no2)

#bug da correggere: alcune serie hanno record duplicati, eliminare righe con pollutant_fk==NA
inquinante %>%
  filter(!is.na(pollutant_fk))->inquinante

stopifnot(nrow(inquinante)!=0)

#codici stazioni dell'inquinante per tutta la regione: sono i codici delle stazioni che hanno dati (altrimenti non comparirebbero
#nel pacchetto R).
unique(inquinante$station_eu_code)->codiciStazioni

#load("brescia.rda")

left_join(inquinante,meteo,by=c("station_eu_code"="station_eu_code","date"="date"))->dati

##creo un oggetto spaziale
#stazioni[,c("station_eu_code","nome_stazione","st_x","st_y","altitudine")]->puntiStazione
#st_as_sf(puntiStazione ,coords = c("st_x","st_y"),crs=4326)->puntiStazione
#st_transform(puntiStazione,crs = 32632)->puntiStazioneUTM

##stazioni che appartengono a "Area" (regione, provincia..)
#lombardia->Area
##st_intersection(puntiStazioneUTM,Area)->stazioniArea

# mgv
stazioniArea <- stazioni
# \mgv
#stazioni che appartengono ad Area e con dati in "inquinante"
stazioniArea[stazioniArea$station_eu_code %in% codiciStazioni,]->stazioniArea


#' quante stazioni
length(unique(stazioniArea$station_eu_code))

#raster("dem.tif")->dem
#raster::crop(dem,lombardia)->demArea
# 
# theme_set(theme_bw()+theme(panel.grid.minor = element_blank(),
#                            panel.grid.major.x = element_blank()))
# 
# ggplot()+
#   ggspatial::layer_spatial(data=demArea)+
#   ggspatial::geom_sf(data=Area,fill="transparent",colour="black")+
#   geom_sf(data=stazioniArea,aes(fill=altitudine),pch=21)+
#   scale_fill_scico(palette="bamako")->mappa
# 

#print(mappa)

#' numero stazioni
nrow(stazioniArea)


#selezioniamo dati con codici in "stazioniArea" (stazioni che cadono in Area e con dati)
#selezioniamo serie nel periodo gennaio-giugno
dati %>%
  filter(station_eu_code %in% stazioniArea$station_eu_code) %>%
  filter(mm %in% 1:6)->subDati

left_join(subDati,stazioniMonitoraggio::stazioni[,c("station_eu_code","nome_stazione","tipo_zona")]) %>%
  mutate(tipo_zona=as.factor(tipo_zona))->subDati

# skimr simile a summary
print(skimr::skim(subDati[,c("yy","mm","dd","value","tipo_zona")] ))

print(skimr::skim(subDati %>% dplyr::select(value,tipo_zona) %>% group_by(tipo_zona)))


ggplot(data=subDati)+
  geom_boxplot(aes(x=station_eu_code,y=value))+
  geom_hline(yintercept = 50,colour="red")+
  facet_wrap(~tipo_zona,scales = "free_x",nrow = 3)

# ----
# MGV correlatione variabili meteo vs inquinante
# ----
# studio della correlazione tra variabili
# subdati
# AQ vs Meteo
# scatterpolot etc etc..
# ----
# Studio correlazioni
# --
subMeteoNames <- c("value" , "t2m" , "tmin2m" ,"tmax2m","tp", "ptp", "rh" , "sp",
                   "nirradiance", "pbl00","pbl12", "pblmin","pblmax","wdir", "wspeed","pwspeed" )
# ----
# Correlazioni complessive
# ----

 
gg <-ggcorr(subDati[,subMeteoNames], palette = "RdBu", label = TRUE) +
  ggplot2::labs(title = paste0("Correlazione per le variabili meteo principali per gli anni 2013-2020 e ",PARAMETRO, " (value)"))


# ggsave(plot = gg, filename = paste0(wrkdir,PARAMETRO,"_correlazione_2013_2020.png"), dpi = 300,  
#        units = "in",width = 13, height = 10) 
# ----
# Correlazioni per anno
# ----
purrr::map(unique(subDati$yy), .f=function(.x){

  dfplot <- subDati[,c("yy",subMeteoNames)] %>% filter(yy == .x)
  
  ggcorr(
        dfplot[,-1],
        geom = "circle",
        max_size = 6,
        size = 3,
        hjust = 0.75,
        angle = -45,
        palette = "PuOr") + ggplot2::labs(title = .x)
  
})->plots

gg <-do.call(grid.arrange,plots)
gg

# ggsave(plot = gg, filename = paste0(wrkdir,PARAMETRO,"_correlazione_per_anno.png"), dpi = 300,  
#        units = "in",width = 13, height = 10) 

# ----
# correlazioni per tipo di stazione_zona
# ----
purrr::map(unique(subDati$tipo_zona), .f=function(.x){

  dfplot <- subDati[,c("tipo_zona",subMeteoNames)] %>% filter(tipo_zona==.x)
  
  ggcorr(
    dfplot,
    geom = "circle",
    max_size = 6,
    size = 3,
    hjust = 0.75,
    angle = -45,
    palette = "PuOr") + ggplot2::labs(title = .x)

})->plots
gg <-purrr::reduce(plots,.f=`+`)
gg
# ggsave(plot = gg, filename = paste0(wrkdir,PARAMETRO,"_correlazione_per_zona.png"), dpi = 300,  
#        units = "in",width = 13, height = 10) 

# ----
# Studio scatterplots & correlazioni
# --
subMeteoNames <- c("yy" ,"tipo_zona","value" , "t2m" , "tmin2m" ,"tmax2m","tp", "ptp", "rh" , "sp",
                   "nirradiance", "pbl00","pbl12", "pblmin","pblmax","wdir", "wspeed","pwspeed" )

# scatterplots temperature-radiazione
nomecolns <- c("value" , "t2m" , "tmin2m" ,"tmax2m","nirradiance")

print(skimr::skim(subDati[,c("tipo_zona",nomecolns)] %>% group_by(tipo_zona) ))
gg <- ggpairs(subDati[,subMeteoNames], columns = nomecolns, ggplot2::aes(colour=tipo_zona, alpha = 0.25))      
gg
# ggsave(plot = gg, filename = paste0(wrkdir,PARAMETRO,"_Temp_scatterplot.png"), dpi = 300,  
#        units = "in",width = 13, height = 10) 

# scatterplots temperature-radiazione
nomecolns <- c("value" , "tp", "ptp", "rh" , "sp")
print(skimr::skim(subDati[,c("tipo_zona",nomecolns)] %>% group_by(tipo_zona) ))
gg <- ggpairs(subDati[,subMeteoNames], columns = nomecolns, ggplot2::aes(color=tipo_zona, alpha = 0.25))      
gg
# ggsave(plot = gg, filename = paste0(wrkdir,PARAMETRO,"_P_wet_scatterplot.png"), dpi = 300,  
#        units = "in",width = 13, height = 10) 

# scatterplots pbl
nomecolns <- c("value" , "pbl00","pbl12", "pblmin","pblmax")
print(skimr::skim(subDati[,c("tipo_zona",nomecolns)] %>% group_by(tipo_zona) ))
gg <- ggpairs(subDati[,subMeteoNames], columns = nomecolns, ggplot2::aes(colour=tipo_zona, alpha = 0.25))      
gg
# ggsave(plot = gg, filename = paste0(wrkdir,PARAMETRO,"_pbl_scatterplot.png"), dpi = 300,  
#        units = "in",width = 13, height = 10) 

# scatterplots wind
nomecolns <- c("value" , "wdir", "wspeed","pwspeed")
print(skimr::skim(subDati[,c("tipo_zona",nomecolns)] %>% group_by(tipo_zona) ))
gg <- ggpairs(subDati[,subMeteoNames], columns = nomecolns, ggplot2::aes(colour=tipo_zona, alpha = 0.25))      
gg
# ggsave(plot = gg, filename = paste0(wrkdir,PARAMETRO,"_wind_scatterplot.png"), dpi = 300,  
#        units = "in",width = 13, height = 10) 

