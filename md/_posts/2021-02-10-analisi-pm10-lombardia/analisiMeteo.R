rm(list=objects())
library("tidyverse")
library("lombardia")
library("datiMeteo")
library("stazioniMonitoraggio")
library("regioniItalia")
library("sf")
library("scico")
library("ggspatial")
library("data.tree")
library("raster")

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

#creo un oggetto spaziale
stazioni[,c("station_eu_code","nome_stazione","st_x","st_y","altitudine")]->puntiStazione
st_as_sf(puntiStazione ,coords = c("st_x","st_y"),crs=4326)->puntiStazione
st_transform(puntiStazione,crs = 32632)->puntiStazioneUTM

#stazioni che appartengono a "Area" (regione, provincia..)
lombardia->Area
st_intersection(puntiStazioneUTM,Area)->stazioniArea

#stazioni che appartengono ad Area e con dati in "inquinante"
stazioniArea[stazioniArea$station_eu_code %in% codiciStazioni,]->stazioniArea


#' quante stazioni
length(unique(stazioniArea$station_eu_code))

raster("dem.tif")->dem
raster::crop(dem,lombardia)->demArea

theme_set(theme_bw()+theme(panel.grid.minor = element_blank(),
                           panel.grid.major.x = element_blank()))

ggplot()+
  ggspatial::layer_spatial(data=demArea)+
  ggspatial::geom_sf(data=Area,fill="transparent",colour="black")+
  geom_sf(data=stazioniArea,aes(fill=altitudine),pch=21)+
  scale_fill_scico(palette="bamako")->mappa


print(mappa)

#' numero stazioni
nrow(stazioniArea)


#selezioniamo dati con codici in "stazioniArea" (stazioni che cadono in Area e con dati)
#selezioniamo serie nel periodo gennaio-giugno
dati %>%
  filter(station_eu_code %in% stazioniArea$station_eu_code) %>%
  filter(mm %in% 1:6)->subDati

left_join(subDati,stazioniMonitoraggio::stazioni[,c("station_eu_code","nome_stazione","tipo_zona")]) %>%
  mutate(tipo_zona=as.factor(tipo_zona))->subDati


print(skimr::skim(subDati[,c("yy","mm","dd","value","tipo_zona")] ))

print(skimr::skim(subDati %>% dplyr::select(value,tipo_zona) %>% group_by(tipo_zona)))


ggplot(data=subDati)+
  geom_bar(aes(x=yy))+
  facet_wrap(~station_eu_code)

ggplot(data=subDati)+
  geom_boxplot(aes(x=station_eu_code,y=value))+
  geom_hline(yintercept = 50,colour="red")+
  facet_wrap(~tipo_zona,scales = "free_x",nrow = 3)



