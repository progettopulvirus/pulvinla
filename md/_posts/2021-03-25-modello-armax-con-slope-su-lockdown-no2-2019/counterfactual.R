rm(list=objects())
library("tidyverse")

readRDS("_dati_staz_no2.RDS")->dati
dati%>%
  dplyr::select(station_eu_code,date,value,predfitted,predcounterf) %>%
  gather(key="tipo_serie",value="value",-station_eu_code,-date) %>%
  mutate(tipo_serie=case_when(tipo_serie=="value"~"Serie osservata",
                              tipo_serie=="predfitted"~"Forecast",
                              TRUE~"Counterfactual"))->dati2

ggplot(data=dati2,aes(x=date,y=value))+
  geom_line(aes(color=tipo_serie))+
  facet_wrap(~station_eu_code,ncol=8,scales = "free_y")+
  theme_bw()
