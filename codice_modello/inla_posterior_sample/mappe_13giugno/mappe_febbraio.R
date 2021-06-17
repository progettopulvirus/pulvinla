rm(list=objects())
library("tidyverse")
library("scico")
library("raster")
library("regioniItalia")
library("patchwork")
library("sf")
library("ggnewscale")
library("bfastSpatial")
library("osmdata")
library("ggspatial")
library("stazioniMonitoraggio")




st_as_sf(stazioni,coords=c("st_x","st_y"),crs=4326)->sfStazioni
st_transform(sfStazioni,crs=32632)->sfStazioni
sfStazioni %>%
  filter(regione=="LOMBARDIA" & no2=="completa 2016-2020")->sfStazioni

MODELLO<-"METEO1920"
SOGLIA<-0
MAX_SOGLIA<-1
MESE<-2
TIPO<-"feriali"


st_union(lombardia,piemonte)->temp
st_union(temp,veneto)->temp
st_union(temp,emiliaromagna)->temp
st_union(temp,toscana)->temp
st_union(temp,trentino)->temp
st_union(temp,valleaosta)->temp
st_union(temp,friuliveneziagiulia)->temp
st_union(temp,liguria)->temp
temp->lombardia

st_transform(italia_senza_isole,crs=32632)->italia_senza_isole
st_intersection(lombardia,italia_senza_isole)->lombardia

filtra<-function(mylayer,soglia=0,nome=NULL,inverti=FALSE,max_soglia=0.8){
  
  if(!inverti){
    mylayer[mylayer<soglia & mylayer>0]<-NA
    mylayer[mylayer > -soglia & mylayer<0]<- NA
  }else{
    mylayer[mylayer> soglia]<-999
    mylayer[mylayer < -soglia]<- 999   
    mylayer[mylayer!=999]<-NA
  }
  
  mylayer[mylayer>  max_soglia]<-  max_soglia
  mylayer[mylayer< -max_soglia]<- -max_soglia
  
  names(mylayer)<-nome
  
  mylayer
  
} #fine filta


### graficoMappe, singolo modello (ovvero senza intersezione modello 2018 vs 2020 con modello 2019 vs 2020)

graficoMappe<-function(intero,mascherato,anni,tipo,max_soglia){
  
  numeroLayers<-nlayers(intero)
  
  purrr::map(1:numeroLayers,.f=function(ll){
    
    #mese mi serve per l'etichetta delle mappe, per distinguere mesi e settimane 
    names(intero[[ll]])->mese
    
    #  
    mascherato[[ll]]->rasterMascherato 
    
    tabularaster::as_tibble(intero[[ll]],xy=TRUE)->df
    
    #il raster delle sole zone statisticamente significative lo vogliamo trasformare in un poligono  
    #ovvero: vogliamo le mappe delle medie e dei poligoni che evidenzino le zone statisticamente significative  
    rasterMascherato[!is.na(rasterMascherato)]<-0
    rasterToPolygons(rasterMascherato,dissolve = TRUE)->poligoni
    
    st_as_sf(poligoni)->sfPoligoni
    st_crs(sfPoligoni)<-32632
    

    st_cast(sfPoligoni,"POLYGON")->tanti_poligoni
    tanti_poligoni[st_area(tanti_poligoni) > units::set_units(50000000,"m^2"),]->subPoligoni
    #nngeo::st_remove_holes(subPoligoni)->subPoligoni
    
    ifelse(grepl("mese1",mese),str_replace(mese,"mese1","January"),str_replace(mese,"mese2","February"))->etichetta_asse_y
    ifelse(grepl("mese3",mese),str_replace(mese,"mese3","March"),str_replace(mese,"mese4","April"))->etichetta_asse_y
    
    ggplot()+
      geom_sf(data=lombardia,fill="transparent",lwd=0.5)+
      geom_tile(data=df,aes(x=x,y=y,fill=cellvalue),alpha=1)+
      scico::scale_fill_scico(palette="broc",limits=c(-max_soglia,max_soglia),na.value=NA,direction=1,guide=guide_colorbar(barheight=unit(10,"cm"),title=str_wrap("% relative change",10)))+
      geom_sf(data=subPoligoni,fill="transparent",color="#333333")+
      #geom_sf(data=sfStazioni,pch=21,fill="firebrick",size=0.1)+
      ylab(etichetta_asse_y)+
      theme_minimal()+
      theme(panel.grid = element_blank(),axis.title.x = element_blank(),axis.text = element_blank(),text = element_text(size=24))->grafico

    
    png(glue::glue("mappe_{MESE}.{ll}_{TIPO}.png"),width=1024,height=1024)
    print(grafico)
    dev.off()    

    grafico
    
  })->listaGrafici
  

  
  purrr::reduce(listaGrafici[1:length(listaGrafici)],.f=`+`)+
    plot_layout(guides="collect",ncol = 2)+plot_annotation(title="",subtitle="")


  
} #fine grafico


preparaRasters<-function(.ffile,max_soglia){ 

  purrr::map(paste0("mese",MESE),.f=function(nomeMese){
    
    .ffile[grepl(nomeMese,.ffile)]->fileMese
    
    #media mascherata: ovvero le zone delle mappe statisticamente significative in base ai percentili calcolati con CDO  
    fileMese[grep("Masked",fileMese)]->fileMasked
    #stessa mappa di sopra ma non mascherata, ovvero contiene incrementi e decrementi significativi o meno  
    fileMese[!grepl("Masked",fileMese)]->fileIntero  
    
    brick(fileMasked)->brickMasked
    brick(fileIntero)->brickIntero
    
    nlayers(brickMasked)->numeroLayers
    stopifnot(nlayers(brickMasked)==nlayers(brickIntero))
    
    purrr::map(1:numeroLayers,.f=function(ll){
      
      brickMasked[[ll]]->mascherato1
      brickIntero[[ll]]->intero1
      
      filtra(mascherato1,soglia=SOGLIA,nome=nomeMese,max_soglia = max_soglia)->mascherato1
      filtra(    intero1,soglia=SOGLIA,nome=nomeMese,max_soglia = max_soglia)->intero1
      
      list("masked"=mascherato1,"not_masked"=intero1)
      
        
    })->lista
    
    purrr::map(lista,"masked") %>% brick->mbrick
    purrr::map(lista,"not_masked") %>% brick->ibrick
    
    list("masked"=mbrick,"not_masked"=ibrick)
    
  })

}#fine funzione preparaRasters   

### medie non mascherate

preparaDati<-function(.x){ 
  
  names(.x)->nomiLayers

  tabularaster::as_tibble(.x) %>%
    mutate(periodo=nomiLayers[dimindex]) %>%
    mutate(dimindex=factor(dimindex,ordered=TRUE,levels=unique(sort(as.integer(dimindex))))) %>%
    mutate(periodo=factor(periodo,ordered=TRUE,levels=nomiLayers))
  
  
}#preparaDati 


istogramma<-function(.x){ 
  
  cut(.x$cellvalue,breaks=seq(-1,0,by=0.2),include.lowest=FALSE,right=TRUE)->.x$decremento
  
  ggplot(data=.x %>% filter(!is.na(cellvalue)))+
    geom_histogram(aes(x=cellvalue,fill=decremento),color="#333333")+
    facet_wrap(~periodo)+
    scale_fill_scico_d()+
    scale_y_continuous(limits=c(0,15000))+
    ylab("Decrementi")+
    xlab("")+
    theme_bw()

}#istogramma 


scatola<-function(.x){ 

  ggplot(data=.x %>% mutate(mese=str_extract(periodo,"^[[:alpha:]]+[0-9]")) %>% filter(cellvalue<0))+
    geom_boxplot(aes(x=periodo,y=cellvalue,fill=mese))+
    scale_fill_scico_d()+
    scale_y_continuous(limits=c(-0.8,0))+
    ylab("Decrementi")+
    xlab("")+
    theme_bw()+
    theme(panel.grid.major.y = element_line(colour="#333333"),axis.text.x.bottom = element_text(angle=90))
  
}#scatola 


### Giorni feriali




#mappe con media mascherata e non
list.files(pattern = glue::glue("^media.*_{TIPO}_modello{MODELLO}.+nc"),full.names = TRUE,include.dirs = TRUE,recursive = TRUE)->ffile

#leggi i rasters e preparali
preparaRasters(.ffile=ffile,max_soglia = MAX_SOGLIA)->listaBrick

#maskedBrick: brick che contiene le mappe di tutti i mesi mascherati (solo medie statisticamente significative)
#not_maskedBrick: mappe non mascherate
purrr::map(listaBrick,"masked") %>% brick->maskedBrick
purrr::map(listaBrick,"not_masked") %>% brick->not_maskedBrick

graficoMappe(intero=not_maskedBrick,mascherato=maskedBrick,anni=MODELLO,tipo=TIPO,max_soglia = MAX_SOGLIA)->mappe_no2

png(glue::glue("mappe_{MESE}_{TIPO}.png"),width=1024,height=1024)
print(mappe_no2)
dev.off()




### Boxplot e istogrammi 


purrr::map(list(not_maskedBrick,maskedBrick),.f=preparaDati)->dfFeriali


purrr::map(dfFeriali,.f=istogramma)->istogramma_feriali
print(purrr::reduce(istogramma_feriali,.f=`+`)+plot_layout(guides="collect")+plot_annotation(title="Giorni feriali, decrementi rispetto al 2020") & theme(legend.position = "bottom"))

purrr::map(dfFeriali,.f=scatola)->scatola_feriali
print(purrr::reduce(scatola_feriali,.f=`+`)+plot_layout(guides="collect")+plot_annotation(title="Giorni feriali, decrementi rispetto al 2020") & theme(legend.position = "bottom"))





