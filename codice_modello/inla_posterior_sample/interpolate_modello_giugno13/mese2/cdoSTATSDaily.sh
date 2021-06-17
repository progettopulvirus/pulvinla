#!/bin/bash

mese=2
modello="METEO1920"
regione="lombardia"
mkdir ../ris${mese}

TIPO=("feriali" "festivi")

for tipo in ${TIPO[@]}; do

fileDaElaborare=$(ls ${tipo}_*.nc)

echo "--------->>>>  ${tipo}"

for ffile in ${fileDaElaborare[@]}; do

	cdo timmean ${ffile} media_${ffile}.nc
	cdo timmin  ${ffile} minimo_${ffile}.nc
	cdo timmax  ${ffile} massimo_${ffile}.nc	
        cdo timpctl,2.5 ${ffile} minimo_${ffile}.nc massimo_${ffile}.nc pctl2.5_${ffile}.nc
        cdo timpctl,97.5 ${ffile} minimo_${ffile}.nc massimo_${ffile}.nc pctl97.5_${ffile}.nc

	rm -rf minimo_${ffile}.nc massimo_${ffile}.nc

	#cdo merge pctl97.5_${ffile}.nc media_${ffile}.nc pctl2.5_${ffile}.nc
done



cdo merge $(ls pctl97.5_${tipo}_*.nc) ../ris${mese}/pctl97.5_${tipo}_modello${modello}_mese${mese}_regione${regione}.nc
cdo merge $(ls pctl2.5_${tipo}_*.nc) ../ris${mese}/pctl2.5_${tipo}_modello${modello}_mese${mese}_regione${regione}.nc
cdo merge $(ls media_${tipo}_*.nc) ../ris${mese}/media_${tipo}_modello${modello}_mese${mese}_regione${regione}.nc

rm -rf media*.nc
rm -rf *97*.nc
rm -rf *2.5*.nc

cdo mul ../ris${mese}/pctl97.5_${tipo}_modello${modello}_mese${mese}_regione${regione}.nc ../ris${mese}/pctl2.5_${tipo}_modello${modello}_mese${mese}_regione${regione}.nc ../ris${mese}/prodotto.nc

cdo gec,0 ../ris${mese}/prodotto.nc ../ris${mese}/maschera.nc

cdo ifthen ../ris${mese}/maschera.nc ../ris${mese}/media_${tipo}_modello${modello}_mese${mese}_regione${regione}.nc ../ris${mese}/mediaMasked_${tipo}_modello${modello}_mese${mese}_regione${regione}.nc

done


