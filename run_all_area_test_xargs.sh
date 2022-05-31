#!/bin/bash +x 

#run once first to generate weight file
#./process_area_weighted.py -f "/groups/ESS/ytang/nei2016v1/IC2016V1_EMIS_PREMERGED_2016fh_GRID/afdust_adj/emis_mole_afdust_adj_20160105_12us1_cmaq_cb6_2016fh_16j.ncf" -o "../../emissions/hemco/NEI2016v1/v2020-07/01/NEI2016v1_0.1x0.1_20160105_afdust_adj.nc" -t 'target.nc'




declare -a sectors=("afdust_adj" "ag" "airports" "nonpt" "nonroad" "np_oilgas" "onroad" "onroad_ca_adj" "onroad_can" "onroad_mex" "othafdust_adj" "othar" "othptdust_adj" "ptnonipm" "pt_oilgas" "rail" "rwc")

#declare -a sectors=("nonpt" "nonroad" "np_oilgas" "onroad" "onroad_ca_adj" "onroad_can" "onroad_mex" "othafdust_adj" "othar" "othptdust_adj" "ptnonipm" "pt_oilgas" "rail" "rwc")


#declare -a sectors=("onroad_mex" "othafdust_adj" "othar" "othptdust_adj" "ptnonipm" "pt_oilgas" "rail" "rwc")
declare -a months=("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12")

outdir=emissions/hemco/NEI2015v1/v2020-07

for sector in "${sectors[@]}"
do


for month in "${months[@]}"
do

mkdir -p ${outdir}/${month}/

	for day in "01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30" "31"; do
        echo ${day}
        done | xargs -I {} --max-procs 32 ./process_area_weighted.py -f "/groups/ESS/ytang/nei2016v1/IC2016V1_EMIS_PREMERGED_2016fh_GRID/${sector}/emis_mole_${sector}_2016${month}{}_12us1_cmaq_cb6_2016fh_16j.ncf" -o "${outdir}/${month}/NEI2016v1_0.1x0.1_2016${month}{}_${sector}.nc" -t 'target.nc' -w 'conservative_normed_299x459_288x640.nc' 


done

done
      
