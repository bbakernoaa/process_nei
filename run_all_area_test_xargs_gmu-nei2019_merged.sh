#!/bin/bash +x 

#Regrid Settings
sourceres=GMU_NEMO_1km
sourcegrid=3177x5397
#0.1 degree
targetres=0.1_degree
targetgrid=317x735
#0.03 degree
#targetres=0.03_degree
#targetgrid=1054x2448
#0.01 degree
#targetres=0.01_degree
#targetgrid=3162x7342

outdir=GMU-NEI2019v1/v2023-03/
mkdir -p $outdir

#May need run once first to generate target and weight files if they don't exist
#./process_area_weighted.py -f "/groups/ESS3/sma8/nei2019/2019ge_19j/premerged/afdust/emis_mole_afdust_20190101_US01_cmaq_cb6ae7_2019ge_19j.ncf" -o "${outdir}/test_target.nc" -t "target_${targetres}.nc"

declare -a sectors=("all")

declare -a months=("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12")

for sector in "${sectors[@]}"
do

indir=/scratch/sma8/data/US01/merged
emiprefix=emis_mole_${sector}
emisuffix=US01_cmaq_cb6ae7_2019ge_19j.ncf

for month in "${months[@]}"
do

mkdir -p ${outdir}/${month}/

	for day in "01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30" "31"; do
        echo ${day}
        done | xargs -I {} --max-procs 10 ./process_area_weighted.py -f "${indir}/${emiprefix}_2019${month}{}_${emisuffix}" -o "${outdir}/${month}/NEI2019v1_${targetres}_2019${month}{}_${sector}.nc" -t "target_${targetres}.nc" -s "${sourceres}_area.nc"  -w "conservative_normed_${sourcegrid}_${targetgrid}.nc" 
#        ./process_area_weighted.py -f "${indir}/${emiprefix}_2019${month}${day}_${emisuffix}" -o "${outdir}/${month}/NEI2019v1_${targetres}_2019${month}${day}_${sector}.nc" -t "target_${targetres}.nc" -s "${sourceres}_area.nc"  -w "conservative_normed_${sourcegrid}_${targetgrid}.nc"
#        done

done

done
      
