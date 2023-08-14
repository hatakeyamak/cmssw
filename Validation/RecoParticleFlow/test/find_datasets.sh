#!/bin/bash

export release=CMSSW_13_3_0_pre
export year=2023

for value in RelValQCD_FlatPt_15_3000HS_14 RelValZEE_14 RelValZMM_14 RelValTenTau_15_500 RelValNuGun
do
    echo $value
    echo $release
    dasgoclient --query="dataset dataset=/"${value}"/*"${release}"*"${year}"*/GEN-SIM-DIGI-RAW"
done



