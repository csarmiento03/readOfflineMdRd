#!/bin/bash

#source /cr/aera02/user/sarmiento/new_apps/apps_Offline/read_ADST/ADST_readout.cc

tag='trunk_r31074'
sim='joint'
#particle='proton'
particle='separated'
#hadronicModel='qgsjet'
#hadronicModel='sibyll'
#hadronicModel='epos-lhc'

source  /cr/aera02/user/sarmiento/offlineTrunkRd/setOfflineTrunk.sh

particle=$1
refdist=$2

#/cr/aera02/user/sarmiento/AppMdRdParameters/AppMdRdReconstruction/ADST_refdist_140_Iron/ADST_event_Iron_DAT104781.root

for file in /cr/aera02/user/sarmiento/AppMdRdParameters/AppMdRdReconstruction/ADST_refdist_${refdist}_${particle}/ADST*; do
  echo $file
  ./ADST_readout_murd $file ${particle}RefDist${refdist} &>> log${particle}${refdist}_murd.txt

#break
done

