#!/bin/bash -x

#load monet modules
#module purge
#module use -a /gpfs/hps3/emc/naqfc/noscrub/Barry.Baker/modulefiles
#module load anaconda3/latest

#get the current date
yyyymmdd=$(date +%Y%m%d)
yyyymmdd='20190203'

#go to data directory
data_dir=/gpfs/hps2/ptmp/Barry.Baker/fengsha/gfs.$yyyymmdd/00
cd data_dir

#make spatial plots for all regions
for i in 'NAU' 'SAU' 'AMZ' 'SSA' 'CAM' 'WNA' 'CNA' 'ENA' 'ALA' 'GRL' 'MED' 'NEU' 'WAF' 'EAF' 'SAF' 'SAH' 'SEA' 'EAS' 'SAS' 'CAS' 'TIB' 'NAS'; do
  echo ${i}
done | xargs -I {} ./fv3_spatial.py -f "gfs.t00z.master.grb2f???.entire_atm.nc" -v pm25aod550 salt25aod550 dust25aod550 sulf25aod550 bc25aod550 oc25aod550

#
