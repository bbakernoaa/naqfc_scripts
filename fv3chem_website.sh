#!/bin/bash -x

# load monet modules
module purge
module use -a /gpfs/hps3/emc/naqfc/noscrub/Barry.Baker/modulefiles
module load anaconda3/latest

export OMP_NUM_THREADS=2
# get the current date
yyyymmdd=$(date +%Y%m%d)
yyyymmdd='20190203'

# go to data directory
data_dir=/gpfs/hps2/ptmp/Barry.Baker/fengsha/gfs.$yyyymmdd/00
cd data_dir

# grib files
files="gfs.t00z.master.grb2f???.entire_atm.nc"

# convert files
ls -t ${files} | xargs -I {} --max-procs 15 fv3grib2nc4.py -f {}

#spatial plotting
spatial=/gpfs/hps3/emc/naqfc/noscrub/Barry.Baker/naqfc_scripts/fv3_spatial.py
pair=/gpfs/hps3/emc/naqfc/noscrub/Barry.Baker/naqfc_scripts/fv3_pair_aeronet.py
box=/gpfs/hps3/emc/naqfc/noscrub/Barry.Baker/naqfc_scripts/fv3_aeronet_box.py
bias=/gpfs/hps3/emc/naqfc/noscrub/Barry.Baker/naqfc_scripts/fv3_aeronet_spatial_bias.py
taylor=/gpfs/hps3/emc/naqfc/noscrub/Barry.Baker/naqfc_scripts/fv3_aeronet_taylor.py

# make spatial plots for all regions
for i in 'global' 'NAU' 'SAU' 'AMZ' 'SSA' 'CAM' 'WNA' 'CNA' 'ENA' 'ALA' 'GRL' 'MED' 'NEU' 'WAF' 'EAF' 'SAF' 'SAH' 'SEA' 'EAS' 'SAS' 'CAS' 'TIB' 'NAS'; do
  echo ${i}
done | xargs -I {} --max-procs 10 ${spatial} -f ${files} -v pm25aod550 salt25aod550 dust25aod550 sulf25aod550 bc25aod550 oc25aod550

# pair the data
${pair} -f ${files} -o fv3chem_aeronet.hdf

# create spatial bias plots
for i in 'global' 'NAU' 'SAU' 'AMZ' 'SSA' 'CAM' 'WNA' 'CNA' 'ENA' 'ALA' 'GRL' 'MED' 'NEU' 'WAF' 'EAF' 'SAF' 'SAH' 'SEA' 'EAS' 'SAS' 'CAS' 'TIB' 'NAS'; do
  echo ${i}
done | xargs -I {} --max-procs 5 ${bias} -p fv3chem_aeronet.hdf -v pm25aod550 -vp aod_550nm

# create box plot
${box} -p fv3chem_aeronet.hdf

# create fv3_aeronet_taylor
for i in 'global' 'NAU' 'SAU' 'AMZ' 'SSA' 'CAM' 'WNA' 'CNA' 'ENA' 'ALA' 'GRL' 'MED' 'NEU' 'WAF' 'EAF' 'SAF' 'SAH' 'SEA' 'EAS' 'SAS' 'CAS' 'TIB' 'NAS'; do
  echo ${i}
done | xargs -I {} --max-procs 5 ${taylor} -p fv3chem_aeronet.hdf -v pm25aod550 -vp aod_550nm

# make GIFS
##########################################################################################################
for i in 'global' 'NAU' 'SAU' 'AMZ' 'SSA' 'CAM' 'WNA' 'CNA' 'ENA' 'ALA' 'GRL' 'MED' 'NEU' 'WAF' 'EAF' 'SAF' 'SAH' 'SEA' 'EAS' 'SAS' 'CAS' 'TIB' 'NAS'; do
  for j in 'pm25aod550' 'salt25aod550' 'dust25aod550' 'sulf25aod550' 'bc25aod550' 'oc25aod550'; do
    echo "${i}.${j}"
  done
done | xargs -I {} --max-procs 10 convert -delay 40 -loop 0 {} FV3CHEM.{}.sp*.jpg FV3CHEM.{}.sp.gif

for i in 'global' 'NAU' 'SAU' 'AMZ' 'SSA' 'CAM' 'WNA' 'CNA' 'ENA' 'ALA' 'GRL' 'MED' 'NEU' 'WAF' 'EAF' 'SAF' 'SAH' 'SEA' 'EAS' 'SAS' 'CAS' 'TIB' 'NAS'; do
  for j in 'pm25aod550' 'salt25aod550' 'dust25aod550' 'sulf25aod550' 'bc25aod550' 'oc25aod550'; do
    echo "${i}.${j}"
  done
done | xargs -I {} --max-procs 10 convert -delay 40 -loop 0 {} FV3CHEM.{}.sb*.jpg FV3CHEM.{}.sb.gif
##########################################################################################################

# Transfer data back to aaqest
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ls -t *.jpg | xargs -I {} --max-procs 4 scp {} barryb@aaqest.arl.noaa.gov:/data/aqf2/testbed-bak/www/testbed/TCHAI/hkim/${yyyymmdd}/
ls -t *.gif | xargs -I {} --max-procs 4 scp {} barryb@aaqest.arl.noaa.gov:/data/aqf2/testbed-bak/www/testbed/TCHAI/hkim/${yyyymmdd}/
