#!/usr/bin/env python

__author__  = 'Patrick Campbell'
__email__   = 'patrick.c.campbell@noaa.gov'
__license__ = 'GPL'



#Simple MONET utility to command line pair model vs. observations

import os
from glob import glob
import sys
#sys.path.append('/home/patrickc/MONET/')
#os.chdir('/home/patrickc/MONET/scripts/')

import subprocess
from distutils.spawn import find_executable
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import monet  
from monet.util.tools import long_to_wide
import pandas as pd

def  pair_point(da,df,sub_map,interp):
     dfpair=da.monet.combine_point(df,sub_map,method=interp,reuse_weights=True)
     print(dfpair)
     return dfpair

def  get_aqs(start,end,datapath=None,species=None,verbose=False):
     dates = pd.date_range(start=start, end=end, freq='H')
     monet.obs.aqs.datadir=datapath
     dfaqs = monet.obs.aqs.add_data(dates,param=species)
     dfwide   = long_to_wide(dfaqs)
     return dfwide

def  get_airnow(start,end,datapath=None,verbose=False):
     dates = pd.date_range(start=start, end=end, freq='H')
     monet.obs.airnow.datadir=datapath
     dfairnow = monet.obs.airnow.add_data(dates)
     dfwide   = long_to_wide(dfairnow)
     #make sure there are no duplicates
     return dfwide.drop_duplicates(subset=['time','siteid'])

def  open_cmaq(finput,verbose=False):
     dset=monet.models.cmaq.open_mfdataset(finput)
     return dset


if __name__ == '__main__':

    parser = ArgumentParser(description='pairs cmaq model data to aqs observations', formatter_class=ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-f', '--files',       help='string input model file directory/names', type=str, required=True)
#    parser.add_argument('-s', '--startdates',  help='string input start date for pairing YYYY-MM-DD HH:MM:SS', type=str, required=True)
#    parser.add_argument('-e', '--enddates',    help='string input end date for pairing YYYY-MM-DD HH:MM:SS', type=str, required=True)
    parser.add_argument('-x', '--species',     help='string input for obs species-variables to pair',type=str,nargs='+', required=False, default=['OZONE','PM2.5'])
    parser.add_argument('-o', '--output',      help='string output path for paired dataframe, stats, plots', type=str, required=False,default='./')
    parser.add_argument('-p', '--path',        help='string path to director of network observations', type=str, required=False, default='/data/aqf2/barryb/5xpm/AQS_DATA/')
    parser.add_argument('-n', '--network',     help='string input data network name: airnow, aqs', type=str, required=False, default='airnow')
    parser.add_argument('-m', '--model',       help='input model: cmaq, fv3, hysplit (not-ready), or camx (not-ready)', type=str, required=False, default='cmaq')
    parser.add_argument('-i', '--interp',      help='xesmf interpolation scheme, bilinear, conservative, nearest_s2d, nearest_d2s, patch', type=str, required=False, default='bilinear')
    parser.add_argument('-v', '--verbose',     help='print debugging information', action='store_true', required=False)
    args = parser.parse_args()

    finput  = args.files
#    start   = args.startdates
#    end     = args.enddates 
    species = args.species
    output  = args.output
    datapath= args.path
    network = args.network
    model   = args.model
    interp  = args.interp
    verbose = args.verbose
    #reads model output (cmaq default)

    if model == 'cmaq':
    	 da=open_cmaq(finput,verbose=verbose)  
    else:
         print('Must enter cmaq model right now')
         raise RuntimeError

#retrieves data observations and formats pandas dataframe (airnow default)
    if network == 'airnow':
         start = da.time.to_index()[0]
         end = da.time.to_index()[1]
         df=get_airnow(start,end,datapath)
    elif network == 'aqs':
         df=get_aqs(start,end,datapath,species)
    else:
         print('Must enter airnow or aqs right now')
         raise RuntimeError
    
    #pairs surface point-type observations with 2D model parameters
    
    if network == 'airnow' and model == 'cmaq': 
         mapping_table = {'OZONE':'O3', 'PM2.5':'PM25_TOT', 'PM10':'PMC_TOT'}
         sub_map = {i: mapping_table[i] for i in species if i in mapping_table}
         dfpair=pair_point(da,df,sub_map,interp)
         print(dfpair)
         dfpair.to_csv('out_dfpair.csv')
         
         
    else:
         print('Must pair airnow right now')
         raise RuntimeError
    
    
    sys.exit(0)
    


