#!/usr/bin/env python

###############################################################
# < next few lines under version control, D O  N O T  E D I T >
# $Date: 2018-03-29 10:12:00 -0400 (Thu, 29 Mar 2018) $
# $Revision: 100014 $
# $Author: Barry.Baker@noaa.gov $
# $Id: nemsio2nc4.py 100014 2018-03-29 14:12:00Z Barry.Baker@noaa.gov $
###############################################################

__author__ = 'Barry Baker'
__email__ = 'Barry.Baker@noaa.gov'
__license__ = 'GPL'

import os
import subprocess
import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

import dask
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import cartopy.crs as ccrs
plt.ioff()


'''
Simple utility to make spatial plots from the NAQFC forecast and overlay observations
'''


def map_projection(f):
    import cartopy.crs as ccrs
    proj = ccrs.LambertConformal(
        central_longitude=f.XCENT, central_latitude=f.YCENT)
    return proj


def chdir(fname):
    dir_path = os.path.dirname(os.path.realpath(fname))
    os.chdir(dir_path)
    return os.path.basename(fname)


def open_cmaq(finput):
    from monet.models import fv3chem
    f = fv3chem.open_mfdataset(finput)
    return f


def load_paired_data(fname):
    return pd.read_hdf(fname)

def make_spatial_plot(da, df, out_name):
    cbar_kwargs = dict(aspect=30,shrink=.8,orientation='horizontal')#dict(aspect=30)
    levels = [0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.5,2,2.5]
    ax = da.where(da > .05).monet.quick_map(cbar_kwargs=cbar_kwargs, figsize=(11.31,  7.67),levels=levels,cmap='plasma')#robust=True)
    date = pd.Timestamp(da.time.values)
    cbar = ax.figure.get_axes()[1]
    vmin, vmax = cbar.get_ybound()
    vars = df.keys()
    varname = [x for x in vars if x not in ['latitude','longitude']][0]
    odf = df.dropna(subset=['latitude','longitude',varname])
    odf.plot.scatter(x='longitude',y='latitude',c=varname,colorbar=False,vmin=0.05,vmax=2.5,s=30,edgecolor='w',linewidth=.08,cmap='plasma',ax=ax)
    plt.tight_layout(pad=0)
    plt.savefig(date.strftime(out_name + da.name + '_sp.%Y%m%d%H.jpg'), dpi=100)
    plt.close()

def make_spatial_bias_plot(df,out_name,**kwargs):
    from monet import plots
    ax = plots.sp_scatter_bias(df,**kwargs)
    date = df.time.min()
    plt.title(date.strftime('time=%Y/%m/%d %H:00 | FV3 - AERONET (AOD)'))
    plt.tight_layout(pad=0)
    plots.savefig(date.strftime(out_name + 'sb.%Y%m%d%H.jpg'),dpi=100),decorate=True)
    plt.close()

def make_plots(finput, paired_data, variable, obs_variable, verbose, region, epa_region, out_name):

    # open the files
    f = open_cmaq(finput)
    # get map projection
#    proj = map_projection(f)
    if paired_data is not None:
        df = load_paired_data(paired_data)
    # loop over varaible list
    plots = []
    for index, var in enumerate(variable):
        obj = f[var]
        # loop over time
        for t in obj.time:
            date = pd.Timestamp(t.values)
            print(date)
            odf = df.loc[df.time == pd.Timestamp(t.values),['latitude','longitude',obs_variable[index]]]
            make_spatial_plot(obj.sel(time=t), odf, out_name)
            make_spatial_bias_plot(df,out_name

def make_boxplot_giorgi(paired_data,savename):
    from monet.util.tools import get_giorgi_region_df as ggrd
    from monet.plots import savefig
    import seaborn as sns
    df = paird_data.copy()
    df = ggrd(df)
    dfa = df.dropna(subset=['aod_550nm','pm25aod550'])['aod_550nm','GIORGI_ACRO']
    dfm = df.dropna(subset=['aod_550nm','pm25aod550'])['pm25aod550','GIORGI_ACRO']
    dfa['Legend'] = 'AERONET'
    dfm['Legend'] = 'FV3CHEM'
    dfa.rename({'aod_550nm':'AOD'},axis=1,inplace=True)
    dfm.rename({'pm25aod550':'AOD'},axis=1,inplace=True)
    dfn = pd.concat([dfa,dfm],ignore_index=True)
    f,ax = plt.subplots(figsize=(12,7))
    sns.boxplot(ax=ax,x='GIORGI_ACRO',y='AOD',hue='Legend',data=dfn)
    sns.despine()
    plt.tight_layout(pad=0)
    savefig(savename,dpi=100)
    plt.close()
    
if __name__ == '__main__':

    parser = ArgumentParser(description='Make Spatial Plots for each time step in files',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-f', '--files', help='input file names', type=str, required=True)
    parser.add_argument(
        '-p', '--paired_data', help='paired data input file names', type=str, required=False)
    parser.add_argument(
        '-v', '--variable', nargs='+', help='variable name to plot', required=True)
    parser.add_argument(
        '-v2', '--obs_variable', nargs='+', help='input file names', required=False)
    parser.add_argument('-vb', '--verbose', help='print debugging information',
                        action='store_true', required=False)
    parser.add_argument(
        '-r', '--region', help='NCEP Region name', nargs='+', required=False)
    parser.add_argument('-e', '--epa_region', help='EPA Air Shed name',
                        nargs='+', required=False)
    parser.add_argument('-d', '--directory', help='Output Directory',
                        type=str, required=False)
    parser.add_argument('-n', '--output_name', help='Output base name',
                        type=str, required=False)
    parser.add_argument('-s', '--suppress_xwindow', help='Suppress X Window',
                        action='store_true', required=False)
    args = parser.parse_args()

    finput = args.files
    paired_data = args.paired_data
    variable = args.variable
    obs_variable = args.obs_variable
    verbose = args.verbose
    out_name = args.output_name
    if args.region is None:
        region = None
    else:
        region = args.region
    if args.epa_region is None:
        epa_region = None
    else:
        epa_region = args.epa_region

    make_plots(finput, paired_data, variable,
               obs_variable, verbose, region, epa_region,out_name)
