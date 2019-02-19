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
mpl.use('agg')


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
    from monet.models import cmaq
    f = cmaq.open_mfdataset(finput)
    return f


def load_paired_data(fname):
    return pd.read_hdf(fname)

def make_spatial_plot(da, df, proj): 
    cbar_kwargs = dict(aspect=30,shrink=.8)#dict(aspect=30)                       
    extent = [da.longitude.min().values+1,da.longitude.max().values-10,da.latitude.min().values+3,da.latitude.max().values-4] 
    ax = da.monet.quick_map(cbar_kwargs=cbar_kwargs, figsize=(15, 8), map_kwarg={'states': True, 'crs': proj,'extent':extent},robust=True) 
    plt.gcf().canvas.draw() 
    plt.tight_layout(pad=-1) 
    date = pd.Timestamp(da.time.values) 
    cbar = ax.figure.get_axes()[1] 
    vmin, vmax = cbar.get_ybound() 
    vars = df.keys() 
    varname = [x for x in vars if x not in ['latitude','longitude']][0] 
    ax.scatter(df.longitude.values,df.latitude.values,s=25,c=df[varname],transform=ccrs.PlateCarree(),edgecolor='w',linewidth=.08,vmin=vmin,vmax=vmax) 
    ax.set_extent(extent,crs=ccrs.PlateCarree())
    plt.savefig(date.strftime('naqfc_' + da.name + '_sp.%Y%m%d%H.jpg'), dpi=100)
    plt.close()

# def make_spatial_plot(da, df, proj):
#     cbar_kwargs = dict(aspect=30,shrink=.8)#dict(aspect=30)
#     ax = da.monet.quick_map(cbar_kwargs=cbar_kwargs, figsize=(
#         12, 8), map_kwarg={'states': True, 'crs': proj})
#     date = pd.Timestamp(da.time.values)
#     cbar = ax.figure.get_axes()[1]
#     vmin, vmax = cbar.get_ybound()
#     vars = df.keys()
#     varname = [x for x in vars if x not in ['latitude','longitude']][0]
#     print(varname)
#     ax.scatter(df.longitude, df.latitude, c=df[varname])
#     plt.tight_layout()
#     plt.savefig(date.strftime(
#         'naqfc_' + da.name + '_sp.%Y%m%d%H.jpg'), dpi=100)
#     plt.close()

def make_plots(finput, paired_data, variable, obs_variable, verbose, region, epa_region):

    # open the files
    f = open_cmaq(finput)
    # get map projection
    proj = map_projection(f)
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
            make_spatial_plot(obj.sel(time=t), odf, proj)
#            plots.append(dask.delayed(make_spatial_plot)
#                         (obj.sel(time=t), odf, proj))
#        plots dask.delayed(make_spatial_plot)(
#            obj.sel(time=t), proj) for t in obj.time]
#    dask.delayed(plots).compute()
    # if paired_data is not None:
    #     ov = obs_variable[[index, 'latitude', 'longitude'].loc[obs_variable.time == t]
    # ax.scatter()

    print(variable)


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
    if args.region is None:
        region = None
    else:
        region = args.region
    if args.epa_region is None:
        epa_region = None
    else:
        epa_region = args.epa_region

    make_plots(finput, paired_data, variable,
               obs_variable, verbose, region, epa_region)
