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

import cartopy.crs as ccrs
import dask
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

import monet

sns.set_context('notebook')

plt.ioff()
'''
Simple utility to make spatial plots from the NAQFC forecast and overlay observations
'''


def chdir(fname):
    dir_path = os.path.dirname(os.path.realpath(fname))
    os.chdir(dir_path)
    return os.path.basename(fname)


def open_fv3chem(finput):
    from monet.models import fv3chem
    f = fv3chem.open_mfdataset(finput)
    return f


def load_paired_data(fname):
    return pd.read_hdf(fname)


def make_spatial_plot(da, df, out_name):
    cbar_kwargs = dict(
        aspect=30, shrink=.8, orientation='horizontal')  # dict(aspect=30)
    levels = [
        0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2,
        2.5
    ]
    ax = da.where(da > .05).monet.quick_map(
        cbar_kwargs=cbar_kwargs,
        figsize=(11, 6.5),
        levels=levels,
        cmap='cividis')  # robust=True)
    date = pd.Timestamp(da.time.values)
    if df is not None:
        cbar = ax.figure.get_axes()[1]
        vmin, vmax = cbar.get_ybound()
        vars = df.keys()
        varname = [x for x in vars if x not in ['latitude', 'longitude']][0]
        odf = df.dropna(subset=['latitude', 'longitude', varname])
        odf.plot.scatter(
            x='longitude',
            y='latitude',
            c=varname,
            colorbar=False,
            vmin=0.05,
            vmax=2.5,
            s=30,
            edgecolor='w',
            linewidth=.08,
            cmap='plasma',
            ax=ax)
    plt.tight_layout(pad=0)
    savename = "{}.{}".format(out_name, date.strftime('sp.%Y%m%d%H.jpg'))
    print(savename)
    monet.plots.savefig(savename, bbox_inches='tight', dpi=100, decorate=True)
    plt.close()


# def make_spatial_bias_plot(df, out_name, **kwargs):
#     ax = monet.plots.sp_scatter_bias(df, **kwargs)
#     date = df.time.min()
#     plt.title(date.strftime('time=%Y/%m/%d %H:00 | FV3 - AERONET (AOD)'))
#     plt.tight_layout(pad=0)
#     savename = "{}.{}".format(out_name, date.strftime('sb.%Y%m%d%H.jpg'))
#     print(savename)
#     monet.plots.savefig(savename, bbox_inches='tight', dpi=100, decorate=True)
#     plt.close()


def make_plots(f, df, variable, obs_variable, out_name):
    # loop over varaible list
    plots = []
    for index, var in enumerate(variable):
        obj = f[var]
        # loop over time
        for t in obj.time:
            date = pd.Timestamp(t.values)
            print(
                "##########################################################################"
            )
            print('Creating Plot:', var, 'at time:', date)
            if df is not None:
                odf = df.loc[df.time == pd.Timestamp(t.values),
                             ['latitude', 'longitude', obs_variable[index]]]
            else:
                odf = None
            name = "{}.{}".format(out_name, obj.name)
            make_spatial_plot(obj.sel(time=t), odf, name)
            # if df is not None:
            #     make_spatial_bias_plot(df, name)


# def make_boxplot_giorgi(paired_data, savename):
#     from monet.util.tools import get_giorgi_region_df as ggrd
#     from monet.plots import savefig
#     import seaborn as sns
#     df=paird_data.copy()
#     df=ggrd(df)
#     dfa=df.dropna(subset = ['aod_550nm', 'pm25aod550'])[
#         'aod_550nm', 'GIORGI_ACRO']
#     dfm= df.dropna(subset=['aod_550nm', 'pm25aod550'])[
#         'pm25aod550', 'GIORGI_ACRO']
#     dfa['Legend']= 'AERONET'
#     dfm['Legend']= 'FV3CHEM'
#     dfa.rename({'aod_550nm': 'AOD'}, axis=1, inplace=True)
#     dfm.rename({'pm25aod550': 'AOD'}, axis=1, inplace=True)
#     dfn= pd.concat([dfa, dfm], ignore_index=True)
#     f, ax= plt.subplots(figsize=(12, 7))
#     sns.boxplot(ax=ax, x='GIORGI_ACRO', y='AOD', hue='Legend', data=dfn)
#     sns.despine()
#     plt.tight_layout(pad=0)
#     savefig(savename, dpi=100)
#     plt.close()


def get_df_region(obj, region):
    from monet.util.tools import get_giorgi_region_df as ggrd
    if region.lower() == 'global':
        obj['GIORGI_ACRO'] = 'global'
        return obj
    else:
        obj = ggrd(region)
        return obj.loc[obj.GIORGI_ACRO == region.upper()]


def get_region(obj, region):
    from monet.util.tools import get_giorgi_region_bounds as ggrb
    if region.lower() == 'global':
        return obj
    else:
        # print(region.upper())
        latmin, lonmin, latmax, lonmax, acro = ggrb(acronym=region.upper())
        out = obj.monet.window(
            lat_min=latmin, lon_min=lonmin, lat_max=latmax, lon_max=lonmax)
        return out


if __name__ == '__main__':

    parser = ArgumentParser(
        description='Make Spatial Plots for each time step in files',
        formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-f', '--files', help='input file names', type=str, required=True)
    parser.add_argument(
        '-p',
        '--paired_data',
        help='paired data input file names',
        type=str,
        required=False)
    parser.add_argument(
        '-v',
        '--variable',
        nargs='+',
        help='variable name to plot',
        required=True)
    parser.add_argument(
        '-vp',
        '--obs_variable',
        nargs='+',
        help='input file names',
        required=False)
    parser.add_argument(
        '-r',
        '--region',
        help='GIORGI Region ACRONYM',
        nargs='+',
        required=False,
        default=['global'])
    parser.add_argument(
        '-n',
        '--output_name',
        help='Output base name',
        type=str,
        required=False,
        default='FV3CHEM')
    parser.add_argument(
        '-d',
        '--daily',
        help='Resample to daily average',
        type=bool,
        required=False,
        default=False)
    args = parser.parse_args()

    finput = args.files
    paired_data = args.paired_data
    variable = args.variable
    obs_variable = args.obs_variable
    out_name = args.output_name
    daily = args.daily
    region = str(args.region[0])

    # open fv3chem
    obj = open_fv3chem(finput)

    # get the region if specified
    # print('region', region)
    ds = get_region(obj, region)
    # print(ds)
    # load the paired data
    if paired_data is not None:
        df = load_paired_data(paired_data)
        df = get_df_region(df, region)  # only the correct region
    else:
        df = None
    if daily:
        ds = ds.resample(time='D').mean()
        outname = "{}.{}.{}".format(out_name, 'daily', region)
        if df is not None:
            df = df.resample('D').mean()
    else:
        outname = "{}.{}".format(out_name, region)

    # make the plots
    make_plots(ds, df, variable, obs_variable, outname)
