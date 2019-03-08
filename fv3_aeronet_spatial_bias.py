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


# def open_fv3chem(finput):
#     from monet.models import fv3chem
#     f = fv3chem.open_mfdataset(finput)
#     return f


def load_paired_data(fname):
    return pd.read_hdf(fname)


def make_spatial_bias_plot(df,
                           out_name,
                           col1='aod_550nm',
                           col2='pm25aod550',
                           date=None,
                           **kwargs):
    ax = monet.plots.sp_scatter_bias(df, col1=col1, col2=col2, **kwargs)
    # date = df.time.min()
    date = pd.Timestamp(date)
    plt.title(date.strftime('time=%Y/%m/%d %H:00 | FV3 - AERONET (AOD)'))
    plt.tight_layout(pad=0)
    savename = "{}.{}".format(out_name, date.strftime('sb.%Y%m%d%H.jpg'))
    print(savename)
    monet.plots.savefig(savename, bbox_inches='tight', dpi=100, decorate=True)
    plt.close()


def make_plots(df, variable, obs_variable, out_name):
    # loop over varaible list
    print(variable, obs_variable)
    for obsv, v in zip(obs_variable, variable):
        # loop over time
        if 'global' in out_name:
            name = "{}.{}".format(out_name, v)
            print(df.keys())
            odf = df.dropna(subset=[obsv, v])
        for t in odf.time.unique():
            date = pd.Timestamp(t)
            print(
                ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            print('Creating Plot:', v, 'at time:', date)
            print(
                ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            odf = odf.loc[df.time ==
                          date, ['time', 'latitude', 'longitude', obsv, v]]
            name = "{}.{}".format(out_name, v)
            print(t)
            if ~odf.empty:
                make_spatial_bias_plot(
                    odf,
                    name,
                    col1=obsv,
                    col2=v,
                    date=t,
                    cmap='RdBu_r',
                    edgecolor='k',
                    linewidth=.8)


def get_df_region(obj, region):
    from monet.util.tools import get_giorgi_region_df as ggrd
    if region.lower() == 'global':
        obj['GIORGI_ACRO'] = 'global'
        return obj
    else:
        obj = ggrd(region)
        return obj.loc[obj.GIORGI_ACRO == region.upper()]


if __name__ == '__main__':

    parser = ArgumentParser(
        description='Make Spatial Plots for each time step in files',
        formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-p',
        '--paired_data',
        help='paired data input file names',
        type=str,
        required=True)
    parser.add_argument(
        '-v', '--variable', nargs='+', help='Varaible 1', required=True)
    parser.add_argument(
        '-vp', '--obs_variable', nargs='+', help='Variable 2', required=False)
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

    paired_data = args.paired_data
    variable = args.variable
    obs_variable = args.obs_variable
    out_name = args.output_name
    daily = args.daily
    region = str(args.region[0])

    df = load_paired_data(paired_data)
    if daily:
        df.index = df.time
        dfnew = df.resample('D').mean()
        outname = "{}.{}.{}".format(out_name, 'daily', region)
    else:
        dfnew = df
        outname = "{}.{}".format(out_name, region)

    dfnew = get_df_region(dfnew, region)  # only the correct region

    # make the plots
    make_plots(df, variable, obs_variable, outname)
