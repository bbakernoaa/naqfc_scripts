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
from glob import glob


'''
Simple utility to convert nemsio file into a netCDF4 file
Utilizes mkgfsnemsioctl utility from NWPROD and CDO
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
    from dask import dataframe as dd
    df = dd.read_hdf(fname)
    return df


def make_spatial_plot(da, proj):
    cbar_kwargs = dict(orientation='horizontal', pad=0.1, aspect=30)
    ax = da.monet.quick_map(cbar_kwargs=cbar_kwargs,
                            map_kwarg={'states': True})


def scatter_obs(df):
    from matplotlib.pyplot import gcf
    cf = gcf()
    cbar = cf.get_axes()[1]


def spatial_scatter(df, col='OZONE'):
    from matplotlib.pyplot import gcf, scatter
    cf = gcf()
    cbar = cf.get_axes()[1]
    s = 25
    scatter(df.longitude, df.latitude, c=df[col], s=s)
#     from .colorbars import cmap_discretize
#     x, y = m(df.longitude.values, df.Latitude.values)
#     s = 20
#     if create_cbar:
#         if discrete:
#             cmap = cmap_discretize(cmap, ncolors)
#             # s = 20
#           if (type(plotargs(vmin)) == None) | (type(plotargs(vmax)) == None):
#                 plt.scatter(x, y, c=df['Obs'].values, **plotargs)
#             else:
#                 plt.scatter(x, y, c=df['Obs'].values, **plotargs)
#         else:
#             plt.scatter(x, y, c=df['Obs'].values, **plotargs)
#     else:
#         plt.scatter(x, y, c=df['Obs'].values, **plotargs)


def make_plots(finput, paired_data, variable, obs_variable, verbose, region, epa_region):

    # open the files
    f = open_cmaq(finput)
    # get map projection
    proj = map_projection(f)
    if paired_data is not None:
        df = load_paired_data(paired_data)
    # loop over varaible list
    for index, var in enumerate(variable):
        obj = f[var]
        # loop over time
        for t in obj.time:
            ax = make_spatial_plot(obj.sel(time=t), proj)
            if paired_data is not None:
                ov = obs_variable[[index, 'latitude', 'longitude'].loc[obs_variable.time == t]
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
                        type=bool, action='store_true' required=False)
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
