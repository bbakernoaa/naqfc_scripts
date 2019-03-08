#!/usr/bin/env python

__author__ = 'Patrick Campbell'
__email__ = 'patrick.c.campbell@noaa.gov'
__license__ = 'GPL'

#Simple MONET utility to command line pair model vs. observations

import os
import subprocess
import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from distutils.spawn import find_executable
from glob import glob

import pandas as pd

import monet
from monet.util.tools import long_to_wide


def pair_point(da, df, sub_map, interp):
    dfpair = da.monet.combine_point(
        df, mapping_table=sub_map, method=interp, cleanup=True)
    #     print(dfpair)
    return dfpair


def get_aeronet(start, end):
    dates = pd.date_range(start=start, end=end, freq='H')
    df = monet.obs.aeronet.add_data(dates, freq='H')
    return df.dropna(subset=['latitude', 'longitude'])


def open_cmaq(finput, verbose=False):
    dset = monet.models.fv3chem.open_mfdataset(finput)
    return dset


if __name__ == '__main__':

    parser = ArgumentParser(
        description='pairs cmaq model data to aqs observations',
        formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        '-f',
        '--files',
        help='string input model file directory/names',
        type=str,
        required=True)
    #    parser.add_argument('-s', '--startdates',  help='string input start date for pairing YYYY-MM-DD HH:MM:SS', type=str, required=True)
    #    parser.add_argument('-e', '--enddates',    help='string input end date for pairing YYYY-MM-DD HH:MM:SS', type=str, required=True)
    parser.add_argument(
        '-x',
        '--species',
        help='string input for obs species-variables to pair',
        type=str,
        nargs='+',
        required=False,
        default=['OZONE', 'PM2.5'])
    parser.add_argument(
        '-o',
        '--output',
        help='string output path for paired dataframe, stats, plots',
        type=str,
        required=False,
        default='paired_output.hdf')
    parser.add_argument(
        '-p',
        '--path',
        help='string path to director of network observations',
        type=str,
        required=False,
        default='/data/aqf2/barryb/5xpm/AQS_DATA/')
    parser.add_argument(
        '-n',
        '--network',
        help='string input data network name: airnow, aqs',
        type=str,
        required=False,
        default='airnow')
    parser.add_argument(
        '-m',
        '--model',
        help='input model: cmaq, fv3, hysplit (not-ready), or camx (not-ready)',
        type=str,
        required=False,
        default='cmaq')
    parser.add_argument(
        '-i',
        '--interp',
        help=
        'xesmf interpolation scheme, bilinear, conservative, nearest_s2d, nearest_d2s, patch',
        type=str,
        required=False,
        default='bilinear')
    parser.add_argument(
        '-v',
        '--verbose',
        help='print debugging information',
        action='store_true',
        required=False)
    args = parser.parse_args()

    finput = args.files
    #    start   = args.startdates
    #    end     = args.enddates
    species = args.species
    output = args.output
    datapath = args.path
    network = args.network
    model = args.model
    interp = args.interp
    verbose = args.verbose
    #reads model output (cmaq default)

    if model == 'cmaq':
        da = open_cmaq(finput, verbose=verbose)
    else:
        print('Must enter cmaq model right now')
        raise RuntimeError

#retrieves data observations and formats pandas dataframe (airnow default)
    start = da.time.to_index().min() - pd.Timedelta(1, unit='h')
    end = da.time.to_index().max() + pd.Timedelta(1, unit='h')
    # print(start.strftime('%Y%m%d%H'), end.strftime('%Y%m%d%H'))
    # print(start)
    # print(end)
    df = get_aeronet(start, end)
    #pairs surface point-type observations with 2D model parameters

    mapping_table = {
        'pm25aod550': 'aod_550nm',
        'dust25aod550': 'aod_550nm',
        'salt25aod550': 'aod_550nm',
        'sulf25aod550': 'aod_550nm',
        'oc25aod550': 'aod_550nm',
        'bc25aod550': 'aod_550nm'
    }
    dfpair = pair_point(da, df, mapping_table, interp)
    # dfpair = dfpair.dropna(subset=['aod_550nm', 'pm25aod550'])
    print(dfpair.head())
    dfpair.to_hdf(output, 'df', format='table')

    sys.exit(0)
