#!/usr/bin/env python

import os
import re
import time
import urllib
import urllib2
import argparse
import warnings

def Calculate_Date(start_date):
    if time.localtime().tm_isdst==0:
        timelag=-5
    else:
        timelag=-4
    sec_start_date=time.mktime(time.strptime(start_date,'%Y%m%d%H'))
    return sec_start_date-5*60,sec_start_date+5*60

def Download_NSSL_Mosiac(start_date,save_dir):
    os.system('mkdir -p '+save_dir+'/'+start_date)
    thresh_min,thresh_max=Calculate_Date(start_date)
    rootpage='https://mrms.ncep.noaa.gov/data/3DReflPlus/MergedReflectivityQC_'
    levels=['00.50','00.75','01.00','01.25','01.50','01.75','02.00','02.25',\
            '02.50','02.75','03.00','03.50','04.00','04.50','05.00','05.50',\
            '06.00','06.50','07.00','07.50','08.00','08.50','09.00','10.00',\
            '11.00','12.00','13.00','14.00','15.00','16.00','17.00','18.00','19.00']
    for cell_level in levels:
        webpage=sorted(list(set(re.findall('MRMS_MergedReflectivityQC_'+cell_level+'_\d{8}-\d{6}.grib2.gz',os.popen('curl '+rootpage+cell_level+'/').read()))))
        for cell in webpage:
            print cell
            date=time.mktime(time.strptime(re.findall('\d{8}-\d{6}',cell)[0],'%Y%m%d-%H%M%S'))
            if thresh_min<=date<=thresh_max and os.path.isfile(save_dir+'/'+start_date+'/'+cell[:-3])==False:
                os.system('wget -t 3 '+rootpage+cell_level+'/'+cell+' -O '+save_dir+'/'+start_date+'/'+cell)
                os.system('gunzip '+save_dir+'/'+start_date+'/'+cell)

"""test!!!test"""
warnings.simplefilter('ignore')
parser=argparse.ArgumentParser(description="Download NSSL Mosiac")
parser.add_argument('a0',metavar='|Start Date|',type=str,help='Init Dir\n')
args=parser.parse_args()
Download_NSSL_Mosiac(args.a0,'/scratch/users/qzhang/Raw_Data/RefD_Mosiac')

