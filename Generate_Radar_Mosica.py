#!/usr/bin/env python

import os
import re
import wrf
import zlib
import time
import glob
import pygrib
import ncepbufr
import argparse
import warnings
import numpy as np
import multiprocessing as mp
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cPickle as pickle
from ctypes import *
from netCDF4 import *
from libs.Radar_Retrieval import *
from scipy.interpolate import interp1d

def make_shared_array(none_shared_array):
    shape=np.shape(none_shared_array)
    shared_array_input_base=mp.Array(c_double,np.prod(shape))
    data=np.ctypeslib.as_array(shared_array_input_base.get_obj())
    data=data.reshape(list(shape))+none_shared_array
    return data

def Read_NSSL_Radar(filename):
    flag=pygrib.open(filename)
    flag.seek(0)
    grb=flag.read(1)[0]
    data,lat,lon=grb.data()
    data[data<0]=np.nan
    return data,lat,lon

def Load_Loc_Map(filename):
    return pickle.load(open(filename))

def interpolate(AGL_Level,NSSL_Mosiac_Dir,loc_map,model_ny,model_nx):
    dbz=list()
    for cell in glob.glob(NSSL_Mosiac_Dir+'/MRMS_MergedReflectivityQC_'+AGL_Level+'_*.grib2'):
        try:
            dbz+=[Read_NSSL_Radar(cell)[0]]
        except:
            pass
    dbz=np.asarray(dbz)
    nssl_nt,nssl_ny,nssl_nx=np.shape(dbz)
    if nssl_nt>0:
        model_mosiac=np.zeros([model_ny,model_nx])
        for i in range(model_ny):
            for j in range(model_nx):
                if loc_map.has_key(i*model_nx+j):
                    nssl_loc_idx=np.asarray(loc_map[i*model_nx+j])
                    nssl_y=list(nssl_loc_idx[:,0].astype(np.int)/nssl_nx)
                    nssl_x=list(nssl_loc_idx[:,0].astype(np.int)%nssl_nx)
                    trans_dbz=dbz[:,nssl_y,nssl_x]
                    trans_mask=list()
                    for cell in trans_dbz:
                        trans=np.zeros_like(nssl_loc_idx[:,1])
                        trans[np.isnan(cell)==False]=1
                        trans_mask+=[trans]
                    trans_mask=np.asarray(trans_mask)
                    model_mosiac[i,j]=np.nansum(trans_dbz/nssl_loc_idx[:,1])/np.nansum(trans_mask/nssl_loc_idx[:,1])
    return model_mosiac

def Interpolate(AGL_Levels,NSSL_Mosiac_Dir,loc_map,model_ny,model_nx,q_):
    rslt=dict()
    for AGL_Level in AGL_Levels:
        ts=time.time()
        rslt[AGL_Level]=interpolate(AGL_Level,NSSL_Mosiac_Dir,loc_map,model_ny,model_nx)
        print "Finish AGL Level @ "+AGL_Level+"km, Time Usage : "+('%09.3f' % (time.time()-ts))
    q_.put(rslt)

def Calculate_Hydrometer_RTV(pressure,height,temp,temp_wb,reflectivity,idx_y,idx_x,nz,q_):
    refd_height=np.array([00.50,00.75,01.00,01.25,01.50,01.75,02.00,02.25,\
                        02.50,02.75,03.00,03.50,04.00,04.50,05.00,05.50,\
                        06.00,06.50,07.00,07.50,08.00,08.50,09.00,10.00,\
                        11.00,12.00,13.00,14.00,15.00,16.00,17.00,18.00])*1000
    trans_total=list()
    ts=time.time()
    for i in idx_y:
        for j in idx_x:
            refd_interp=np.interp(height[:,i,j],refd_height,reflectivity[:,i,j],left=0,right=0)
            #trans_rain,trans_snow=Radar_Hydrometer_Thompson(pressure[:,i,j],temp[:,i,j],temp_wb[:,i,j],refd_interp,nz)
            trans_rain,trans_snow,trans_graup=Radar_Hydrometer_Ferrier(pressure[:,i,j],temp[:,i,j],temp_wb[:,i,j],refd_interp,nz)
            trans_total+=[[i,j,trans_rain,trans_snow,trans_graup]]
    q_.put(trans_total)

def write_radarbufr(bufr,ny,nx,data,obs_time):
    hdrstr='SID XOB YOB DHR TYP'
    obsstr='HREF'
    idate=int(obs_time)
    subset='ADPUPA'
    ludx=22
    lendian_in=10
    for i in range(nx):
        for j in range(ny):
            bufr.open_message(subset,idate)
            hdr=bufr.missing_value*np.ones(len(hdrstr.split(' ')),np.float)
            obs=bufr.missing_value*np.ones((len(obsstr.split(' ')),31),np.float)
            hdr[0]=np.fromstring(('%08d' % 99999999),dtype=np.float)[0]
            hdr[1]=float(nx+1)
            hdr[2]=float(ny+1)
            hdr[3]=0.0
            hdr[4]=500
            obs[:]=data[:,j,i]
            bufr.write_subset(hdr,hdrstr)
            bufr.write_subset(obs,obsstr,end=True)
            bufr.close_message()

"""test!!!test"""
warnings.simplefilter('ignore')
parser=argparse.ArgumentParser(description="Convert NSSL Mosiac")
parser.add_argument('a0',metavar='|Start Date|',type=str,help='Start Date\n')
parser.add_argument('a1',metavar='|Init Dir|',type=str,help='Init Date\n')
args=parser.parse_args()

NSSL_Mosiac_Dir='/scratch/users/qzhang/Raw_Data/RefD_Mosiac/'+args.a0
flag_Init=Dataset(args.a1+'/wrfinput_d01.'+args.a0[:4]+'-'+args.a0[4:6]+'-'+args.a0[6:8]+'_'+args.a0[8:]+':00:00')
lat=flag_Init.variables['XLAT'][0,:,:]
model_ny,model_nx=np.shape(lat)
flag_Init.close()
AGL_Levels=np.asarray(['00.50','00.75','01.00','01.25','01.50','01.75','02.00','02.25',\
                        '02.50','02.75','03.00','03.50','04.00','04.50','05.00','05.50',\
                        '06.00','06.50','07.00','07.50','08.00','08.50','09.00','10.00',\
                        '11.00','12.00','13.00','14.00','15.00','16.00','17.00','18.00'])
cpu_count=mp.cpu_count()
if cpu_count>len(AGL_Levels):
    cpu_count=len(AGL_Levels)
q_=mp.Queue(cpu_count)
M=mp.Manager()
loc_map=M.dict()
loc_map=Load_Loc_Map('/data/users/qzhang/HWT_3km_System/Radar_Process/Grid/Loc_Map.pickle')
for cell in np.array_split(AGL_Levels,cpu_count):
    p=mp.Process(target=Interpolate,args=(cell,NSSL_Mosiac_Dir,loc_map,model_ny,model_nx,q_))
    p.start()
rslt=dict()
for i in range(cpu_count):
    rslt.update(q_.get())
q_.close()
radar_mosiac=list()
for cell in AGL_Levels:
    radar_mosiac+=[rslt[cell]]
radar_mosiac=np.asarray(radar_mosiac)
flag_save=open(NSSL_Mosiac_Dir+'/RefD.pickle','wb')
pickle.dump(radar_mosiac,flag_save)
flag_save.close()
#os.system('rm -f '+NSSL_Mosiac_Dir+'/MRMS_MergedReflectivityQC_*')
