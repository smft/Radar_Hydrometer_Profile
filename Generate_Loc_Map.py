#!/usr/bin/env python

import glob
import math
import pygrib
import numpy as np
import cPickle as pickle
import multiprocessing as mp
from ctypes import *
from netCDF4 import Dataset
from scipy.spatial import cKDTree

def make_shared_array(none_shared_array):
    shape=np.shape(none_shared_array)
    shared_array_input_base=mp.Array(c_double,np.prod(shape))
    data=np.ctypeslib.as_array(shared_array_input_base.get_obj())
    data=data.reshape(list(shape))+none_shared_array
    return data

def distance(origin,destination):
    lat1=origin[0]
    lon1=origin[1]
    lat2=destination[0]
    lon2=destination[1]
    radius=6371  # earth radius
    dlat=math.radians(lat2-lat1)
    dlon=math.radians(lon2-lon1)
    a=(math.sin(dlat/2)*math.sin(dlat/2)+math.cos(math.radians(lat1))*math.cos(math.radians(lat2))*\
        math.sin(dlon/2)*math.sin(dlon/2))
    c=2*math.atan2(math.sqrt(a),math.sqrt(1-a))
    d=radius*c
    return d

def ckd_search(radar_loc_1d,model_loc_1d,radius,idx,q_):
    # radar_loc_1d:[lat,lon]; model_loc_1d:[lat,lon]
    obs_tree=cKDTree(radar_loc_1d)
    print "Start"
    index_inuse=dict()
    for i in idx:
        dist,radar_idx=obs_tree.query(model_loc_1d[i,:],k=30)
        trans=list()
        if len(dist)>1:
            for cell in radar_idx:
                dist_km=distance(model_loc_1d[i,:],radar_loc_1d[cell,:])
                if dist_km<=radius:
                    trans+=[[cell,dist_km]]
            if len(trans)!=0:
                index_inuse[i]=trans
            print i,len(trans)
        else:
            dist_km=distance(model_loc_1d[i,:],radar_loc_1d[idx])
            if dist_km<=radius:
                trans+=[[idx,dist_km]]
            if len(trans)!=0:
                index_inuse[i]=trans
            print i,len(trans)
    q_.put(index_inuse)

def Interpolation(radar_lat,radar_lon,model_lat,model_lon,radius):
    radar_ny,radar_nx=np.shape(radar_lat)
    model_ny,model_nx=np.shape(model_lat)
    radar_loc=list()
    for i in range(radar_ny):
        for j in range(radar_nx):
            radar_loc+=[[radar_lat[i,j],radar_lon[i,j]]]
    radar_loc_share=make_shared_array(np.asarray(radar_loc))
    del radar_loc
    model_loc=list()
    for i in range(model_ny):
        for j in range(model_nx):
            model_loc+=[[model_lat[i,j],model_lon[i,j]]]
    model_loc_share=make_shared_array(np.asarray(model_loc))
    del model_loc
    cpu_count=mp.cpu_count()
    q_=mp.Queue(cpu_count)
    for cell in np.array_split(np.arange(0,model_ny*model_nx,1),cpu_count):
        p=mp.Process(target=ckd_search,args=(radar_loc_share,model_loc_share,radius,cell,q_))
        p.start()
    rslt=dict()
    for i in range(cpu_count):
        trans=q_.get()
        if len(trans.keys())>0:
            rslt.update(trans)
    return rslt

def Read_NSSL_Radar(filename):
    flag=pygrib.open(filename)
    flag.seek(0)
    grb=flag.read(1)[0]
    data,lat,lon=grb.data()
    data[data<0]=np.nan
    return data,lat,lon

def Read_WRF_Domain(filename):
    flag=Dataset(filename)
    lat=flag.variables['XLAT_M'][0,:,:]
    lon=flag.variables['XLONG_M'][0,:,:]
    print np.nanmin(lon)
    lon[lon<0]+=360
    return lat,lon

"""test!!!test"""
dbz,radar_lat,radar_lon=Read_NSSL_Radar('/scratch/users/qzhang/Raw_Data/RefD_Mosiac/2020040717/MRMS_MergedReflectivityQC_05.00_20200407-170234.grib2')
wrf_lat,wrf_lon=Read_WRF_Domain('/data/users/qzhang/HWT_3km_System/Radar_Process/Grid/geo_em.d01.nc')
loc_map=Interpolation(radar_lat,radar_lon,wrf_lat,wrf_lon,2)
flag_save=open('/data/users/qzhang/HWT_3km_System/Radar_Process/Grid/Loc_Map.pickle','wb')
pickle.dump(loc_map,flag_save)
flag_save.close()

