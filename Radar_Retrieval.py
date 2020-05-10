#!/usr/bin/env python
import wrf
import time
import zlib
import warnings
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import cPickle as pickle
import matplotlib.pyplot as plt
from mpi4py import MPI
from netCDF4 import Dataset

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# var                   description
# temp                  temperature in K
# pressure              pressure in hPa
# temp_wb               wet_bulb temperature in K
# reflectivity          radar echo in model levels
# cloud phase           0: clear
#                       1: rain
#                       2: snow
#                       3: freeze rain
#                       4: wet snow (sleet)
#                       5: graupel
# psfc                  surface pressure in Pa
# rh                    relative humidity from 0 to 100
# watervapor            mixing ratio in g/kg
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

def tmr(watervapor,pressure):
#!
#!   this function returns the temperature (celsius) on a mixing
#!   ratio line w (g/kg) at pressure p (mb). the formula is given in
#!   table 1 on page 7 of stipanuk (1973).
#!
    c1=.0498646455;c2=2.4082965;c3=7.07475
    c4=38.9114;c5=.0915;c6=1.2035
    y=watervapor*pressure/(622.+watervapor)
    x=np.log10(y)
    tmrk= 10.**(c1*x+c2)-c3+c4*((10.**(c5*x)-c6)**2.)
    return tmrk-273.15

def tda(temp,pressure):
#!
#!   this function returns the temperature tda (celsius) on a dry adiabat
#!   at pressure p (millibars). the dry adiabat is given by
#!   potential temperature o (celsius). the computation is based on
#!   poisson's equation.
#!
    return (temp+273.15)*((pressure*.001)**.286)-273.15

def tsa(os,pressure):
#!
#!   this function returns the temperature tsa (celsius) on a saturation
#!   adiabat at pressure p (millibars). os is the equivalent potential
#!   temperature of the parcel (celsius). sign(a,b) replaces the
#!   algebraic sign of a with that of b.
#!   b is an empirical constant approximately equal to 0.001 of the latent
#!   heat of vaporization for water divided by the specific heat at constant
#!   pressure for dry air.
#!
    b=2.6518986
    a=os+273.15
    tq=253.15
    d=120.
#!   iterate to obtain sufficient accuracy....see table 1, p.8
#!   of stipanuk (1973) for equation used in iteration.
    for i in range(1,13):
        tqk=tq-273.15
        d=d/2.
        x=a*np.exp(-b*w(tqk,pressure)/tq)-tq*((1000./pressure)**.286)
        if np.abs(x)<1e-7:
            tsa_=tq-273.15
            break;
        else:
            if x>=0:
                tq+=np.abs(d)
            else:
                tq-=np.abs(d)
    tsa_=tq-273.15
    return tsa_

def esat(temp):
#!
#!   this function returns the saturation vapor pressure over
#!   water (mb) given the temperature (celsius).
#!   the algorithm is due to nordquist, w.s.,1973: "numerical approxima-
#!   tions of selected meteorlolgical parameters for cloud physics prob-
#!   lems," ecom-5475, atmospheric sciences laboratory, u.s. army
#!   electronics command, white sands missile range, new mexico 88002.
#!
    tk = temp+273.15
    p1 = 11.344-0.0303998*tk
    p2 = 3.49149-1302.8844/tk
    c1 = 23.832241-5.02808*np.log10(tk)
    return 10.**(c1-1.3816e-7*10.**p1+8.1328e-3*10.**p2-2949.076/tk)

def w(temp,pressure):
#!
#!  this function returns the mixing ratio (grams of water vapor per
#!  kilogram of dry air) given the dew point (celsius) and pressure
#!  (millibars). if the temperture  is input instead of the
#!  dew point, then saturation mixing ratio (same units) is returned.
#!  the formula is found in most meteorological texts.
#!
  return 622.*esat(temp)/(pressure-esat(temp))

def Calculate_Qvapor(temp,rh):
    t00=273.15
    a=10.79574
    b=-5.02800
    c=1.50475e-4
    d=-8.2969
    e=0.42873e-3
    f=4.76955
    g=0.78614
    tmp_above_freezing=a*(1-t00/temp)+\
                        b*np.log10(temp/t00)+\
                        c*(1-10**(d*(temp/t00-1)))+\
                        e*(10**(f*(1-t00/temp))-1)+g
    a=-9.09685
    b=-3.56654
    c=0.87682
    d=0.78614
    tmp_below_freezing=a*(t00/temp-1)+\
                        b*np.log10(t00/temp)+\
                        c*(1-temp/t00)+0.78614
    if temp>=273.15:
        return (10**tmp_above_freezing)*(rh/100.0)
    else:
        return (10**tmp_below_freezing)*(rh/100.0)

def Calculate_Wetbulb(temp,pressure,rh,qvapor):
    watervapor=Calculate_Qvapor(temp,rh)
    ao=temp*((1000./pressure)**.286)-273.15
    pi=pressure.copy()
    for i in range(1,11):
        x=0.02*(tmr(watervapor,pi)-tda(ao,pi))
        if abs(x)<0.01:
            break
        pi=pi*(2.0**x)
    ti=tda(ao,pi)
    aos=(ti+273.15)*((1000./pi)**.286)*(np.exp(2.6518986*watervapor/(ti+273.15)))-273.15
    temp_wb=tsa(aos,pressure)
    return temp_wb+273.15,watervapor-qvapor

def Radar_Hydrometer_Ferrier(pressure,temp_,watervapor,qrain,qsnow,qgraupel,reflectivity):
    rmse_list=list()
    rtv_qvapor=list()
    rtv_qrain=list()
    rtv_qsnow=list()
    rtv_graupel=list()
    for delta_rh in np.arange(-5,5,0.1):
        trans_qrain=0
        trans_qsnow_wet=0
        trans_qsnow_dry=0
        trans_qgraupel=0
        ki2=0.176
        kw2=0.93
        m3todBZ=1.0e+18
        Zefact=720.0
        lg10div=0.10
        pi=3.1415926
        N0r=8.0E+06
        N0s=3.0E+06
        N0h=4.0E+04
        N0xpowf=3.0/7.0
        K2powf=4.0/7.0
        zkpowf=4.0/7.0
        zepowf=4.0/7.0
        zehpowf=(4.0/7.0)*1.0526
        rhoi=917.0
        rhor=1000.0
        rhos=100.0
        rhoh=913.0
        rhoipowf=8.0/7.0
        rhospowf=1.0/7.0
        rd=287.0
        thresh_ref=0.0
        zkconst=(Zefact*m3todBZ)**zkpowf
        zerf=1000*(pi*(N0r**N0xpowf)*rhor)/zkconst
        zesnegf=1000*(pi*(kw2**K2powf)*(N0s**N0xpowf)*(rhoi**rhoipowf))/(zkconst*(ki2**K2powf)*(rhos**rhospowf))
        zesposf=1000*(pi*(N0s**N0xpowf)*rhos)/zkconst
        zehf=1000*(pi*(N0h**N0xpowf)*rhoh)/zkconst
        temp=temp_
        rho=pressure*100/(rd*temp)
        ze=10.0**(0.1*reflectivity)
        temp_wb,qvapor_bias=Calculate_Wetbulb(temp,pressure,95+delta_rh,watervapor)
        if temp>=273.15+5:
            #print "rain"
            trans_qrain=zerf*(ze**zepowf)/rho
        elif temp<273.15+5:
            if temp_wb>=273.15+1.3:
                #print "mix (rain & wet snow)"
                snow_frac=(273.15+5-temp)/5
                trans_qsnow_wet=zesnegf*((ze*snow_frac)**zepowf)/rho
                trans_qrain=zerf*((ze*(1-snow_frac))**zepowf)/rho
            else:
                #print "wet snow"
                trans_qsnow_wet=zesnegf*(ze**zepowf)/rho
        elif 273.15+5>temp>=273.15:
            if temp_wb<273.15+1.3:
                #print "wet snow"
                trans_qsnow_wet=zesnegf*(ze**zepowf)/rho
            else:
                #print "mix (rain & wet snow)"
                rain_frac=(-273.15+temp)/5
                trans_qsnow_wet=zesnegf*((ze*(1-rain_frac))**zepowf)/rho
                trans_qrain=zerf*((ze*rain_frac)**zepowf)/rho
        elif temp<273.15:
            if temp_wb>=273.15:
                #print "mix (wet & dry snow)"
                dry_snow_frac=(273.15-temp)
                wet_snow_frac=(273.15-temp_wb)
                trans_qsnow_wet=zesnegf*((ze*(wet_snow_frac/(dry_snow_frac+wet_snow_frac)))**zepowf)/rho
                trans_qsnow_dry=zesnegf*((ze*(dry_snow_frac/(dry_snow_frac+wet_snow_frac)))**zepowf)/rho
            elif reflectivity<55:
                #print "dry snow"
                trans_qsnow_dry=zesnegf*(ze**zepowf)/rho
            else:
                #print "graupel"
                trans_qgraupel=zehf*(zeh**zehpowf)/rho
        if np.isnan(trans_qgraupel):
            trans_qgraupel=0
        if np.isnan(trans_qsnow_dry):
            trans_qsnow_dry=0
        if np.isnan(trans_qsnow_wet):
            trans_qsnow_wet=0
        if np.isnan(trans_qrain):
            trans_qrain=0
        rmse_list+=[np.abs(trans_qrain-qrain)+np.abs(trans_qsnow_dry+trans_qsnow_wet-qsnow)+\
                    np.abs(trans_qgraupel-qgraupel)+np.abs(qvapor_bias)]
        rtv_qrain+=[trans_qrain];rtv_qsnow+=[trans_qsnow_dry+trans_qsnow_wet];rtv_graupel+=[trans_qgraupel]
        rtv_qvapor+=[watervapor+qvapor_bias]
    rmse_list=np.asarray(rmse_list)
    idx=np.unravel_index(rmse_list.argmin(),rmse_list.shape)[0]
    if abs(qvapor_bias/watervapor)<1e9:
        return rtv_qrain[idx],rtv_qsnow[idx],rtv_graupel[idx],rtv_qvapor[idx]
    else:
        return np.nan,np.nan,np.nan,np.nan

def Refd_DA(data,rank):
    rslt=list()
    refd_height=np.array([00.50,00.75,01.00,01.25,01.50,01.75,02.00,02.25,\
                        02.50,02.75,03.00,03.50,04.00,04.50,05.00,05.50,\
                        06.00,06.50,07.00,07.50,08.00,08.50,09.00,10.00,\
                        11.00,12.00,13.00,14.00,15.00,16.00,17.00,18.00])*1000
    for cell in data:
        ts=time.time()
        refd_interp=np.interp(cell['height'],refd_height,cell['reflectivity'],left=0,right=0)
        nz=len(refd_interp)
        for i in range(nz):
            if refd_interp[i]>0:
                rain,snow,graup,watervapor_=Radar_Hydrometer_Ferrier(cell['pressure'][i],cell['temp'][i],\
                                                    cell['watervapor'][i],cell['qrain'][i],cell['qsnow'][i],\
                                                    cell['qgraupel'][i],refd_interp[i])
                if np.isnan(rain)==False and np.isnan(snow)==False and \
                    np.isnan(graup)==False and np.isnan(watervapor_)==False:
                    rslt+=[[i,cell['idx_ny'],cell['idx_nx'],rain,snow,graup,watervapor_]]
        #print "Rank:",("%04d" % rank),"Index:",("%05d" % cell['idx_ny']),("%05d" % cell['idx_nx']),\
        #        "Time Usage:",('%09.5f' % (time.time()-ts))
    return rslt

"""test!!!test"""

comm=MPI.COMM_WORLD
size=comm.Get_size()
rank=comm.Get_rank()
if rank==0:
    warnings.simplefilter("ignore")
    parser=argparse.ArgumentParser(description="Download NSSL Mosiac")
    parser.add_argument('a0',metavar='|Initial Dir|',type=str,help='Init Dir\n')
    parser.add_argument('a1',metavar='|Initial Date|',type=str,help='Init Date\n')
    args=parser.parse_args()
    flag=Dataset(args.a0+'/wrfinput_d01.'+args.a1[:4]+'-'+args.a1[4:6]+'-'+args.a1[6:8]+'_'+args.a1[8:]+':00:00','r+')
    pressure=wrf.getvar(flag,'pressure',timeidx=0,meta=False)
    temp=wrf.getvar(flag,'tk',timeidx=0,meta=False)
    watervapor=flag.variables['QVAPOR'][0,:,:,:]*1000
    qrain=flag.variables['QRAIN'][0,:,:,:]*1000
    qsnow=flag.variables['QSNOW'][0,:,:,:]*1000
    qgraupel=flag.variables['QGRAUP'][0,:,:,:]*1000
    height=wrf.getvar(flag,'height',timeidx=0,meta=False)
    reflectivity=pickle.load(open('/scratch/users/qzhang/Raw_Data/RefD_Mosiac/'+args.a1+'/RefD.pickle'))
    reflectivity[reflectivity<0]=0
    nz,ny,nx=np.shape(watervapor)
    data=list()
    rslt=list()
    for i in range(ny):
        for j in range(nx):
            if np.nansum(reflectivity[:,i,j])>0:
                data+=[{'idx_ny':i,'idx_nx':j,\
                    'temp':temp[:,i,j],\
                    'pressure':pressure[:,i,j],\
                    'watervapor':watervapor[:,i,j],\
                    'qrain':qrain[:,i,j],\
                    'qsnow':qsnow[:,i,j],\
                    'qgraupel':qgraupel[:,i,j],\
                    'height':height[:,i,j],\
                    'reflectivity':reflectivity[:,i,j]}]
    rank_counts=np.linspace(0,len(data),size).astype(int)
    print 'Retrievals to be Calculated:',len(data)
    for i in range(1,size):
        if i!=size-1:
            comm.send(data[rank_counts[i-1]:rank_counts[i]],dest=i)
        else:
            comm.send(data[rank_counts[i-1]:],dest=i)
else:
    warnings.simplefilter("ignore")
    cell_data=comm.recv()
    ts=time.time()
    rslt=Refd_DA(cell_data,rank)
    print "CPU: ",('%04d' % rank),\
            "Finishes ",('%09d' % len(cell_data)),\
            "Reflectivity Profiles",\
            "Time Use: ",('%09.3f' % (time.time()-ts))
rslt=comm.gather(rslt,root=0)
if rank==0:
    rtv_qvapor=np.zeros([nz,ny,nx])
    rtv_qrain=np.zeros([nz,ny,nx])
    rtv_qsnow=np.zeros([nz,ny,nx])
    rtv_qgraupel=np.zeros([nz,ny,nx])
    min_graup=np.min(qgraupel)
    min_rain=np.min(qrain)
    min_snow=np.min(qsnow)
    for cell in rslt:
        for idx in cell:
            i=idx[0];j=idx[1];k=idx[2]
            rtv_qgraupel[i,j,k]=idx[-2] if qgraupel[i,j,k]==min_graup else qgraupel[i,j,k]
            rtv_qsnow[i,j,k]=idx[-3] if qsnow[i,j,k]==min_snow else qsnow[i,j,k]
            rtv_qrain[i,j,k]=idx[-4] if qrain[i,j,k]==min_rain else qrain[i,j,k]
    flag.variables['QRAIN'][0,:,:,:]=rtv_qrain/1000.0
    flag.variables['QSNOW'][0,:,:,:]=rtv_qsnow/1000.0
    flag.variables['QGRAUP'][0,:,:,:]=rtv_qgraupel/1000.0
    #flag.variables['QVAPOR'][0,:,:,:]=rtv_qvapor/1000.0
    flag.close()
    #flag_save=open(args.a0+'/hydro_rtv.pickle','wb')
    #pickle.dump({'QRAIN':zlib.compress(pickle.dumps(rtv_qvapor)),\
    #            'QGRAUP':zlib.compress(pickle.dumps(rtv_qgraupel)),\
    #            'QSNOW':zlib.compress(pickle.dumps(rtv_qsnow)),\
    #            'QRAIN':zlib.compress(pickle.dumps(rtv_qrain))},flag_save)
    #flag_save.close()
    #qvapor_bias=rtv_qvapor-watervapor;qvapor_bias[rtv_qvapor==0]=np.nan
    #qrain_bias=rtv_qrain-qrain;qrain_bias[rtv_qrain==0]=np.nan
    #qsnow_bias=rtv_qsnow-qsnow;qsnow_bias[rtv_qsnow==0]=np.nan
    #qgraupel_bias=rtv_qgraupel-qgraupel;qgraupel_bias[rtv_qgraupel==0]=np.nan
    #print "Watervapor Bias:\n",np.nanmean(np.nanmean(qvapor_bias,axis=-1),axis=-1)
    #print "Rain Bias:\n",np.nanmean(np.nanmean(qrain_bias,axis=-1),axis=-1)
    #print "Snow Bias:\n",np.nanmean(np.nanmean(qsnow_bias,axis=-1),axis=-1)
    #print "Graupel Bias:\n",np.nanmean(np.nanmean(qgraupel_bias,axis=-1),axis=-1)
    #for i in range(nz):
    #    plt.subplot(2,3,1);plt.imshow(rtv_qrain[i,:,:],vmin=0,vmax=2,cmap='rainbow');plt.colorbar()
    #    plt.subplot(2,3,2);plt.imshow(qrain[i,:,:],vmin=0,vmax=2,cmap='rainbow');plt.colorbar()
    #    plt.subplot(2,3,3);plt.imshow(qrain[i,:,:]-rtv_qrain[i,:,:],vmin=-2,vmax=2,cmap='bwr');plt.colorbar()
    #    plt.subplot(2,3,4);plt.imshow(rtv_qsnow[i,:,:],vmin=0,vmax=2,cmap='rainbow');plt.colorbar()
    #    plt.subplot(2,3,5);plt.imshow(qsnow[i,:,:],vmin=0,vmax=2,cmap='rainbow');plt.colorbar()
    #    plt.subplot(2,3,6);plt.imshow(qsnow[i,:,:]-rtv_qsnow[i,:,:],vmin=-2,vmax=2,cmap='bwr');plt.colorbar()
    #    plt.savefig(args.a0+'/'+("%02d" % i)+'.png')
    #    plt.close()
