#!/usr/bin/env python
import numpy as np

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
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

def Radar_Hydrometer_Thompson(pressure,temp,temp_wb,reflectivity,nz):
    qrain=np.zeros(nz)
    qsnow=np.zeros(nz)
    # constant parameter
    min_ref=0.0
    max_ref_snow=28.0
    max_ref_rain=55.0
    n0r_mp=8.0e6
    rd=287.0
    am_s=0.069
    bm_s=2.0
    PI=3.1415926536
    rho_i=890.0
    rho_w=1000.0
    a_min=1.0e-5
    sa=[5.065339,-0.062659,-3.032362,0.029469,-0.000285,0.31255,0.000204,0.003199,0.0,-0.015952]
    sb=[0.476221,-0.015896,0.165977,0.007468,-0.000141,0.060366,0.000079,0.000594,0.0,-0.003577]
    istatus=0
    # parameters
    f=(0.176/0.93)*(6.0/PI)*(6.0/PI)*(am_s/rho_i)*(am_s/rho_i)
    cse=[bm_s+1.0,bm_s+2.0,bm_s*2.0]
    oams=1.0/am_s
    crg=[24.0,1.0,24.0,5040.0]
    am_r=PI*rho_w/6.0
    for i in range(nz):
        if reflectivity[i]>=min_ref:
            rho=pressure[i]*100/(rd*temp[i])
            tc=temp[i]-273.15
            if temp[i]<=273.15:
                rfract=0.0
            elif temp[i]>=273.15+5:
                rfract=1.0
            else:
                rfract=(temp[i]-273.15)/5
            zes=(10.0**(0.1*np.minimum(reflectivity[i],max_ref_snow)))*(1.0-rfract)*1.0e-18
            zer=(10.0**(0.1*np.minimum(reflectivity[i],max_ref_snow)))*rfract*1.0e-18
            tc0=np.minimum(-0.1,tc)
            if bm_s<1.999 or bm_s>2.001:
                break
            loga_=sa[0]+sa[1]*tc0+sa[2]*cse[2]+sa[3]*tc0*cse[2]+sa[4]*tc0*tc0+\
                    sa[5]*cse[2]*cse[2]+sa[6]*tc0*tc0*cse[2]+\
                    sa[7]*tc0*cse[2]*cse[2]+sa[8]*tc0*tc0*tc0+\
                    sa[9]*cse[2]*cse[2]*cse[2]
            a_=np.maximum(10.0**loga_,a_min)
            b_=sb[0]+sb[1]*tc0+sb[2]*cse[2]+sb[3]*tc0*cse[2]+\
                sb[4]*tc0*tc0+sb[5]*cse[2]*cse[2]+\
                sb[6]*tc0*tc0*cse[2]+sb[7]*tc0*cse[2]*cse[2]+\
                sb[8]*tc0*tc0*tc0+sb[9]*cse[2]*cse[2]*cse[2]
            qsnow[i]=((zes/(f*a_))**(1.0/b_))/(rho*oams)
            qrain[i]=n0r_mp*am_r*crg[2]/rho*(zer/(n0r_mp*crg[3]))**(4.0/7.0)
    return qrain,qsnow

def Radar_Hydrometer_Ferrier(pressure,temp,temp_wb,reflectivity,nz):
    trans_qrain=np.zeros(nz)
    trans_qsnow_wet=np.zeros(nz)
    trans_qsnow_dry=np.zeros(nz)
    trans_qgraupel=np.zeros(nz)
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
    for i in range(nz):
        rho=pressure[i]*100/(rd*temp[i])
        ze=10.0**(0.1*reflectivity[i])
        if temp[i]>=273.15+5:
            #print "rain"
            trans_qrain[i]=zerf*(ze**zepowf)/rho
        elif temp[i]<273.15+5:
            if temp_wb[i]>=273.15+1.3:
                #print "mix (rain & wet snow)"
                snow_frac=(273.15+5-temp[i])/5
                trans_qsnow_wet[i]=zesnegf*((ze*snow_frac)**zepowf)/rho
                trans_qrain[i]=zerf*((ze*(1-snow_frac))**zepowf)/rho
            else:
                #print "wet snow"
                trans_qsnow_wet[i]=zesnegf*(ze**zepowf)/rho
        elif 273.15+5>temp[i]>=273.15:
            if temp_wb[i]<273.15+1.3:
                #print "wet snow"
                trans_qsnow_wet[i]=zesnegf*(ze**zepowf)/rho
            else:
                #print "mix (rain & wet snow)"
                rain_frac=(-273.15+temp[i])/5
                trans_qsnow_wet[i]=zesnegf*((ze*(1-rain_frac))**zepowf)/rho
                trans_qrain[i]=zerf*((ze*rain_frac)**zepowf)/rho
        elif temp[i]<273.15:
            if temp_wb[i]>=273.15:
                #print "mix (wet & dry snow)"
                dry_snow_frac=(273.15-temp[i])
                wet_snow_frac=(273.15-temp_wb[i])
                trans_qsnow_wet[i]=zesnegf*((ze*(wet_snow_frac/(dry_snow_frac+wet_snow_frac)))**zepowf)/rho
                trans_qsnow_dry[i]=zesnegf*((ze*(dry_snow_frac/(dry_snow_frac+wet_snow_frac)))**zepowf)/rho
            elif reflectivity[i]<55:
                #print "dry snow"
                trans_qsnow_dry[i]=zesnegf*(ze**zepowf)/rho
            else:
                #print "graupel"
                trans_qgraupel[i]=zehf*(zeh**zehpowf)/rho
    return trans_qrain,trans_qsnow_dry+trans_qsnow_wet,trans_qgraupel
