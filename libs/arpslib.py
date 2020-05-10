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
    return temp*((pressure*.001)**.286)-273.15

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
            tsa=tq-273.15
            break;
        if x>=0:
            tq+=np.abs(d)
        else:
            tq-=np.abs(d)
    tsa=tq-273.15
    return tsa

def esat(temp):
#!
#!   this function returns the saturation vapor pressure over
#!   water (mb) given the temperature (celsius).
#!   the algorithm is due to nordquist, w.s.,1973: "numerical approxima-
#!   tions of selected meteorlolgical parameters for cloud physics prob-
#!   lems," ecom-5475, atmospheric sciences laboratory, u.s. army
#!   electronics command, white sands missile range, new mexico 88002.
#!
    tk = t
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
  return 622.*esat(t)/(p-esat(t))

def Calculate_Wetbulb(temp,pressure,watervapor):
    ao=temp*((1000./p)**.286)-273.15
    pi=pressure.copy()
    for i in range(1,11):
        x=0.02*(tmr(watervapor,pi)-tda(ao,pi))
        if abs(x)<0.01:
            break
        pi=pi*(2.0**x)
    ti=tda(ao,pi)
    aos=(ti+273.15)*((1000./pi)**.286)*(exp(2.6518986*watervapor/(ti+273.15)))-273.15
    temp_wb=tsa(aos,p)
    return temp_wb+273.15
