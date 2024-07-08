# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 08:52:19 2023

@author: equin
"""

import time
from nrlmsise_00_header import *
from nrlmsise_00 import *


def cal_gtd7(alt,doy):
    """
    This function is based on the NRLMSISE-00 model to obtain relevant atmospheric parameters
    
    According to estimates, F107 changes from 60 to 240, resulting in an error of approximately 0.8% for atmospheric density near 30km,
    The AP changes from 1 to 100, resulting in an error of approximately 1% in the atmospheric density near 30km,
    UT varies from 0 to 24, and the error in atmospheric density near 30km is also within 1%,
    Therefore, for both F107 and AP, the default values are 150 and 4
    We take UT at 20:00 local time, which means UT=12
    """    
    
    output = nrlmsise_output()
    Input = nrlmsise_input()
    flags = nrlmsise_flags()
    aph = ap_array()

    for i in range(7):
        aph.a[i]=4
    flags.switches[0] = 0
    for i in range(1, 24):
        flags.switches[i]=1
        
    
    Input.doy=doy;
    Input.year=0; #/* without effect */
    #/* seconds in day (UT) */
    Input.sec=12*3600;
    # alt: unit km
    Input.alt=alt;
    Input.g_lat=19.5;
    Input.g_long=109.1;
    Input.lst=Input.sec/3600 + Input.g_long/15
    Input.f107A=150;
    Input.f107=150;
    Input.ap=4;
    
    gtd7(Input, flags, output)
    # output.t[1] tem: K; 
    # output.d[5]: total mass density(gm/cm3)
    Den_num=output.d[0]+output.d[1]+output.d[2]+output.d[3]+output.d[4]+output.d[6]+output.d[7]
    return output.t[1], output.d[5], Den_num

def cal_atom_press_temp(alt,doy):
    # This function obtains atmospheric related parameters
    tem,den_mass, Den_num=cal_gtd7(alt,doy)
    # p=nRT/V
    # 28.959 Atmospheric molar mass, g/mol
    # 8.3145 universal gas constant,  J/mol/K
    # press J/cm3 => 10^6Pa [1mPa] => hPa
    # Den_num atmospheric number density, cm^-3
    press=den_mass/28.959*8.314*tem* 10**4
    return tem,press,Den_num*10**6

def caltem_gtd7(doy, sec, alt, f107a, f107, ap):
    """
    This function is based on the NRLMSISE-00 model to obtain relevant atmospheric temperature
    """    
    
    output = nrlmsise_output()
    Input = nrlmsise_input()
    flags = nrlmsise_flags()
    aph = ap_array()

    for i in range(7):
        aph.a[i]=4
    flags.switches[0] = 0
    for i in range(1, 24):
        flags.switches[i]=1
        
    
    Input.doy=doy;
    Input.year=0; #/* without effect */
    #/* seconds in day (UT) */
    Input.sec=sec;
    # alt: unit km
    Input.alt=alt;
    Input.g_lat=19.5;
    Input.g_long=109.1;
    Input.lst=Input.sec/3600 + Input.g_long/15
    Input.f107A=f107a;
    Input.f107=f107;
    Input.ap=ap;
    
    gtd7(Input, flags, output)
    # output.t[1] tem: K; 
    # output.d[5]: total mass density(gm/cm3)
    Den_num=output.d[0]+output.d[1]+output.d[2]+output.d[3]+output.d[4]+output.d[6]+output.d[7]
    return output.t[1], output.d[5], Den_num


if __name__ == '__main__':
    #start = time.clock()
    tem,press=cal_atom_press_temp(alt=50,doy=172)
    #print(time.clock() - start)













