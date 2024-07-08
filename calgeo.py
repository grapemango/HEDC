import numpy as np
import os
import datetime 
import pickle
from astropy.time import Time
import pretem
from scipy.integrate import quad
import inversion as al



def geofactor_single(R1, R2, d):
    # This function computes the geometric factor for a single height, 
    # and the input can be arrays of the same size.
    # sigma = R1/np.sqrt(2)
    sigma = R1/2
    def f(x):
        r21 = d*np.cos(x)-np.sqrt( R2**2-d**2*np.sin(x)**2 )
        r22 = d*np.cos(x)+np.sqrt( R2**2-d**2*np.sin(x)**2 )
        y = np.exp( -r22**2/(2*sigma**2) )-np.exp( -r21**2/(2*sigma**2) )
        return y

    def g(x):
        r22 = d*np.cos(x)+np.sqrt( R2**2-d**2*np.sin(x)**2 )
        y = np.exp( -r22**2/(2*sigma**2)  )-1
        return y
    
    if d>=R2:
        beta = np.arcsin(R2/d)
        result, error = quad(f, -beta, beta)
    elif d<R2:
        beta = np.pi
        result, error = quad(g, -beta, beta)
    
    # print('error:', error)
    # print(R1)

    g = -result/2/np.pi
    return g

def geometric(alt, ftype = 0):
    # This function calculates the geometric overlap factor for different heights
    # alt, km
    # return, geometric factors
    alt = np.array([alt])
    n_alt = np.size(alt)
    alt = alt.reshape(n_alt)
    G = np.zeros_like(alt, dtype=float)
    for k in range(n_alt):
        para_laser, para_atm, para_tele=al.para_lidar()
        rad_laser=para_laser[2]     # unit, rad
        alt_mm = alt[k]*1e6   # unit, mm

        if ftype == 1:
            if k>=1:
                print('R1:', R_laser)
    
       
        O_laser = np.array([[np.cos(np.pi*0/3),np.sin(np.pi*0/3)],
                            [np.cos(np.pi*2/3),np.sin(np.pi*2/3)],
                            [np.cos(np.pi*4/3),np.sin(np.pi*4/3)]])*250   # unit, mm
        
        R_laser = 150 + alt_mm * rad_laser/2   # unit, mm
        
        rad_tele = para_tele[1]*para_tele[2]/para_tele[0]*2    # unit, rad
        
        O_tele = np.array([np.cos(np.pi*0/3),np.sin(np.pi*0/3)])*1750  # unit, mm
        
        R_tele = 500 + alt_mm * rad_tele/2       # unit, mm
    
        if ftype == 1:
            print('R1:', R_laser)
            print('R2:', R_tele)
            print('d1:', np.sqrt(np.sum((O_laser[0,:]-O_tele)**2)) )
            print('d2:', np.sqrt(np.sum((O_laser[1,:]-O_tele)**2)) )
            print('d3:', np.sqrt(np.sum((O_laser[2,:]-O_tele)**2)) )
            print('G1:', geofactor_single(R_laser, R_tele, np.sqrt(np.sum((O_laser[0,:]-O_tele)**2)) ))
            print('G2:', geofactor_single(R_laser, R_tele, np.sqrt(np.sum((O_laser[1,:]-O_tele)**2)) ))
            print('G3:', geofactor_single(R_laser, R_tele, np.sqrt(np.sum((O_laser[2,:]-O_tele)**2)) ))
    
        O_tele1=np.zeros_like(O_tele)
        theta = np.linspace(0,np.pi*2/3,50)
        ratio = np.zeros( len(theta) )
        for j in range(len(theta)):
            x = O_tele[0] * np.cos(theta[j]) - O_tele[1] * np.sin(theta[j])
            y = O_tele[0] * np.sin(theta[j]) + O_tele[1] * np.cos(theta[j])
            O_tele1[0] = x
            O_tele1[1] = y
            
            g = np.zeros( 3 )
            for i in range(3):
                d = np.sqrt(np.sum((O_laser[i,:]-O_tele1)**2))
                g[i] = g[i] + geofactor_single(R_laser, R_tele, d)
            ratio[j] = np.mean(g)
        G[k] = np.mean(ratio)
    return G

if  __name__=='__main__':
    alt = np.arange(0,1510,1)+0.5
    G = geometric(alt)
    print(G)
    geo_data = np.hstack( 
        (
            alt.reshape(alt.shape[0],1),
            G.reshape(G.shape[0],1)
         )
          )
    # define saveed filename
    dir_current = os.path.abspath(".")
    fname = os.path.join(dir_current,"geo_data.txt")
    
    head='      alt[km]    geometric factor'
    np.savetxt(fname,geo_data,fmt='%13.6E %20.6E',delimiter=',',newline='\n',header=head,footer='',comments='')
    

    # # read geometric factor data
    # dir_current = os.path.abspath(".")
    # fname = os.path.join(dir_current,"geo_data.txt")
    # geo_data = np.loadtxt(fname, skiprows=1)
    # print(geo_data[(50*3-30):(50*(3 + 1)-30),1].shape)




