
import numpy as np
import os
import datetime 
import pickle
from astropy.time import Time
import pretem

def ReadTxtUTC(filename_abs, getMore = 0):
    # This program reads data from TXT data files
    with open(filename_abs,'r') as f:
        lines = f.readlines()
        distanceList = []
        Count = []
        for line in lines[20:]:
            datalist=line.strip().split(' ')
            datalist =[i for i in datalist if i !='']
            distanceList.append(float(datalist[0]))
            Count.append(int(datalist[1]))
        distanceList = np.array(distanceList)
        Count = np.array(Count)

        start_time_utc = datetime.datetime.strptime(lines[8][16:35],"%Y-%m-%dT%H:%M:%S")
        end_time_utc = datetime.datetime.strptime(lines[9][14:33],"%Y-%m-%dT%H:%M:%S")
        
        
    if getMore == 1:
        return [distanceList,Count,start_time_utc,end_time_utc]
    else:
        return [distanceList,Count]

def createfolder(path_create):
    # The program creates a new folder if it does not exist
    if not os.path.isdir(path_create):
        os.makedirs(path_create)
        print('created a floder'+path_create)
    return None

def get_spec_file(root_dir,suffix):
    # This program obtains all files of the specified type in the root directory
    specfile=[]
    if not os.path.exists(root_dir):
        return []
    
    for filename in os.listdir(root_dir):   
        if os.path.splitext(filename)[1]=='.'+suffix:
            specfile.append(filename)
    return(specfile)  


def ISOtoSZA(UTC, Loc=[19.5,109.1,70]):
    # This program converts UTC time to solar zenith angle (SZA)
    # UTC, format=datetime.datetime, scale=utc
    # Loc, [0]:glot[deg], [1]:glong[deg], [2]:alt[m]
    # return, SZA[deg]

    from astropy.time import Time   
    from astropy import units as u
    from astropy.coordinates import EarthLocation, AltAz
    from astropy.coordinates import get_sun
    
    loc = EarthLocation(lat=Loc[0]*u.deg, lon=Loc[1]*u.deg, height=Loc[2]*u.m)
    timestamp = Time(UTC.strftime("%Y-%m-%dT%H:%M:%S"), format='isot', scale='utc')
    frame = AltAz(obstime = timestamp, location = loc)
    sunaltaz = get_sun(timestamp).transform_to(frame)
    
    return 90-sunaltaz.alt.value

def Save_all_data(dir_data, date, dir_cache):
    """
    This program integrates a day's data into an array and 
    stores the array in a Python specific format using pickle
    
    dir_data:  The absolute address for storing all data with date as the file name
    date:      yyyymmdd
    dir_cache: Pickle data storage address
    return:    None       
    data_raw:  [nfile, nhigh]
    time_raw:  [nfile], s, unix
    SZA_raw：  [nfile], deg
    alt：      [nhigh], km
    """

    dir_date = os.path.join(dir_data, date)
    dir_raw = os.path.join(dir_date, 'raw')
    
    if not os.path.exists(dir_raw):
        return None

    num_high=3000
    data_raw=np.zeros((0, num_high))
    time_raw=np.array([])
    SZA_raw=np.array([])

    for j,filename_read in enumerate(get_spec_file(dir_raw, 'TXT')):          

        alt, counts,start_time_utc,end_time_utc = ReadTxtUTC(
            os.path.join(dir_raw,filename_read),getMore=1)
        
        """ Ensure a basic altitude resolution of 1km"""
        altreso_base=alt[2]-alt[1]
        n_stand=np.ceil(len(alt)/(1/altreso_base))
        alt_stand=np.arange(n_stand)
        counts_stand=np.arange(n_stand)
        for i in range(int(n_stand)):
            if len(alt)>i*(1/altreso_base):
                num=1/altreso_base
            else:
                num=len(alt)-(i-1)*(1/altreso_base)
            alt_stand[i]=sum(alt[i*round(1/altreso_base):round(i*(1/altreso_base)+num)])/num
            counts_stand[i]=sum(counts[i*round(1/altreso_base):round(i*(1/altreso_base)+num)])
        counts = counts_stand

        data_raw = np.vstack(( data_raw, 
                            np.hstack((  counts.reshape(1, len(counts)), np.zeros((1,num_high-len(counts)))
                                        )) 
                            ))
        time_data=start_time_utc+(end_time_utc-start_time_utc)/2
        timestamp = Time(time_data.strftime("%Y-%m-%dT%H:%M:%S.%f"), format='isot', scale='utc')
        time_raw = np.hstack((time_raw, timestamp.unix))
        SZA_raw = np.hstack((SZA_raw, ISOtoSZA(time_data) ))

    alt=np.arange(3000)+0.5

    if not os.path.exists(dir_cache):
        os.makedirs(dir_cache)
    absfile_save=os.path.join(dir_cache, date+'.pickle')
    with open( absfile_save,'wb') as f:
        pickle.dump([data_raw, time_raw, SZA_raw, alt], f)
        
    return None

def fac_saturation(times, alt):
    # This function calculates the factor of saturation effect
    # alt: unit,km

    alt=alt*10**3           # unit, m
    h=6.62607015*10**(-34);       # Js, planck constant
    c0=2.99792458*10**(8);        # m/s, Lightspeed

    para_laser, para_atm, para_tele=para_lidar()

    Aki=1.0216*10**7   # s-1
    tau_R=1./Aki    # unit, s 
    Delta_tL=10.*1e-9    # unit, s
    N_L=para_laser[3]/3*1e-3*para_laser[0]/h/c0    #   unit, counts
    T=0.85     # unit, %
    sigma_eff=mehe_scatter_cross(times, alt)   # unit, m^2
    theta_L=para_laser[2]     
    Omega=np.pi*theta_L**2/4;
    
    t_S=alt**2.*Omega*Delta_tL/(2*sigma_eff*N_L*T)
    Sat=1./(1.+tau_R/t_S)*(1.-tau_R/Delta_tL*(tau_R/t_S)/(1.+(tau_R/t_S))*
        (np.exp(-Delta_tL/tau_R*(1.+tau_R/t_S))-1.)) 
    return Sat

def rolling_window(a, window):
    # This function calculates the sliding window based on the given array
    # a: One-dimensional array
    # window: scalar data 
    pad = np.ones(len(a.shape), dtype=np.int32)
    pad[-1] = window-1
    pad = list(zip(pad - int(pad/2), np.zeros(len(a.shape),
                             dtype=np.int32) + int(pad/2)))
    
    a = np.pad(a, pad,mode='reflect')
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)

def filter_data(filename_abs, location = 'OP'):
    # This function acquires the solar position for each time interval, then performs several
    # filters. It returns the filtered data, corresponding timestamps and solar zenith angles
    # and the index of local solar midnight.
    
    with open(filename_abs, 'rb') as f:
            data_all = pickle.load(f)

    data = data_all[0]
    time = data_all[1]
    SZA = data_all[2]
    alt = data_all[3]


    bg = np.mean(data[:,1500:2500], axis=1)
        
    N = max(np.where(SZA > 96)[0][0], 1)
    M = min(np.where(SZA > 96)[0][-1], len(bg) - 2)
    ds = bg[N:M]
    times_section = time[N:M]
    data_section = data[N:M,:]
    SZA_section = SZA[N:M]

    if N==M:
        return np.zeros((0,3000)),np.zeros(0),np.zeros(0),180
    
    ds_median = np.median(rolling_window(ds, 60), axis=-1)
    mask = np.where(ds < 2*ds_median)[0]
    
    if np.sum(mask)==0:
        return np.zeros((0,3000)),np.zeros(0),np.zeros(0),180


    he = np.mean(data_section[mask, 300:1500], axis=1) - ds[mask]
    he_median = np.median(rolling_window(he, 60), axis=-1)
    std = np.std(rolling_window(he, 60), axis=-1)
    mask_temp = np.where(np.abs(he - he_median) < 3*std)[0]
    mask = mask[mask_temp]
    
    if np.sum(mask)==0:
        return np.zeros((0,3000)),np.zeros(0),np.zeros(0),180

    metric = (np.mean(data_section[mask,46:60], axis=1)
              - np.mean(data_section[mask,1500:2500], axis=1))
    threshhold = 0.4*np.max(metric)
    mask_temp = np.where(metric > threshhold)[0]
    mask = mask[mask_temp]

    if np.sum(mask)==0:
        return np.zeros((0,3000)),np.zeros(0),np.zeros(0),180
    
    split = np.where(SZA_section[mask] == 
                     np.max(SZA_section[mask]))[0][0]
    
    return data_section[mask,:], times_section[mask], SZA_section[mask], split

def para_lidar():
    # This function defines unified system parameters for use by other functions
    # some laser parameter
    lambda_laser=1083.034*10**(-9)      # atmosphere lambda, unit: m
    rmswid_laser=62*10**6               # uit: Hz
    div_laser = 1*1e-3/30     # Laser divergence angle after beam expansion, unit rad
    energy_laser = 200      # Single pulse Laser energy, unit 200mJ
    frequ_laser = 50       # Hz, Laser frequency
    para_laser=[lambda_laser, rmswid_laser, div_laser, energy_laser, frequ_laser]
    # some atmosphere parameter
    tem_atm450=1500    # K
    para_atm=tem_atm450
    # some telescope parameter
    NA_fiber = 0.27      
    diam_fiber = 62.5*1e-6    # unit, m
    diam_tele = 1000*1e-3     # unit, m
    para_tele = [diam_tele, diam_fiber, NA_fiber]
    return para_laser, para_atm, para_tele

def tem_atm(times_unix, alt):
    # This function obtains the atmospheric temperature based on time, 
    # altitude and the atmospheric model.
    alt = np.array([alt])
    n_alt = np.size(alt)
    alt=alt.reshape(n_alt)
    tem_atm = np.zeros_like(alt, dtype=float)

    times = datetime.datetime.fromtimestamp(times_unix)
    doy = times.timetuple().tm_yday  
    seconds_since_midnight = (times - times.replace(hour=0, minute=0, second=0)
                             ).total_seconds()  
    
    dir_current = os.path.abspath(".")
    fname = os.path.join(dir_current,"Geomagnetic and solar indices.txt")
    data_ap = np.loadtxt(fname, skiprows= 40 )
    t = seconds_since_midnight/(3600*3)
    ap = data_ap[:,15+int(t)]
    mm = times.month
    dd = times.day
    year = times.year
    mask_date = (data_ap[:,0]==year) & (data_ap[:,1]==mm) & (data_ap[:,2]==dd)
    # index = np.where(mask_date == True)[0][0]
    index=100
    ap = data_ap[mask_date,15+int(t)]
    f107 = data_ap[mask_date, 25]
    f107a = np.mean(data_ap[(index-40):(index+41), 25])
    
    for i in range(n_alt):
        tem_atm[i],pre_atm,Den_num=pretem.caltem_gtd7(doy, seconds_since_midnight, 
                                                      alt[i], f107a, f107, ap)
    return tem_atm


def mehe_scatter_cross(times, alt):
    # This function calculates the effective scattering cross-section of metastable helium

    # some const.
    c0=2.99792458*10**(8);      #m/s, Lightspeed
    kb=1.380649*10**(-23);      #J/K, Boltzmann constant
    mp=1.672621637*10**(-27);   #kg, Proton mass
    mn=1.674927211*10**(-27);   #kg, neutron mass
    me=9.10938215*10**(-31);    #kg,  Electronic quality

    # some parameter
    para_laser, para_atm, para_tele=para_lidar()
    
    # some atmosphere parameter
    T= tem_atm(times, alt)     # K
    # some atom parameter
    A=1.0216*10**7          # Einstein's A-coefficient, unit s-1
    
    # atmosphere lambda of he(2^3S)
    lambda_1=1082.909114*10**(-9)    # unit: m
    lambda_2=1083.025010*10**(-9)    # unit: m
    lambda_3=1083.033977*10**(-9)    # unit: m
    # frequen of return laser
    nu_1=c0/lambda_1     # unit: Hz
    nu_2=c0/lambda_2     # unit: Hz
    nu_3=c0/lambda_3     # unit: Hz
    # some laser parameter
    lambda_laser=para_laser[0]
    rmswid_laser=para_laser[1]
    nu_L=c0/lambda_laser
    sigma_L=rmswid_laser
    
    
    sigma_D1=nu_1*np.sqrt(kb*T/(2*mp+2*mn+2*me)/c0**2)
    sigma_D2=nu_2*np.sqrt(kb*T/(2*mp+2*mn+2*me)/c0**2)
    sigma_D3=nu_3*np.sqrt(kb*T/(2*mp+2*mn+2*me)/c0**2)
    
    sigma_1=1/np.sqrt(2*np.pi)/sigma_D1*(c0/nu_1)**2/8/np.pi*3*A
    sigma_2=1/np.sqrt(2*np.pi)/sigma_D2*(c0/nu_2)**2/8/np.pi*3*A
    sigma_3=1/np.sqrt(2*np.pi)/sigma_D3*(c0/nu_3)**2/8/np.pi*3*A
    
    sigma_eff1=sigma_D1*sigma_1/np.sqrt(sigma_D1**2+sigma_L**2)*\
        np.exp(-(nu_1-nu_L)**2/2/(sigma_D1**2+sigma_L**2))
    sigma_eff2=sigma_D2*sigma_2/np.sqrt(sigma_D2**2+sigma_L**2)*\
        np.exp(-(nu_2-nu_L)**2/2/(sigma_D2**2+sigma_L**2))
    sigma_eff3=sigma_D3*sigma_3/np.sqrt(sigma_D3**2+sigma_L**2)*\
        np.exp(-(nu_3-nu_L)**2/2/(sigma_D3**2+sigma_L**2))
    sigma_eff=sigma_eff1/11+sigma_eff2*10/33+sigma_eff3*20/33
    return sigma_eff

def get_He_profiles(data, times, windows_rolling_alt=0, windows_rolling_time=0):
    # This function takes filtered data and calculates raw photon return
    # profiles and metastable helium density profiles, each with associated
    # uncertainties, over 50 km vertical bins.
    # data    : [n_time, n_alt]
    # times   : [n_time]
    
    # read geometric factor data
    dir_current = os.path.abspath(".")
    fname = os.path.join(dir_current,"geo_data.txt")
    geo_data = np.loadtxt(fname, skiprows=1)

    time = np.mean(times)

    if data.shape[0]==0:
        return np.zeros(30),np.zeros(30),np.zeros(30),np.zeros(30)

    bg = np.sum(data[:,1500:2500], axis=1)
    rayleigh = np.sum(np.mean(data[:,46:60], axis=1) - np.mean(data[:,1500:2500],
                      axis=1))
    
    bin_size = 50 
    bg_bin = 1000 
    
    n = int(bin_size / 1.)   
    N = int(1500  / n)  
    binned_data = np.zeros((data.shape[0], N), dtype=np.float64)
    z_squared_sighe = np.zeros(N)
    
    binned_bg_sub = np.zeros((data.shape[0], N), dtype=np.float64)
    binned_bg_sub_error = np.zeros((data.shape[0], N), dtype=np.float64)

    for i in range(N):
        binned_data[:,i] = np.sum(data[:,n*i:n*(i + 1)], axis=1)

        z_squared_sighe[i] = np.mean( (np.arange(0,1500,1)[n*i:n*(i + 1)]+0.5)**(-2)*
                               geo_data[n*i:n*(i + 1),1]*
                               fac_saturation(time, np.arange(0,1500,1)[n*i:n*(i + 1)]+0.5) 
                               * mehe_scatter_cross(time, np.arange(0,1500,1)[n*i:n*(i + 1)]+0.5)
                               )  # km^2
    
    for i in range(data.shape[0]):
        if windows_rolling_alt:
            binned_data[i,:]=np.mean(rolling_window(binned_data[i,:].reshape(len(binned_data[i,:])), windows_rolling_alt), axis=-1)
        binned_bg_sub[i,:] = binned_data[i,:] - n/bg_bin*bg[i]
        binned_bg_sub_error[i,:] = np.sqrt(binned_data[i,:] + (n/bg_bin)**2*bg[i])
  
    integ_bin_bg_sub = np.sum(binned_bg_sub, axis = 0)
    integ_bin_bg_sub_error = np.sqrt(np.sum(binned_bg_sub_error**2,axis = 0))



    para_laser, para_atm, para_tele=para_lidar()
    lam_laser=para_laser[0]     # unit: m
    sig_back_R=3.462*10**(-33)    # unit: m^2 sr^-1
    nR_z_R2=np.zeros(14)
    for i in np.arange(46,60,1):
        # some atmosphere parameter
        # pre_atm(alt_ray)  Atmospheric pressure,   unit mPa
        # tem_atm(alt_ray)  atmospheric temperature,   unit K
        # Den_num   atmospheric number density   unit m^-3
        tem_atm,pre_atm,Den_num=pretem.cal_atom_press_temp(alt=i+0.5,doy=6)
        
        nR_z_R2[i-46]=Den_num*geo_data[i,1]/((i+0.5)*1000.)**2      # m^-1
        

    n_R_over_z_R_squared=np.mean(nR_z_R2)


    delta_z_R = 1000               # m
    conversion_factor = (rayleigh  * bin_size / n_R_over_z_R_squared 
                         * z_squared_sighe / (sig_back_R*4*np.pi) / delta_z_R * 10**3) # counts cm^3
    

    density = integ_bin_bg_sub / conversion_factor # cm^-3
    density_error = np.sqrt((integ_bin_bg_sub_error / conversion_factor)**2
                            + (0.3*density)**2) # cm^-3
    
    return integ_bin_bg_sub, integ_bin_bg_sub_error, density, density_error

def data_cache(dir_cache, dates, dir_data):
    # This program preprocesses data and 
    # places the preprocessed data in a specific folder
    for filename in os.listdir(dir_data):
        if os.path.isdir(os.path.join(dir_data, filename)):
            if len(filename)==8 and int(filename)>20231111:
                if filename in dates:
                    fname_pickle = os.path.join(dir_cache, filename+'.pickle')
                    if not os.path.exists(fname_pickle):
                        Save_all_data(dir_data, filename, dir_cache)

if  __name__=='__main__':
    para_laser, para_atm, para_tele=para_lidar()
    rad_tele = para_tele[1]*para_tele[2]/para_tele[0]    # unit, rad
    print(rad_tele)

