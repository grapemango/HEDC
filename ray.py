"""检查滤波效果"""
import inversion as al
import numpy as np
import os
import datetime 
from astropy.time import Time
from matplotlib import pyplot as plt



def plot_ray(data, times, accutype, highview, sza_edges, twilight_type, figname_prefix, dir_save):
    # This program draws Line graphics
    # data: List, with each internal item containing data within each SZA range
    # data[i]: i=0-3, raw data array, N * 3000, where N represents the number of raw data and 3000 represents the number of height units
    # times: List, with each internal item containing times(unit) data within each SZA range
    # times[i]: Original data array, (N,), N represents the number of original data
    # accutype: Data accumulation type, 0 represents 1 day of data, and 1 represents the accumulation of all. pickle files in dir_cache
    # highview: Display height range when drawing, integer, must be divisible by 50, 250-1450
    # sza_edges: sza coordinate data corresponding to data, 2 * N array, first row min, second row max
    # figtype: Image type, 0 represents the inversion number density, and 1 represents the original number of photons with background removed
    # twilight_type: Morning and dusk type, 0 represents Dusk, 1 represents Dawn
    # figname_prefix: The part of the image name that represents the date
    # dirsave: Path to save folder
    # return: None
    
    
    nametype_twi = ['dusk','dawn']
    nametype_file = '_'+nametype_twi[twilight_type]+'_counts.png'

    
    i=0
    ray = np.sum(data[i][0:20,:], axis=0)
    bg = np.mean(ray[1500:2500], axis=-1)
    time = np.mean(times[0][0:20])
    SZA = al.ISOtoSZA(datetime.datetime.strptime(Time(time, format='unix').iso,'%Y-%m-%d %H:%M:%S.%f'), Loc=[19.5,109.1,70])
    SZA_round = round(SZA, 2)

    plt.rcParams['font.size'] = 14
    plt.figure()
    plt.plot(ray,np.arange(3000)+0.5,
                label='SZA:'+ "{:.2f}".format(SZA_round) )
    plt.fill_betweenx(np.arange(3000)+0.5,
                    ray-np.sqrt(ray),
                    ray+np.sqrt(ray),
                    alpha=0.2)
    plt.plot([bg,bg],[min(np.arange(3000)+0.5),max(np.arange(3000)+0.5)], color='black')
    plt.ylabel('Altitude [km]')
    plt.xlabel('Counts')
    plt.ylim([30,highview])
    plt.xscale('log')
    if accutype == 1:
        plt.xlim([10**1,10**4])
    elif accutype == 0:
        plt.xlim([10**1,10**4])

    plt.legend()
    plt.title(figname_prefix+' (10min, 1km)')
    plt.savefig(os.path.join(dir_save,figname_prefix+nametype_file),dpi=300)    
    # plt.show()
    plt.close()
    return None



if  __name__=='__main__':

    accutype = 0
    highview_even=800

    
    dir_current = os.path.dirname(os.path.abspath(__file__))
    dir_data = os.path.split(dir_current)[0]

    dir_cache = os.path.join(dir_data, 'data_rayl')
    al.createfolder(dir_cache)
    

    date='20240106'
    al.data_cache(dir_cache, date, dir_data)

    figname=date

    dir_save = os.path.join(dir_cache, figname)
    al.createfolder(dir_save)


    #Filters the data, gets solar elevation data and separates into morning/evening
    data_evenings = []
    data_mornings = []
    times_evenings = []
    times_mornings = []
    SZA_evenings = []
    SZA_mornings = []

    nday = 0
    for i,fname in enumerate(al.get_spec_file(dir_cache, 'pickle')):
        if (accutype == 0) and (fname == date+'.pickle'):
            data, time, SZA, split = al.filter_data(os.path.join(dir_cache, fname))
            data_evenings.append(data[:split,:])
            data_mornings.append(data[split:,:])
            times_evenings.append(time[:split])
            times_mornings.append(time[split:])
            SZA_evenings.append(SZA[:split])
            SZA_mornings.append(SZA[split:])
            nday += 1
        elif accutype == 1:
            data, time, SZA, split = al.filter_data(os.path.join(dir_cache, fname))
            data_evenings.append(data[:split,:])
            data_mornings.append(data[split:,:])
            times_evenings.append(time[:split])
            times_mornings.append(time[split:])
            SZA_evenings.append(SZA[:split])
            SZA_mornings.append(SZA[split:])
            nday += 1


    


    #Sorts the data based on solar elevation angles.
    sza_edges_min = np.arange(96,140,10)
    sza_edges_max = np.arange(106,150,10)


    data_szasplit_evening = []
    data_szasplit_morning = []
    for i in range(len(sza_edges_min)):
        data_szasplit_evening.append(np.zeros((0,3000)))
        data_szasplit_morning.append(np.zeros((0,3000)))
    
    for i in range(nday):
        for j in range(len(sza_edges_min)):
            data_szasplit_evening[j] = np.vstack((data_szasplit_evening[j],
                            data_evenings[i][(SZA_evenings[i] < sza_edges_max[j])
                            & (SZA_evenings[i] > sza_edges_min[j])]))
            data_szasplit_morning[j] = np.vstack((data_szasplit_morning[j],
                            data_mornings[i][(SZA_mornings[i] < sza_edges_max[j])
                            & (SZA_mornings[i] > sza_edges_min[j])]))
            

    times_szasplit_evening = []
    times_szasplit_morning = []
    for i in range(len(sza_edges_min)):
        times_szasplit_evening.append(np.zeros((0,)))
        times_szasplit_morning.append(np.zeros((0,)))
    
    for i in range(nday):
        for j in range(len(sza_edges_min)):
            times_szasplit_evening[j] = np.hstack((times_szasplit_evening[j],
                            times_evenings[i][(SZA_evenings[i] < sza_edges_max[j])
                            & (SZA_evenings[i] > sza_edges_min[j])]  ))
            times_szasplit_morning[j] = np.hstack((times_szasplit_morning[j],
                            times_mornings[i][(SZA_mornings[i] < sza_edges_max[j])
                            & (SZA_mornings[i] > sza_edges_min[j])]  ))
            
    
    
    plot_ray(data_szasplit_morning, times_szasplit_morning,
              accutype, 100, np.vstack((sza_edges_min, sza_edges_max)), 1, date, dir_save)




