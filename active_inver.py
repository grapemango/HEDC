import inversion as al
import numpy as np
import os
from matplotlib import pyplot as plt
from matplotlib import cm
import datetime


def plot_2D_2(data, accutype, highview, sza_edges, figtype, twilight_type, figname_prefix, dir_save, fig, ax, index_ax):
    # This program draws contour maps
    # data: List, with each internal item containing data within each SZA range
    # data[i]: i=0-3, List, representing the four return values of the function get_He_profile internally
    # accutype: Data accumulation type, 0 represents 1 day of data, and 1 represents the accumulation of all. pickle files in dir_cache
    # highview: Display height range when drawing, integer, must be divisible by 50, 250-1450
    # sza_edges: sza coordinate data corresponding to data, 2 * N array, first row min, second row max
    # figtype: Image type, 0 represents the inversion number density, and 1 represents the original number of photons with background removed
    # twilight_type: Morning and dusk type, 0 represents Dusk, 1 represents Dawn
    # figname_prefix: The part of the image name that represents the date
    # dirsave: Path to save folder
    # return: None
    # ax: Axes of image
    local_figtype = [2, 0]
    nametype_co=['Number Density of MetaHelium [cm-3]',
          'Counts after removing background']
    nametype_twi = ['dusk','dawn']
    nametype_file = ['_2D_'+nametype_twi[twilight_type]+'_den.png',
                '_2D_'+nametype_twi[twilight_type]+'_counts.png']

    num_high = int(highview/50-2)

    den=np.zeros((data[0][0].shape[0],0))
    for i in range(len(sza_edges[0,:])):
        den=np.hstack((den,data[i][local_figtype[figtype]].reshape(30,1)))
    
    x_sza=(sza_edges[0,:]+sza_edges[1,:])/2.
    X, Y = np.meshgrid(x_sza, (np.arange(num_high)+3)*50.+25.)
    den[den == 0] = np.nan

    ax2=ax[index_ax]

    Z = den[3:(3+num_high),:]
    if figtype == 0:
        if accutype == 0:
            Z[Z>0.6]=0.6
            Z[Z<-0.6]=-0.6
            
        
    cset1 = ax2.contourf(X, Y, Z, levels=10,cmap=cm.jet)
   
    ax3 = ax2.twiny()
    ax3.set_xlim(ax2.get_xlim()[0], ax2.get_xlim()[1])
    if twilight_type == 1:
        xticks = [al.ISOtoSZA( datetime.datetime(2024,1,6,4)-datetime.timedelta(hours = 109.1/15)),
                        al.ISOtoSZA( datetime.datetime(2024,1,6,5)-datetime.timedelta(hours = 109.1/15)),
                        ]
        print(xticks, ax3.get_xlim())
        time_str = ['5:00', '4:00']
    elif twilight_type == 0:
        if accutype == 0:
            xticks = [
                            al.ISOtoSZA( datetime.datetime(2024,4,2,20)-datetime.timedelta(hours = 109.1/15)),
                            al.ISOtoSZA( datetime.datetime(2024,4,2,21)-datetime.timedelta(hours = 109.1/15)),
                            al.ISOtoSZA( datetime.datetime(2024,4,2,22)-datetime.timedelta(hours = 109.1/15))
                            ]
            time_str = ['20:00', '21:00', '22:00']
        elif accutype == 1:
            xticks = [
                            al.ISOtoSZA( datetime.datetime(2024,2,20,20)-datetime.timedelta(hours = 109.1/15)),
                            al.ISOtoSZA( datetime.datetime(2024,2,20,21)-datetime.timedelta(hours = 109.1/15)),
                            al.ISOtoSZA( datetime.datetime(2024,2,20,22)-datetime.timedelta(hours = 109.1/15))
                            ]
            time_str = ['20:00', '21:00', '22:00']
    mask = (xticks>ax2.get_xlim()[0])&(xticks<ax2.get_xlim()[1])
    xticks_err = np.where(mask, xticks, np.nan)
    
    time_str_lim = np.where(mask, time_str, np.nan)
    ax3.set_xticks(xticks_err)
    ax3.set_xticklabels(time_str_lim)
    ax3.tick_params(axis='x', pad=-4)
    ax3.set_xlabel('LT:', loc = 'left', labelpad=-7)


    cbar = fig.colorbar(cset1, cax=ax[index_ax+2], location='bottom')
    cbar.set_label(nametype_co[figtype], size=9)

    if figtype == 0:
        if twilight_type == 0:
            ax2.set_title('(b) '+figname_prefix, fontsize=10, pad=10)  
        elif twilight_type == 1:
            ax2.set_title('(d) '+figname_prefix, fontsize=10, pad=10)  
    elif figtype == 1:
        if twilight_type == 0:
            ax2.set_title('(a) '+figname_prefix, fontsize=10, pad=10)  
        elif twilight_type == 1:
            ax2.set_title('(c) '+figname_prefix, fontsize=10, pad=10) 

    if twilight_type:
        ax2.invert_xaxis()

   
    if figtype == 0:
        if twilight_type == 1:
            if accutype == 0:
                ticks = np.arange(-0.3, 0.6, 0.2) 
                ax[index_ax+2].set_xticks(ticks)
            elif accutype == 1:
                ticks = np.arange(-0.1, 0.3, 0.1) 
                ax[index_ax+2].set_xticks(ticks)
        elif twilight_type == 0:
            if accutype == 1:
                # ticks = np.arange(-0.3, 0.6, 0.15)  
                ticks = np.array([-0.1,0.1 ,0.3])
                ax[index_ax+2].set_xticks(ticks)
            elif accutype == 0:
                ticks = np.arange(-0.6, 0.7, 0.3)
                ax[index_ax+2].set_xticks(ticks)

    ax2.set_ylabel('Altitude [km]',labelpad=0)
    ax2.set_ylim([200,highview])
    ax2.set_xlabel('SZA [deg]',labelpad=0)
    if figtype ==1:
        if twilight_type == 0:
            fig.savefig(os.path.join(dir_save,figname_prefix+nametype_file[figtype]),dpi=600)
    return None

def plot_line_2(data, accutype, highview, sza_edges, figtype, twilight_type, figname_prefix, dir_save, fig, ax, index_ax):
    # This program draws Line graphics
    # data: List, with each internal item containing data within each SZA range
    # data[i]: i=0-3, List, representing the four return values of the function get_He_profile internally
    # accutype: Data accumulation type, 0 represents 1 day of data, and 1 represents the accumulation of all. pickle files in dir_cache
    # highview: Display height range when drawing, integer, must be divisible by 50, 250-1450
    # sza_edges: sza coordinate data corresponding to data, 2 * N array, first row min, second row max
    # figtype: Image type, 0 represents the inversion number density, and 1 represents the original number of photons with background removed
    # twilight_type: Morning and dusk type, 0 represents Dusk, 1 represents Dawn
    # figname_prefix: The part of the image name that represents the date
    # dirsave: Path to save folder
    # return: None
    local_figtype = [2, 0]
    nametype_co=['Number Density of MetaHelium [cm^-3]',
          'Counts after removing background']
    nametype_twi = ['dusk','dawn']
    nametype_file = ['_'+nametype_twi[twilight_type]+'_den.png',
                '_'+nametype_twi[twilight_type]+'_counts.png']

    ax1=ax[index_ax]

    for i in range(len(sza_edges[0,:])):
        ax1.plot(data[i][local_figtype[figtype]],np.arange(30)*50+25,
                 label=str(sza_edges[0,i])+'-'+str(sza_edges[1,i]))
        ax1.fill_betweenx(np.arange(30)*50+25,
                      data[i][local_figtype[figtype]]-data[i][local_figtype[figtype]+1],
                      data[i][local_figtype[figtype]]+data[i][local_figtype[figtype]+1],
                      alpha=0.2)
    ax1.plot([0,0],[min(np.arange(30)*50+25),max(np.arange(30)*50+25)], color='black')
    ax1.set_ylim([200,highview])
    ax1.set_yticklabels([])
    ax1.yaxis.set_ticks_position('none')

    if figtype == 0:
        if accutype == 1:
            if twilight_type == 0:
                ax1.set_xlim([-0.25,0.5])
                ax1.set_xticks([0.0,0.5])
            elif twilight_type == 1:
                ax1.set_xlim([-0.15,0.3])
                ax1.set_xticks([0.0,0.3])
        elif accutype == 0:
            if twilight_type == 0:
                ax1.set_xlim([-0.4,0.8])
                ax1.set_xticks([0.0,0.8])
            elif twilight_type == 1:
                ax1.set_xlim([-0.3,0.6])
                ax1.set_xticks([0.0,0.6])
    elif figtype == 1:
        if accutype == 1:
            if twilight_type == 0:
                ax1.set_xlim([-500,1000])
            elif twilight_type == 1:
                ax1.set_xlim([-2500,5000])
        elif accutype == 0:
            if twilight_type == 0:
                ax1.set_xlim([-250,500])
            elif twilight_type == 1:
                ax1.set_xlim([-500,1000])

    leg = ax1.legend(bbox_to_anchor=(-2.3, 1.22), 
                     loc='upper left', 
                     frameon=False,   
                     handlelength=0.5, 
                    handletextpad=0.1,  
                    fontsize=8,   
                    columnspacing=0.5, 
                    ncol = 4
                    )
    ax1_xlim = ax1.get_xlim()
    ax1_ylim = ax1.get_ylim()
    ax1.text(-(ax1_xlim[1]-ax1_xlim[0])*3, (ax1_ylim[1]-ax1_ylim[0])*1.33
             , 'SZA [deg]:', fontsize=8, color='black', ha='center', va='bottom')
    
    return None


def main(accutype, date, dir_cache, figname, twilight_type, fig, ax, sza_end_2D):
    # This function filters and inverts the original data and draws a graph

    dir_save = os.path.join(dir_cache, date)
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
        if (accutype == 0) and (fname == date+'.pickle') :
            data, time, SZA, split = al.filter_data(os.path.join(dir_cache, fname))
            data_evenings.append(data[:split,:])
            data_mornings.append(data[split:,:])
            times_evenings.append(time[:split])
            times_mornings.append(time[split:])
            SZA_evenings.append(SZA[:split])
            SZA_mornings.append(SZA[split:])
            nday += 1
        elif (accutype == 1) and (int(os.path.splitext(fname)[0])!=20240101):
            data, time, SZA, split = al.filter_data(os.path.join(dir_cache, fname))
            data_evenings.append(data[:split,:])
            data_mornings.append(data[split:,:])
            times_evenings.append(time[:split])
            times_mornings.append(time[split:])
            SZA_evenings.append(SZA[:split])
            SZA_mornings.append(SZA[split:])
            nday += 1

    sza_edges_min = np.arange(96,130,10)
    sza_edges_max = np.arange(106,140,10)

    # Integrate data according to the predetermined range of solar zenith angles
    data_szasplit_evening = []
    data_szasplit_morning = []
    times_szasplit_evening = []
    times_szasplit_morning = []

    # The format of elements in the reservation list
    for i in range(len(sza_edges_min)):
        data_szasplit_evening.append(np.zeros((0,3000)))
        data_szasplit_morning.append(np.zeros((0,3000)))
        times_szasplit_evening.append(np.zeros((0,)))
        times_szasplit_morning.append(np.zeros((0,)))
    for i in range(nday):
        for j in range(len(sza_edges_min)):
            data_szasplit_evening[j] = np.vstack( (
                            data_szasplit_evening[j],
                            data_evenings[i][ (SZA_evenings[i] < sza_edges_max[j])
                            & (SZA_evenings[i] > sza_edges_min[j]) ]
                            ) )
            data_szasplit_morning[j] = np.vstack( (
                            data_szasplit_morning[j],
                            data_mornings[i][ (SZA_mornings[i] < sza_edges_max[j])
                            & (SZA_mornings[i] > sza_edges_min[j]) ]
                            ) )
            times_szasplit_evening[j] = np.hstack( (
                            times_szasplit_evening[j],
                            times_evenings[i][(SZA_evenings[i] < sza_edges_max[j])
                            & (SZA_evenings[i] > sza_edges_min[j])]
                            ) )
            times_szasplit_morning[j] = np.hstack( (
                            times_szasplit_morning[j],
                            times_mornings[i][(SZA_mornings[i] < sza_edges_max[j])
                            & (SZA_mornings[i] > sza_edges_min[j])]
                            ) )
    
    # Accumulate and invert metastable helium density from data within a certain solar zenith angle range
    profile_szasplit_evening = []
    profile_szasplit_morning = []
    # print(times_szasplit_evening[0].shape)
    for j in range(len(sza_edges_min)):
        profile_szasplit_evening.append(al.get_He_profiles(data_szasplit_evening[j], times_szasplit_evening[j], windows_rolling_alt=3)) 
        profile_szasplit_morning.append(al.get_He_profiles(data_szasplit_morning[j], times_szasplit_morning[j], windows_rolling_alt=3)) 
 
    if twilight_type==0:
        index_ax=[4,1]
        for figtype in range(2):
            plot_line_2(profile_szasplit_evening, accutype, 1000, 
                        np.vstack((sza_edges_min, sza_edges_max)), figtype, 0, figname, dir_save, fig, ax, index_ax[figtype])
    elif twilight_type==1: 
        index_ax=[10,7]
        for figtype in range(2):
            plot_line_2(profile_szasplit_morning, accutype, 1000, 
                    np.vstack((sza_edges_min, sza_edges_max)), figtype, 1, figname, dir_save, fig, ax, index_ax[figtype])

    
    # Processing data for drawing contour maps
    sza_edges_min = np.arange(96,sza_end_2D,1)
    sza_edges_max = np.arange(106,sza_end_2D+10,1)

    # Integrate data according to the predetermined range of solar zenith angles
    data_szasplit_evening = []
    data_szasplit_morning = []
    times_szasplit_evening = []
    times_szasplit_morning = []
    for i in range(len(sza_edges_min)):
        data_szasplit_evening.append(np.zeros((0,3000)))
        data_szasplit_morning.append(np.zeros((0,3000)))
        times_szasplit_evening.append(np.zeros((0,)))
        times_szasplit_morning.append(np.zeros((0,)))
    
    for i in range(nday):
        for j in range(len(sza_edges_min)):
            data_szasplit_evening[j] = np.vstack((data_szasplit_evening[j],
                            data_evenings[i][(SZA_evenings[i] < sza_edges_max[j])
                            & (SZA_evenings[i] > sza_edges_min[j])]))
            data_szasplit_morning[j] = np.vstack((data_szasplit_morning[j],
                            data_mornings[i][(SZA_mornings[i] < sza_edges_max[j])
                            & (SZA_mornings[i] > sza_edges_min[j])]))
            times_szasplit_evening[j] = np.hstack( (
                            times_szasplit_evening[j],
                            times_evenings[i][(SZA_evenings[i] < sza_edges_max[j])
                            & (SZA_evenings[i] > sza_edges_min[j])]
                            ) )
            times_szasplit_morning[j] = np.hstack( (
                            times_szasplit_morning[j],
                            times_mornings[i][(SZA_mornings[i] < sza_edges_max[j])
                            & (SZA_mornings[i] > sza_edges_min[j])]
                            ) )
            
    # Accumulate and invert metastable helium density from data within a certain solar zenith angle range
    profile_szasplit_evening = []
    profile_szasplit_morning = []
    for j in range(len(sza_edges_min)):
        profile_szasplit_evening.append(al.get_He_profiles(data_szasplit_evening[j], times_szasplit_evening[j],windows_rolling_alt=3)) 
        profile_szasplit_morning.append(al.get_He_profiles(data_szasplit_morning[j], times_szasplit_morning[j],windows_rolling_alt=3)) 

    if twilight_type==0:
        index_ax=[3,0]
        for figtype in range(2):
            plot_2D_2(profile_szasplit_evening, accutype, highview_even, 
                    np.vstack((sza_edges_min, sza_edges_max)), figtype, 0, figname, dir_save, fig, ax, index_ax[figtype])
    elif twilight_type==1:
        index_ax=[9,6]
        for figtype in range(2):
            plot_2D_2(profile_szasplit_morning, accutype, 1000, 
                    np.vstack((sza_edges_min, sza_edges_max)), figtype, 1, figname, dir_save, fig, ax, index_ax[figtype])

    return None

def fig_axes():
    # This function establishes the axes of the figure
    ax=[]

    width_fig = 6.5
    heigh_fig = 6
    fig = plt.figure(figsize=(width_fig,heigh_fig))    
    size_font = 10
    heigh_font = size_font/72/heigh_fig

    width_ticks = 0.01
    width_ticks_labels = heigh_font*2.5
    width_axis_label = 0.03
    width_label = width_ticks+width_ticks_labels+width_axis_label

    heigh_ticks = 0.01
    heigh_ticks_labels = 0.035
    heigh_axis_label = 0.04
    heigh_label = heigh_ticks+heigh_ticks_labels+heigh_axis_label

    left_bar = width_label+(0.5-width_label)/8
    bottom_bar = heigh_label
    width_bar = (0.5-width_label)*6/8
    heigh_bar = 0.02

    heigh_title = heigh_font + 10/72/heigh_fig

    heigh_axes = 0.5-heigh_label-heigh_bar-bottom_bar-heigh_font-heigh_title
    bottom_axes = heigh_label+heigh_bar+bottom_bar
    width_axes = 0.5-width_label-0.02
    width_2D_axes = width_axes*3/4
    width_line_axes = width_axes*1/4

    axes_position = [width_label, 0.5+bottom_axes, width_2D_axes, heigh_axes]  # left, bottom, width, height
    ax1 = fig.add_axes(axes_position)
    ax.append(ax1)

    axes_position = [width_2D_axes+width_label, 0.5+bottom_axes, width_line_axes, heigh_axes]  # left, bottom, width, height
    ax2 = fig.add_axes(axes_position)
    ax2.set_yticklabels([])  
    ax2.yaxis.set_ticks_position('none') 
    ax.append(ax2)

    axes_position = [left_bar, 0.5+bottom_bar, width_bar, heigh_bar]  # left, bottom, width, height
    ax1_1 = fig.add_axes(axes_position)
    ax.append(ax1_1)

    axes_position = [0.5+width_label, 0.5+bottom_axes, width_2D_axes, heigh_axes]  # left, bottom, width, height
    ax3 = fig.add_axes(axes_position)
    ax.append(ax3)

    axes_position = [0.5+width_2D_axes+width_label, 0.5+bottom_axes, width_line_axes, heigh_axes]  # left, bottom, width, height
    ax4 = fig.add_axes(axes_position)
    ax4.set_yticklabels([])  
    ax4.yaxis.set_ticks_position('none') 
    ax.append(ax4)

    axes_position = [0.5+left_bar, 0.5+bottom_bar, width_bar, heigh_bar]  # left, bottom, width, height
    ax3_1 = fig.add_axes(axes_position)
    ax.append(ax3_1)

    axes_position = [width_label, bottom_axes, width_2D_axes, heigh_axes]  # left, bottom, width, height
    ax5 = fig.add_axes(axes_position)
    ax.append(ax5)

    axes_position = [width_2D_axes+width_label, bottom_axes, width_line_axes, heigh_axes]  # left, bottom, width, height
    ax6 = fig.add_axes(axes_position)
    ax6.set_yticklabels([])  
    ax6.yaxis.set_ticks_position('none') 
    ax.append(ax6)

    axes_position = [left_bar, bottom_bar, width_bar, heigh_bar]  # left, bottom, width, height
    ax5_1 = fig.add_axes(axes_position)
    ax.append(ax5_1)

    axes_position = [0.5+width_label, bottom_axes, width_2D_axes, heigh_axes]  # left, bottom, width, height
    ax7 = fig.add_axes(axes_position)
    ax.append(ax7)

    axes_position = [0.5+width_2D_axes+width_label, bottom_axes, width_line_axes, heigh_axes]  # left, bottom, width, height
    ax8 = fig.add_axes(axes_position)
    ax8.set_yticklabels([])  
    ax8.yaxis.set_ticks_position('none') 
    ax.append(ax8)

    axes_position = [0.5+left_bar, bottom_bar, width_bar, heigh_bar]  # left, bottom, width, height
    ax7_1 = fig.add_axes(axes_position)
    ax.append(ax7_1)

    return fig, ax





if  __name__=='__main__':

    # 0, 0ne day; 1, multi days
    accutype = 1
    highview_even = 1000


    dir_current = os.path.dirname(os.path.abspath(__file__))
    dir_data = os.path.split(dir_current)[0]

    dir_cache_dawn = os.path.join(dir_data, 'data_dawn')
    dir_cache_dusk = os.path.join(dir_data, 'data_dusk')
    
    al.createfolder(dir_cache_dawn)
    al.createfolder(dir_cache_dusk)

    dates_dawn = ['20231116','20231117','20231118','20231120','20231121','20231209','20231210',
                '20240102','20240104','20240106','20240107','20240108','20240109','20240114',
                '20240115','20240211']
    dates_dusk = ['20231120','20240107','20240222','20240402']

    al.data_cache(dir_cache_dawn, dates_dawn, dir_data)
    al.data_cache(dir_cache_dusk, dates_dusk, dir_data)

    if accutype == 1:
        date_dawn='2023-11-16TO2024-02-11'
        date_dusk='2023-11-20TO2024-04-02'
    if accutype == 0:
        date_dawn='20240106'
        date_dusk='20240402'
        
    figname=['night','morning']


    fig, ax = fig_axes()    

    plt.rcParams['font.family'] = 'Times New Roman'

    twilight_type=1
    main(accutype, date_dawn, dir_cache_dawn, date_dawn+figname[twilight_type], twilight_type, fig, ax, sza_end_2D=130)   

    twilight_type=0
    main(accutype, date_dusk, dir_cache_dusk, date_dusk+figname[twilight_type], twilight_type, fig, ax, sza_end_2D=140)

    plt.close(fig)
        

    