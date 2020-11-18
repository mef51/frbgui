"""
Originally written by Kenzie Nimmo 2020
https://github.com/KenzieNimmo/FRB_filterbank_tools
Adapted by Mark Snelders for VERY specific purposes
"""
"""
WARNING: EVEN THOUGH THIS SCRIPT LOOKS LIKE IT IS FLEXIBLE, IT ONLY WORKS
IN IT'S CURRENT STATE. CHANGING THINGS WILL EITHER BREAK IT OR GIVE IN-
CORRENT VALUES. THIS IS BECAUSE SOME VALUES ARE HARDCODED.
"""

# import modules
import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from matplotlib.ticker import MultipleLocator
import pickle
import sys
import matplotlib.patheffects as PathEffects

# set matplotlib variables
mpl.rcParams['font.size'] = 10
mpl.rcParams['axes.linewidth'] = 1
# mpl.rcParams['ps.useafm'] = True
# mpl.rcParams['pdf.use14corefonts'] = True
# mpl.rcParams['text.usetex'] = True
mpl.rcParams['legend.fontsize'] = 10
mpl.rcParams['axes.labelsize'] = 10
mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10
# mpl.rcParams['text.latex.preamble'] = '\\usepackage[T1]{fontenc}'

def plot(pkl,bind,favg,tavg,width,fig,plot_grid_idx,index,threshold,plot_spectrum=False,addmask=None,stokesI=None,linear=None,circular=None,tavgpol=1, acfwidth=None, acfshift_lc=0, acfshift_dc=0, scalefwhm=1, plot_PPA=False, PPA=None):
    """
    Function which generates the figure for each individual burst. pkl is the pklfile containing
    the burst information. bind = 1 for B1, and 2 for B2. In our case favg should be one since we
    do not downsample in frequency. tavg = 4 since we downsample from 8 us data to 32 us data.
    width was the width as determined by the 2D gaussian method, but is no longer being used.
    tavgpol = 1 since the polarization profiles (I/L/V_B1/B2.npy) have 32 us time res. acfshift_lc
    and dc (light cyan and dark cyan) are integers and denote how many bins they need to be shifted with
    respect to the peak of the burst profile (= stokes I).
    """
    print(pkl)
    with open(pkl,'rb') as f:
        pklfile=pickle.load(f, encoding='latin-1')
    f.close()

    arr_corr=pklfile['array_corrected']
    mask = pklfile['mask']
    res_t = pklfile['t_samp']
    freqs = pklfile['freqs']
    nchan_tot = len(freqs)
    res_f = np.abs(freqs[0] - freqs[1])
    print(res_f)
    fmin = freqs[0]-(res_f/2.)
    fmax = freqs[-1]+(res_f/2.)
    print(fmin, fmax)
    print(len(freqs), freqs)

    t_cent = pklfile['centre_bin']
    t_fwhm = pklfile['width_bin']
    prof_flux = pklfile['profile']
    spec_flux = pklfile['spectrum']
    spec_flux[mask]=0
    nbins = np.shape(arr_corr)[1]

    conv = 2.3548 # convert FWHM to sigma, 2 * sqrt(2 * (ln(2)))

    if favg == None:
        favg = int(1)
    if tavg == None:
        tavg = int(1)

    # use the corrected array (bandpass and subbandedges)
    arr=arr_corr

    spectrum=np.mean(arr,axis=1)
    if addmask != None:
        for i in range(len(addmask)):
            mask=np.append(mask,addmask[i])
    arr[mask]=0

    # make copy of the array so that I can safely calculate the S/N stuff for time
    tarr = arr.copy()
    # make copy of the array so that I can safely calculate the S/N stuff for frequency
    farr = arr.copy()

    # determine the bounds of the offtimes
    for i, d in enumerate(np.array(pklfile["offtimes"])):
        if np.array(pklfile["offtimes"][i+1]) - d > 1:
            lb, rb = d, pklfile["offtimes"][i+1]
            break
    lb = int(lb / tavg - 1)
    rb = int(rb / tavg + 1)
    timeseries = tarr.sum(axis=0)
    timeoff = np.concatenate((timeseries[:lb], timeseries[rb:]))
    timeseries = timeseries - np.mean(timeoff)
    timeseries = timeseries / np.std(timeoff)

    # determine where the peak of the timeseries (8 us data)
    pkprofile = np.argmax(timeseries)
    # downsample time series
    timeseries = timeseries.reshape(-1, tavg).mean(axis=1)

    # switch t_cent from 2d gauss to peak profile
    t_cent = pkprofile

    maskarray=np.ones_like(spectrum)
    maskarray[mask]=0

    #downsample in time
    if tavg>1:
        tavg=float(tavg)
        binstoremove=(nbins/tavg)-np.int(nbins/tavg)
        binstoremove*=tavg
        endind=int(nbins-binstoremove)
        arr_trim=arr[:,:endind]
        prof_flux_trim=prof_flux[:endind]
        newnbins=np.shape(arr_trim)[1]
        tavgbins=newnbins/tavg
        arr=np.array(np.column_stack([np.mean(subint, axis=1) for subint in np.hsplit(arr_trim,tavgbins)]))
        tavg=int(tavg)
        prof_flux = np.sum(prof_flux_trim.reshape(-1, tavg), axis=1)

    #downsample in frequency
    if favg>1:
        favg=float(favg)
        if (nchan_tot/favg)-int(nchan_tot/favg)!=0:
            print("The total number of channels is %s, please choose an fscrunch value that divides the total number of channels.")
            sys.exit()
        else:
            newnchan=nchan_tot/favg
            arr=np.array(np.row_stack([np.mean(subint, axis=0) for subint in np.vsplit(arr,newnchan)]))
            favg=int(favg)
            maskarray=np.mean(maskarray.reshape(-1,favg),axis=1)
            spec_flux=np.mean(spec_flux.reshape(-1,favg),axis=1)

    mask=np.where(maskarray==0)[0]
    arr[mask]=0

    #determine the time and frequency resolution of the plot based on downsampling
    res_t*=tavg
    res_f*=favg

    f_l_bin=0
    f_h_bin=nchan_tot-1

    extent=np.zeros(4) #time_min, time_max, freq_min, freq_max
    peak = np.int(np.ceil(t_cent))/tavg #time bin containing peak of burst
    extent[0] = - width / 2. # width is time window around the burst
    extent[1] = width / 2.
    extent[2] = fmin
    extent[3] = fmax
    centre = np.ceil(width / 2. / (res_t*1e3)) #bin number of the centre of window

    # cut out the relevant part of the array (within extent)
    t_fwhm/=tavg
    spectrum=spec_flux
    maskspec=maskarray
    spectrum=np.ma.masked_where(maskspec==False,spectrum)
    arr=arr[:,int(np.ceil(peak-centre)):int(np.ceil(peak+centre))]
    timeseries = timeseries[int(np.ceil(peak-centre)):int(np.ceil(peak+centre))]

    tidx = np.argmax(timeseries)
    deltatidx = tidx - len(timeseries) // 2
    ts = timeseries

    if stokesI!=None and linear!=None and circular!=None:
        polI=np.load(stokesI)
        polL=np.load(linear)
        polV=np.load(circular)
        if tavgpol>1:
            tavgpol=float(tavgpol)
            nbins = len(polI)
            binstoremove=(nbins/tavgpol)-np.int(nbins/tavgpol)
            binstoremove*=tavgpol
            endind=int(nbins-binstoremove)
            polItrim=polI[:endind]
            polLtrim=polL[:endind]
            polVtrim=polV[:endind]
            newnbins=len(polItrim)
            tavgbins=newnbins/tavgpol
            polI=np.sum(polItrim.reshape(-1,tavgpol),axis=1)
            polL=np.sum(polLtrim.reshape(-1,tavgpol),axis=1)
            polV=np.sum(polVtrim.reshape(-1,tavgpol),axis=1)

        peakpolI = np.argmax(polI) - deltatidx
        Its = polI[int(np.ceil(peakpolI-centre)):int(np.ceil(peakpolI+centre))]
        linearts = polL[int(np.ceil(peakpolI-centre)):int(np.ceil(peakpolI+centre))]
        circularts = polV[int(np.ceil(peakpolI-centre)):int(np.ceil(peakpolI+centre))]

        # Normalize polarizations
        Itsmaxval = np.max(Its)
        Its = Its / Itsmaxval
        linearts = linearts / Itsmaxval
        circularts = circularts / Itsmaxval

    if PPA != None:
        PPA_array = np.load(PPA)
        PPA_cut = PPA_array[:, int(np.ceil(peakpolI-centre)):int(np.ceil(peakpolI+centre))]

    # creating the figure structure
    rows=2
    cols=1
    if plot_spectrum:
        cols += 1

    if plot_PPA:
        rows += 1

    if plot_PPA:
        plot_grid = gridspec.GridSpecFromSubplotSpec(rows, cols, plot_grid_idx, wspace=0., hspace=0.,\
                height_ratios=[101,288,578],width_ratios=[5,]+[1,]*(cols-1))
    else:
        plot_grid = gridspec.GridSpecFromSubplotSpec(rows, cols, plot_grid_idx, wspace=0., hspace=0.,\
                height_ratios=[1,]*(rows-1)+[2,],width_ratios=[5,]+[1,]*(cols-1))

    ax1 = plt.Subplot(fig, plot_grid[rows-1,0])
    ax2 = plt.Subplot(fig, plot_grid[rows-2,0], sharex=ax1)
    if plot_spectrum:
        ax3 = plt.Subplot(fig, plot_grid[rows-1,1], sharey=ax1)
    if plot_PPA:
        ax4 = plt.Subplot(fig, plot_grid[rows-3,0], sharex=ax1)

    units = ("GHz", "ms")

    # plotting the waterfall plot with the burst at the centre, arr
    cm1 = mpl.colors.ListedColormap(['black','red'])
    cm2 = mpl.colors.ListedColormap(['black','blue'])

    # set colours limits to 1st and 99th percentile of array
    vmin = np.percentile(arr, 1)
    vmax = np.percentile(arr, 99)
    zapped = mask #identifying the frequencies that were masked

    cmap = plt.cm.viridis
    cmap.set_bad((1, 0, 0, 1)) #set color for masked values
    zap_size = int(arr.shape[1]/18)
    arr[zapped,:zap_size] = vmin-1000
    mask1=arr<vmin-600
    mask1 = np.ma.masked_where(mask1==False,mask1)
    cmap.set_over(color='white') #set color used for high out-of-range values (the zapped values have been set to NaN)
    arr[zapped,zap_size:] = vmax+1000.

    ax1.imshow(arr, cmap=cmap, origin='lower', aspect='auto', interpolation='nearest', vmin=vmin, vmax=vmax, extent=extent)

    np.save(pkl.split('.')[0] + '_downsampled', arr)

    ax1.imshow(mask1, cmap=cm1, origin='lower', aspect='auto', interpolation='nearest', vmin=0, vmax=1, extent=extent)

    if width: ax1.set_xlim(-width/2.-0.001, width/2.+0.001)
    ax1.set_ylim(fmin, fmax)

    #Label only edge plots
    if index % ncols == 0:
        ax1.set_ylabel('${\\rm Frequency}\\ ({\\rm GHz})$')#.format(units[0]))
        yticks=np.array([1.28, 1.30, 1.32, 1.34, 1.36, 1.38])
        ax1.set_yticks(yticks*1000.)
        ax1.set_yticklabels(['$1.28$','$1.30$','$1.32$','$1.34$','$1.36$','$1.38$'])
        ax2.tick_params(axis='both', direction='in')
        ax1.tick_params(axis='both', direction='in')
    else:ax1.tick_params(axis='y', labelleft='off', direction='in')

    if (index <threshold) and width:
        ax1.tick_params(axis='x', labelbottom='off', direction='in')
        ax1.tick_params(axis='both', direction='in')
    else:
        #ax1.set_xlabel('${\rm Time}\\ ({\rm ms})$')#.format(units[1]))
        xticks=([(-width/2.),(-width/4.),0,(width/4.),(width/2.)])
        ax1.set_xticks(xticks)
        ax1.tick_params(axis='x', direction='in')

    ax1.xaxis.set_minor_locator(MultipleLocator(0.5))
    ax1.tick_params(axis='both', direction='in', which='both')

    #plot pulse profile (time series)
    x = np.linspace(extent[0], extent[1], ts.shape[0])
    if stokesI!=None and linear!=None and circular!=None:
        ax2.plot(x, Its, 'k-',alpha=0.8,zorder=1.2,label=r"${\rm Total}$" + "\n" + r"${\rm Intensity}$")
        ax2.plot(x,linearts,'r-', zorder=1.1, alpha=0.6, label=r"${\rm Linear}$")
        ax2.plot(x,circularts,'b-', zorder=1.0, alpha=0.6, label=r"${\rm Circular}$")
        ax2.set_yticks(np.array([0, round(Its.max(), 1)]))
        if index % ncols == 0:
            ax2.set_ylabel('${\\rm Normalised}$ ${\\rm Flux}$ ')
        arraymaxval = [Its.max(), linearts.max(), circularts.max()]
        arrayminval = [Its.min(), linearts.min(), circularts.min()]
        y_range = np.max(arraymaxval) - np.min(arrayminval)
        ymaxax2 = np.max(arraymaxval)
        legax2 = ax2.legend(loc=1, ncol=1, frameon=True, handletextpad=0.5, handlelength=0.7)
        legax2.get_frame().set_edgecolor('k')
    else:
        ax2.plot(x, ts/np.max(ts), 'g-',alpha=1.0,zorder=1, label="TS")
        ax2.set_yticks(np.array([0,round(ts.max(),1)]))
        if index % ncols == 0:
            ax2.set_ylabel('${\\rm S/N}$')
        y_range = ts.max() - ts.min()
        ymaxax2 = ts.max()

    box = ax2.get_position()
    point = ax2.scatter(x[0], ts[0], facecolors='none', edgecolors='none')
    legend_x = 0.82
    legend_y = 0.72

    # txtb=ax4.text(0.9,0.70,('${\rm {\bf B%s}}$')%bind, ha='center', va='center', fontsize=20, transform = ax4.transAxes)
    # txtb.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='w')])

    if bind == '1':
        panels = ['a', 'c', 'e', 'f']
    elif bind == '2':
        panels = ['b', 'd', 'g', 'h']

    # txtp2=ax2.text(0.08,0.90,('${\rm {\bf %s}}$')%panels[1], ha='center', va='center', fontsize=10, transform = ax2.transAxes)
    # txtp2.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='w')])

    # txtp4=ax4.text(0.08,0.80,('${\rm {\bf %s}}$')%panels[0], ha='center', va='center', fontsize=10, transform = ax4.transAxes)
    # txtp4.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='w')])

    # txtp1=ax1.text(0.08,0.95,('${\rm {\bf %s}}$')%panels[2], ha='center', va='center', fontsize=10, transform = ax1.transAxes)
    # txtp1.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='w')])

    # txtp3=ax3.text(0.29,0.95,('${\rm {\bf %s}}$')%panels[3], ha='center', va='center', fontsize=10, transform = ax3.transAxes)
    # txtp3.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='w')])

    ax2.tick_params(axis='x', labelbottom='off', top='off', direction='in')
    ax2.set_ylim(-y_range/4.5, ymaxax2*1.1)

    if acfwidth is None:
        ax2.hlines(y=-y_range/4.5,xmin=(-t_fwhm*2./conv*res_t*1e3),xmax=(t_fwhm*2./conv*res_t*1e3), lw=10,color='#48D1CC',zorder=0.8,alpha=0.3) # light cyan
        ax2.hlines(y=-y_range/4.5, xmin=(-t_fwhm/2.*res_t*1e3),xmax=(t_fwhm/2.*res_t*1e3), lw=10, color='#48D1CC',zorder=0.9) # dark cyan
    elif acfwidth != None:
        ax2.hlines(y=-y_range/4.5,xmin=(-acfwidth*scalefwhm*res_t*1e3+acfshift_lc*1e3*res_t), xmax=(acfwidth*scalefwhm*res_t*1e3+acfshift_lc*1e3*res_t), lw=10,color='#48D1CC',zorder=0.8,alpha=0.3)  # light cyan
        ax2.hlines(y=-y_range/4.5, xmin=(-acfwidth*res_t*1e3+acfshift_dc*1e3*res_t), xmax=(acfwidth*res_t*1e3+acfshift_dc*1e3*res_t), lw=10, color='#48D1CC',zorder=0.9) # dark cyan

    ax2.tick_params(direction='in', which='both', labelbottom=False, labeltop=False, labelleft=True, labelright=False, bottom=True, top=True, left=True, right=True)

    ax2.xaxis.set_minor_locator(MultipleLocator(0.5))

    #plot spectrum (amplitude vs freq) only if plot_spectrum=True
    if plot_spectrum:
        ax3.tick_params(axis='both', direction='in')
        y = np.linspace(extent[2], extent[3], spectrum.size)
        ax3.plot(spectrum, y, 'k-',zorder=2)
        ax3.tick_params(axis='x', which='both', top='off', bottom='on', labelbottom='on')
        ax3.tick_params(axis='y', labelleft='off')
        ax3.set_ylim(fmin, fmax)
        x_range = spectrum.max() - spectrum.min()
        ax3.set_xticks(np.array([0,round(np.max(spectrum),2)]))
        ax3.set_xticklabels(['$0$','$%d$'%int(np.max(spectrum))])

        if (index <threshold) and width:
            ax3.set_xlim(-x_range/3., x_range*6./5.)
        else:
            ax3.tick_params(axis='x', pad=6.0, length=6.0)
            ax3.set_xlabel('${\\rm S/N}$')
            ax3.set_xlim(-x_range/2., x_range*6./5.)

    fig.add_subplot(ax1)
    fig.add_subplot(ax2)
    if plot_spectrum:
        fig.add_subplot(ax3)

    if plot_PPA:
        cmap = plt.cm.gist_yarg
        ax4.imshow(PPA_cut, cmap=cmap, origin="lower", aspect="auto", interpolation="nearest", extent=[extent[0], extent[1], -90, 90])
        ax4.tick_params(axis='both', which='both', direction='in', top='on', right='on', labelbottom='False')
        ax4yticks=np.array([-80, -40, 0, 40, 80])
        ax4.set_yticks(ax4yticks)
        ax4.set_yticklabels(['$-80$','$-40$','$0$','$40$','$80$'])
        ax4.yaxis.set_minor_locator(MultipleLocator(20))
        if index % ncols == 0:
            ax4.set_ylabel('${\\rm PPA}$ ${\\rm (deg)}$')

    if plot_PPA:
        fig.add_subplot(ax4)

    return

if __name__ == '__main__':
    #adapt the figure parameters as needed:
    nrows = 1
    ncols = 2 # change this to suit the number of bursts you want to plot
    threshold=(nrows*ncols)-ncols #for axes on plots
    width = 6. #Time window around the burst in ms.
    plot_grid = gridspec.GridSpec(nrows, ncols, wspace=0.1, hspace=0.1) # grid of burst plots
    fig = plt.figure(figsize=[8, 8]) # defines the size of the plot

    IDs_ordered = [1,2]

    td = 4
    fd = 1

    bursts = {"1": {
                'pkl': 'B1.pkl',
                'favg':fd,
                'tavg':td,
                'addmask':[],
                'stokesIprof':'I_B1.npy',
                'linearprof':'L_B1.npy',
                'circprof':'V_B1.npy',
                'tavgpol':1,
                'acfwidth':14,
                'acfshift_lc':13,
                'acfshift_dc':5,
                'scalefwhm':2,
                'PPA':'PPA_B1.npy'
            },
            "2": {
                'pkl': 'B2.pkl',
                'favg':fd, 'tavg':td,
                'addmask':[],
                'stokesIprof':'I_B2.npy',
                'linearprof':'L_B2.npy',
                'circprof':'V_B2.npy',
                'tavgpol':1,
                'acfwidth':15,
                'acfshift_lc':13,
                'acfshift_dc':9,
                'scalefwhm':1.5,
                'PPA':'PPA_B2.npy'}
            }
    idx = 0
    for burst in IDs_ordered:
        burst=str(burst)
        plot(bursts[burst]['pkl'],burst,bursts[burst]['favg'],bursts[burst]['tavg'],width,fig,plot_grid[idx],idx,threshold,plot_spectrum=True,addmask=bursts[burst]['addmask'],stokesI=bursts[burst]['stokesIprof'],linear=bursts[burst]['linearprof'],circular=bursts[burst]['circprof'],tavgpol=bursts[burst]['tavgpol'], acfwidth=bursts[burst]['acfwidth'], acfshift_lc=bursts[burst]['acfshift_lc'], acfshift_dc=bursts[burst]['acfshift_dc'], scalefwhm=bursts[burst]['scalefwhm'], plot_PPA=True,PPA=bursts[burst]['PPA'])
        idx+=1

    fig.subplots_adjust(hspace=0.08, wspace=0.05, left=0.09,right=.96,bottom=.1,top=.97)
    fig.savefig("bursts2.pdf".format(td), format='pdf', dpi=1200, bbox_inches='tight')
    # plt.show()
