import numpy as np
import matplotlib.pyplot as plt
import your
import driftrate
from matplotlib.gridspec import GridSpec
import matplotlib.colors as colors
import scipy
from tqdm import tqdm
import pandas as pd

# Based on https://github.com/mef51/subdriftlaw/blob/master/ArrivalTimes.ipynb
filename = '/Users/mchamma/dev/SurveyFRB20121102A/data/snelders2023/figure_1/burst_B30.npz'
# filename = '/Users/mchamma/dev/SurveyFRB20121102A/data/snelders2023/figure_1/burst_B43.npz'
# filename = '/Users/mchamma/dev/SurveyFRB20121102A/data/snelders2023/figure_1/burst_B07.npz'
bname = filename.split('/')[-1].split('.')[0]
data = np.load(filename, allow_pickle=True)
wfall = data['wfall']

freqs_bin0 = min(data['dfs'])
res_freq = data['bandwidth'] / wfall.shape[0] # MHz
res_time = data['duration'] / wfall.shape[1] # seconds
duration = wfall.shape[1]*res_time

freqs = np.linspace(freqs_bin0, max(data['dfs']), num=wfall.shape[0]) # channel width/2 is already added
times = np.linspace(0, duration, num=wfall.shape[1]) # array of timestamps, in seconds

def line(nu, dtdnu):
	return dtdnu * nu

def gaussian(x, a, xo, sigma):
	return a*np.exp(-(x-xo)**2/(2*(sigma**2)))

def fitgaussian(data, duration):
	# use curve-fit (non-linear leastsq)
	data = data / np.max(data) # normalize
	x = np.linspace(0, duration, num=wfall.shape[1])*1000 # times in ms
	xo = sum(x*data)/sum(data)
	sigma = np.sqrt(abs(sum(data*(x-xo)**2)/sum(data)))
	guess = [np.max(data), xo, sigma]
	popt, pcov = scipy.optimize.curve_fit(gaussian, x, data, p0=guess)
	return popt, pcov

####### Fitting
print(f"{wfall.shape = }")
fitdata = np.zeros((wfall.shape[0], 9))
for i, row in enumerate(wfall):
	try:
		popt, pcov = fitgaussian(row, duration)
		perr = np.sqrt(np.diag(pcov))
		sigma = abs(popt[2])
		tstart = (popt[1]-np.sqrt(2)*sigma)
		tstart_err = np.sqrt(perr[1]**2 + 2*perr[2]**2)
		tend   = (popt[1]+np.sqrt(2)*sigma)
		fitdata[i,:] = [freqs[i], tstart, tend, popt[0], popt[1], sigma, *perr]
	except RuntimeError:
		continue

arrtimesdf = pd.DataFrame(data=fitdata, columns=[
	'freqs',
	'tstart',
	'tend',
	'amp',
	'xo',
	'sigma',
	'amp_err',
	'xo_err',
	'sigma_err'
])
print(arrtimesdf.head(wfall.shape[0]-1)[['xo']])
arrtimesdf.to_csv('test.csv')
arrtimesdf = arrtimesdf[(arrtimesdf.amp > 0)]
pktime = np.nanargmax(np.nanmean(wfall, axis=0))*res_time*1000
window = 0.02 # ms
arrtimesdf = arrtimesdf[(pktime-window < arrtimesdf.xo) & (arrtimesdf.xo < pktime+window)] #B30
# arrtimesdf = arrtimesdf[(0.16 < arrtimesdf.xo) & (arrtimesdf.xo < 0.1975)] #B43
# arrtimesdf = arrtimesdf[(0.25 < arrtimesdf.xo) & (arrtimesdf.xo < 0.3)] #B07

xobar = arrtimesdf['xo'].mean() # The average xo AFTER filtering

dtdnu_start, pcov = scipy.optimize.curve_fit(
	line,
	arrtimesdf['freqs'],
	arrtimesdf['tstart'] - xobar,
	sigma=arrtimesdf['xo_err'],
	absolute_sigma=True,
)
dtdnu_start, dtdnu_start_err = dtdnu_start[0], np.sqrt(np.diag(pcov))[0]
print(f"{dtdnu_start = } +/- {dtdnu_start_err} ms/MHz")

dtdnu_end, pcov = scipy.optimize.curve_fit(
	line,
	arrtimesdf['freqs'],
	arrtimesdf['tend'] - xobar,
	sigma=arrtimesdf['xo_err'],
	absolute_sigma=True,
)
dtdnu_end, dtdnu_end_err = dtdnu_end[0], np.sqrt(np.diag(pcov))[0]
print(f"{dtdnu_end = } +/- {dtdnu_end_err} ms/MHz")

dtdnu, pcov = scipy.optimize.curve_fit(
	line,
	arrtimesdf['freqs'],
	arrtimesdf['xo'] - xobar,
	sigma=arrtimesdf['xo_err'],
	absolute_sigma=True,
)
dtdnu, dtdnu_err = dtdnu[0], np.sqrt(np.diag(pcov))[0]
print(f"{dtdnu = } +/- {dtdnu_err} ms/MHz")

##### Plotting
extent = [
	-xobar,
	res_time*1000*wfall.shape[1]-xobar,
	freqs_bin0,
	freqs_bin0 + res_freq*wfall.shape[0]
]

fig, axs = plt.subplot_mosaic(
	'''
	AB
	AC
	AD
	EE
	''',
	figsize=(10, 8),
	width_ratios=[2,1],
)

### Waterfall
ax_wfall = axs['A']
ax_wfall.imshow(
	wfall,
	aspect='auto',
	origin='lower',
	interpolation='none',
	extent=extent,
	norm='linear',
	vmax=np.quantile(wfall, 0.999),
)
pt_colors = [(1, 1, 1, alpha) for alpha in np.clip(arrtimesdf['amp'], 0, 1)]
ax_wfall.scatter(arrtimesdf['xo']-xobar, arrtimesdf['freqs'], c=pt_colors, marker='o', s=25)
ax_wfall.set_xlabel("Time (ms)")
ax_wfall.set_ylabel("Frequency (MHz)")
ax_wfall.plot(times*1000-xobar, (1/dtdnu)*(times*1000 - xobar), 'w--', alpha=0.75)

m = (5613 - 6452)/(0.0774 - 0.0756)
b = 6452 - 0.0756*m
# ax_wfall.plot(times*1000, (m)*times*1000+b, 'w--')

### Individual Profiles
profiles = arrtimesdf.sort_values('amp', ascending=False).head(3)
chans = profiles.index
lines = [np.interp(chan, range(0, wfall.shape[0]), freqs) for chan in chans]
for line in lines:
	ax_wfall.axhline(
		y=line,
		c='w',
		linewidth=0.75,
		ls='-.'
	)

for irow, panel in zip(profiles.iterrows(), ['B', 'C', 'D']):
	chan, row = irow
	x = np.linspace(0, duration, num=wfall.shape[1])*1000 # times in ms
	xm = np.linspace(0, duration, num=2000)*1000 # times in ms
	axs[panel].plot(x-xobar, wfall[chan,:]/np.max(wfall)) # match fit normalization
	axs[panel].axvline(x=row['xo']-xobar, ls='--', c='k')
	axs[panel].set_title(f'{row["freqs"]:.1f} MHz Spectrum', fontsize='medium')
	# axs[panel].plot(xm, row['amp']*np.exp(-(xm-row['xo'])**2/(2*(row['sigma']**2))), 'k--')

axs['A'].sharex(axs['B'])
axs['B'].sharex(axs['C'])
axs['C'].sharex(axs['D'])
# axs['B'].set_xticklabels([])
# axs['C'].set_xticklabels([])
axs['D'].set_xlabel('Time (ms)')
axs['B'].set_xlim(row['xo']-xobar - 0.015, row['xo']-xobar + 0.015)
ax_wfall.set_xlim(extent[0], extent[1])
ax_wfall.set_ylim(extent[2], extent[3])

### Slope measurement plot
ax_slope = axs['E']
ax_slope.scatter(arrtimesdf['freqs'], arrtimesdf['xo']-xobar, c='k', s=20)
ax_slope.errorbar(
	arrtimesdf['freqs'],
	arrtimesdf['xo']-xobar,
	yerr=arrtimesdf['xo_err'],
	xerr=None,
	fmt='none',
	zorder=-1,
	color='#888888'
)
ax_slope.plot(freqs, dtdnu*freqs, 'k--')
ax_slope.annotate(
	f"$dt/d\\nu = $ {dtdnu:.1e} $\\pm$ {dtdnu_err:.1e} ms/MHz",
	xy=(0.675, 0.8),
	xycoords='axes fraction',
	color='white',
	weight='black',
	size=10,
	bbox={"boxstyle":"round"}
)
ax_slope.set_xlabel("Frequency (MHz)")
ax_slope.set_ylabel("Time (ms)")
ax_slope.set_xlim(min(arrtimesdf['freqs']), max(arrtimesdf['freqs']))
plt.tight_layout()
plt.show()
