import numpy as np
import matplotlib.pyplot as plt
import your
import driftrate
from driftrate import scilabel
from matplotlib.gridspec import GridSpec
import matplotlib.colors as colors
from matplotlib.patches import Rectangle, Ellipse
import scipy, glob
from itertools import zip_longest
from tqdm import tqdm
import pandas as pd
# from sklearn.mixture import GaussianMixture

# Based on https://github.com/mef51/subdriftlaw/blob/master/ArrivalTimes.ipynb

def line_model(nu, dtdnu):
	return dtdnu * nu

def gauss_model(x, a, xo, sigma):
	return a*np.exp(-(x-xo)**2/(2*(sigma**2)))

def gaussmix_model(x, *p):
	n = len(p)/3
	if n == 1:
		raise "Use gauss_model instead of gauss_model_n for n=1"
	if n == 4:
		a1, a2, a3, a4, xo1, xo2, xo3, xo4, sigma1, sigma2, sigma3, sigma4 = p
		return (gauss_model(x, a1, xo1, sigma1) +
		gauss_model(x, a2, xo2, sigma2) +
		gauss_model(x, a3, xo3, sigma3) +
		gauss_model(x, a4, xo4, sigma4))
	# etc..

###### N component models. Can these be generated?
# See GaussianMixture models (sklearn, lmfit, general)
def gauss2_model(
	x,
	a1, a2,
	xo1, xo2,
	sigma1, sigma2
):
	return (gauss_model(x, a1, xo1, sigma1) + gauss_model(x, a2, xo2, sigma2))

def gauss3_model(
	x,
	a1, a2, a3,
	xo1, xo2, xo3,
	sigma1, sigma2, sigma3
):
	return (gauss_model(x, a1, xo1, sigma1) +
	gauss_model(x, a2, xo2, sigma2) +
	gauss_model(x, a3, xo3, sigma3))

def gauss4_model(
	x,
	a1, a2, a3, a4,
	xo1, xo2, xo3, xo4,
	sigma1, sigma2, sigma3, sigma4
):
	# n = 4 => 3*n parameters
	return (gauss_model(x, a1, xo1, sigma1) +
	gauss_model(x, a2, xo2, sigma2) +
	gauss_model(x, a3, xo3, sigma3) +
	gauss_model(x, a4, xo4, sigma4))

def gauss5_model(
	x,
	a1, a2, a3, a4, a5,
	xo1, xo2, xo3, xo4, xo5,
	sigma1, sigma2, sigma3, sigma4, sigma5
):
	return (gauss_model(x, a1, xo1, sigma1) +
	gauss_model(x, a2, xo2, sigma2) +
	gauss_model(x, a3, xo3, sigma3) +
	gauss_model(x, a4, xo4, sigma4) +
	gauss_model(x, a5, xo5, sigma5))

def gauss6_model(
	x,
	a1, a2, a3, a4, a5, a6,
	xo1, xo2, xo3, xo4, xo5, xo6,
	sigma1, sigma2, sigma3, sigma4, sigma5, sigma6
):
	return (gauss_model(x, a1, xo1, sigma1) +
	gauss_model(x, a2, xo2, sigma2) +
	gauss_model(x, a3, xo3, sigma3) +
	gauss_model(x, a4, xo4, sigma4) +
	gauss_model(x, a5, xo5, sigma5) +
	gauss_model(x, a6, xo6, sigma6))
########### End n-component models

def fitgauss(data, duration):
	# use curve-fit (non-linear leastsq)
	if np.max(data) != 0:
		data = data / np.max(data) # normalize
	x = np.linspace(0, duration, num=len(data))*1000 # times in ms
	xo = sum(x*data)/sum(data)
	popt, pcov = scipy.optimize.curve_fit(
		gauss_model,
		x,
		data,
		p0=[
			np.max(data),
			xo,
			np.sqrt(abs(sum(data*(x-xo)**2)/sum(data))) # sigma
		]
	)
	return popt, pcov

def fit4gauss(data, duration, xos):
	if len(xos) != 4:
		raise "Please specify 4 peak times"
	if np.max(data) != 0:
		data = data / np.max(data) # normalize
	x = np.linspace(0, duration, num=wfall.shape[1])*1000 # times in ms
	sigma = np.sqrt(abs(sum(data*(x-np.mean(xos))**2)/sum(data)))/4
	guess = [*[np.max(data)]*4, *xos, *[sigma]*4]
	popt, pcov = scipy.optimize.curve_fit(gauss4_model, x, data, p0=guess)
	return popt, pcov

def plotburst(data, band, retfig=False, extent=None):
	fig, axs = plt.subplot_mosaic(
		'''
		T.
		WB
		''',
		figsize=(8, 7),
		width_ratios=[3,1],
		height_ratios=[1,3],
		gridspec_kw={
			'hspace': 0,
			'wspace': 0
		}
	)
	axs['W'].imshow(
		data,
		aspect='auto',
		origin='lower',
		interpolation='none',
		extent=extent
	)
	if not extent:
		extent = [0, data.shape[1], 0, data.shape[0]]
	axs['T'].plot(np.linspace(*extent[:2], num=data.shape[1]), np.nanmean(data, axis=0))
	axs['B'].stairs(band, np.linspace(*extent[2:], num=len(band)+1), orientation='horizontal')
	if extent:
		axs['W'].set_xlabel('Time (ms)')
		axs['W'].set_ylabel('Frequency (MHz)')

	axs['B'].yaxis.set_tick_params(labelleft=False)

	axs['T'].sharex(axs['W'])
	axs['B'].sharey(axs['W'])
	if retfig:
		return fig, axs
	else:
		plt.show()
		plt.close()

def measureburst():
	""" Compute the inverse sub-burst slope (dt/dnu) using the per-row arrival time method

	If multiple components are present, split them up and measure individually
	"""
	pass

def measuresubbursts():
	""" Mult-component version of `measureburst` """
	pass

if __name__ == '__main__':
	# files = glob.glob('/Users/mchamma/dev/SurveyFRB20121102A/data/snelders2023/figure_1/*.npz')
	# files = [
	# 	'/Users/mchamma/dev/SurveyFRB20121102A/data/snelders2023/figure_1/burst_B43.npz',
	# 	'/Users/mchamma/dev/SurveyFRB20121102A/data/snelders2023/figure_1/burst_B44.npz',
	# 	'/Users/mchamma/dev/SurveyFRB20121102A/data/snelders2023/figure_1/burst_B06_a.npz',
	# 	'/Users/mchamma/dev/SurveyFRB20121102A/data/snelders2023/figure_1/burst_B06_b.npz',
	# 	'/Users/mchamma/dev/SurveyFRB20121102A/data/snelders2023/figure_1/burst_B30.npz',
	# 	'/Users/mchamma/dev/SurveyFRB20121102A/data/snelders2023/figure_1/burst_B31_a.npz',
	# 	'/Users/mchamma/dev/SurveyFRB20121102A/data/snelders2023/figure_1/burst_B31_b.npz',
	# 	'/Users/mchamma/dev/SurveyFRB20121102A/data/snelders2023/figure_1/burst_B10.npz',
	# 	'/Users/mchamma/dev/SurveyFRB20121102A/data/snelders2023/figure_1/burst_B38.npz',
	# 	'/Users/mchamma/dev/SurveyFRB20121102A/data/snelders2023/figure_1/burst_B07.npz',
	# ]

	files = [
		'/Users/mchamma/dev/frbdata/FRB20220912A/sheikh2023/npzs/fil_59883_37405_244758056_FRB20220912a_0001_lob_10sec_crop.npz',
	]
	# files = glob.glob('/Users/mchamma/dev/frbdata/FRB20220912A/sheikh2023/npzs/*.npz')

	masks = [
		# [(356,371), (928,1247)],
	]

	results = []
	for filename, masklist in zip_longest(files, masks):
		bname = filename.split('/')[-1].split('.')[0]
		data = np.load(filename, allow_pickle=True)
		wfall = data['wfall']
		wfall = driftrate.subtractbg(wfall, 0, int(wfall.shape[1]*0.1))
		wfall = driftrate.subsample(wfall, wfall.shape[0]//6, wfall.shape[1]//2) # for 4 component burst

		# masking
		if masklist != None:
			for mask in masklist:
				wfall[mask[0]:mask[1],...] = 0

		freqs_bin0 = min(data['dfs'])
		res_freq = data['bandwidth'] / wfall.shape[0] # MHz
		res_time = data['duration'] / wfall.shape[1] # seconds
		res_time_ms = res_time*1000 # convenience
		duration = wfall.shape[1]*res_time

		freqs = np.linspace(freqs_bin0, max(data['dfs']), num=wfall.shape[0]) # channel width/2 is already added
		times = np.linspace(0, duration, num=wfall.shape[1]) # array of timestamps, in seconds
		times_ms = times*1000 # convenience
		tseries = np.nanmean(wfall, axis=0)

		def fitrows(wfall):
			fitdata = np.zeros((wfall.shape[0], 10))
			for i, row in enumerate(wfall):
				try:
					popt, pcov = fitgauss(row, wfall.shape[1]*res_time)
					perr = np.sqrt(np.diag(pcov))
					sigma = abs(popt[2])
					tstart = (popt[1]-np.sqrt(2)*sigma)
					tstart_err = np.sqrt(perr[1]**2 + 2*perr[2]**2)
					tend   = (popt[1]+np.sqrt(2)*sigma)
					fitdata[i,:] = [freqs[i], tstart, tend, popt[0], popt[1], tstart_err, sigma, *perr]
				except RuntimeError:
					continue

			return pd.DataFrame(data=fitdata, columns=[
				'freqs',
				'tstart',
				'tend',
				'amp',
				'xo',
				'tstart_err',
				'sigma',
				'amp_err',
				'xo_err',
				'sigma_err'
			])

		def measuredtdnu(pointdf, xo, tpoint='tstart'):
			dtdnu, pcov = scipy.optimize.curve_fit(
				line_model,
				pointdf['freqs'],
				pointdf[tpoint] - xo,
				sigma=pointdf[f'{tpoint}_err'],
				absolute_sigma=True,
			)
			dtdnu, dtdnu_err = dtdnu[0], np.sqrt(np.diag(pcov))[0]
			return dtdnu, dtdnu_err

		tpoint = 'tstart' # 'tend', 'xo'

		# 1 burst:
		pktime = np.nanargmax(np.nanmean(wfall, axis=0))*res_time*1000
		t_popt, _ = fitgauss(tseries, duration)
		window = 2*abs(t_popt[2]) # 2*stddev of a guassian fit to the integrated time series

		##### 4 component model: use 4 1d gaussians to make 4 windows in time,
		# then use the time windows to make frequency windows
		xos = [9.26, 15.75, 22.61, 27.94] # ms
		xos_chans = np.floor(np.array(xos)/res_time_ms)
		t4_popt, _ = fit4gauss(tseries, duration, xos=xos)
		windows4 = 2*np.abs(t4_popt[8:]) # 2*sigma
		subfalls = []
		subbands = []
		bandpass = np.zeros(wfall.shape[0])
		for xoi, s in zip(xos_chans, t4_popt[8:]):
			s4 = np.floor(4*np.abs(s)/res_time_ms)
			s2 = np.floor(4*np.abs(s)/res_time_ms)
			subfall = wfall[..., int(xoi-s4):int(xoi+s4)] # 4sigma window around burst
			subband = wfall[..., int(xoi-s2):int(xoi+s2)].mean(axis=1) # sum only 2sigma from burst peak
			subband = subband/np.max(subband) # Normalize
			bandpass += subband
			subfall = driftrate.subtractbg(subfall, 0, int(subfall.shape[1]*0.1)) # again
			subfalls.append(subfall)
			subbands.append(subband)

		#### Fitting (4 subburst model)
		dtdnus, subdfs = [], []
		subpeaks, subbandmodels = [], []
		for subfall, subband, xosi, win4i in zip(subfalls, subbands, xos, windows4):
			subdf = fitrows(subfall) # Fit a 1d gaussian in each row of the waterfall
			subpktime = np.nanargmax(np.nanmean(subfall, axis=0))*res_time*1000

			# Fit 1d gauss to burst spectrum
			fo = sum(freqs*subband)/sum(subband)
			subband_popt, _ = scipy.optimize.curve_fit(
				gauss_model,
				freqs,
				subband,
				p0=[
					np.max(subband),
					fo,
					np.sqrt(abs(sum(subband*(freqs-fo)**2)/sum(subband))) # sigma
				]
			)
			bandwin = 3*subband_popt[2] # 3sigma spectral window, to leave room for uncertainty
			pkfreq = subband_popt[1]

			## Apply time and spectral filters to points
			subdf = subdf[(subdf.amp > 0)]
			subdf = subdf[(subpktime-win4i < subdf[tpoint]) & (subdf[tpoint] < subpktime+win4i)]
			subdf = subdf[(pkfreq-bandwin < subdf['freqs']) & (subdf['freqs'] < pkfreq+bandwin)]

			dtdnu, dtdnu_err = measuredtdnu(subdf, subpktime)

			# Sub-burst plot
			subfig, subaxs = plotburst(
				subfall,
				subband.reshape(-1, 4).mean(axis=1),
				retfig=True,
				extent=[
					0,
					res_time_ms*subfall.shape[1],
					freqs_bin0,
					freqs_bin0 + res_freq*wfall.shape[0]
				]
			)
			subcolors = [(1, 1, 1, alpha) for alpha in np.clip(subdf['amp'], 0, 1)]
			subaxs['W'].scatter(
				(subdf[tpoint]),
				(subdf['freqs']),
				c=subcolors,
				marker='o',
				s=25
			)
			subaxs['W'].set_xlim(0, res_time_ms*subfall.shape[1])
			subaxs['W'].set_ylim(freqs_bin0, freqs_bin0 + res_freq*wfall.shape[0])
			subtimes = np.linspace(0, res_time_ms*subfall.shape[1], num=1000)
			subaxs['W'].plot(subtimes, (1/dtdnu)*(subtimes-subpktime), 'w--', label=f'{dtdnu=:.2e} ms/MHz')
			subaxs['W'].legend()

			subaxs['B'].plot(
				gauss_model(freqs, *subband_popt),
				freqs
			)
			subbandmodels.append(gauss_model(freqs, *subband_popt))
			# plt.show()
			plt.close()

			subdf[tpoint] = subdf[tpoint] + (xosi-2*win4i) # transform to full waterfall times

			dtdnus.append((dtdnu, dtdnu_err))
			subpeaks.append(subpktime)
			subdfs.append(subdf)
			print(f"{dtdnu = } +/- {dtdnu_err = }")
		subdf = pd.concat(subdfs)

		####### Fitting (1 burst)
		print(f"{bname} {wfall.shape = }")
		arrtimesdf = fitrows(wfall)
		print(arrtimesdf.head(wfall.shape[0]-1)[[tpoint]])
		arrtimesdf.to_csv('snelders/test.tmp.csv')
		arrtimesdf = arrtimesdf[(arrtimesdf.amp > 0)]
		arrtimesdf = arrtimesdf[(pktime-window < arrtimesdf[tpoint]) & (arrtimesdf[tpoint] < pktime+window)]

		# xobar = arrtimesdf[tpoint].mean() # The average xo AFTER filtering
		xobar = pktime
		print(f"{pktime = } {window = } {xobar = }")

		dtdnu_start, dtdnu_start_err = measuredtdnu(arrtimesdf, xobar)
		print(f"{dtdnu_start = } +/- {dtdnu_start_err} ms/MHz")

		results.append([bname, data['DM'], dtdnu_start, dtdnu_start_err])

		##### Plotting
		extent = [
			-xobar,
			res_time*1000*wfall.shape[1]-xobar,
			freqs_bin0,
			freqs_bin0 + res_freq*wfall.shape[0]
		]

		fig, axs = plt.subplot_mosaic(
			'''
			T.
			AS
			AS
			AS
			EE
			''',
			figsize=(10, 8),
			width_ratios=[2,1],
			# gridspec_kw={'hspace':0.464}
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
		ax_wfall.scatter(
			arrtimesdf[tpoint]-xobar,
			arrtimesdf['freqs'],
			c=pt_colors,
			marker='o',
			s=25,
			alpha=0.0
		)
		subcolors = [(1, 1, 1, alpha) for alpha in np.clip(subdf['amp'], 0, 1)]
		ax_wfall.scatter( # component fit points
			subdf[tpoint]-xobar,
			subdf['freqs'],
			c=subcolors,
			marker='o',
			s=25
		)
		ax_wfall.set_xlabel("Time (ms)")
		ax_wfall.set_ylabel("Frequency (MHz)")
		# ax_wfall.plot(times_ms-xobar, (1/dtdnu_start)*(times_ms - xobar), 'w--', alpha=0.75)

		# component lines
		for (dtdnu, dtdnu_err), xoi in zip(dtdnus, xos):
			ax_wfall.plot(
				times_ms-xobar,
				(1/dtdnu)*(times_ms-xoi),
				'w--',
				alpha=0.75,
				# label=f'$dt/d\\nu = $ {dtdnu:.2e} $\\pm$ {dtdnu_err:.2e}'
				label=f'{xos.index(xoi)+1}. $dt/d\\nu =$ {scilabel(dtdnu, dtdnu_err)} ms/MHz'
			)

		ax_wfall.set_title(f"{bname}")
		ax_wfall.legend(loc=1, handlelength=0)

		ax_tseries = axs['T']
		ax_tseries.plot(times_ms-xobar, tseries)
		ax_tseries.add_patch(Rectangle(
			(-window, ax_tseries.get_ylim()[0]),
			width=2*window,
			height=np.max(tseries)*0.075,
			color='tomato',
			alpha=0.0
		))

		for w, xoi in zip(windows4, xos):
			ax_tseries.add_patch(Rectangle(
				(xoi-xobar-w, ax_tseries.get_ylim()[0]),
				width=2*w,
				height=np.max(tseries)*0.075,
				color='tomato',
				alpha=0.5
			))

		ax_tseries.plot(
			times_ms-xobar,
			gauss_model(
				times_ms-xobar,
				np.max(tseries)*t_popt[0],
				t_popt[1]-xobar,
				t_popt[2]
			),
			'k--',
			alpha=0.1
		)

		# 4 Gaussian model
		t4_popt[0] *= np.max(tseries)
		t4_popt[1] *= np.max(tseries)
		t4_popt[2] *= np.max(tseries)
		t4_popt[3] *= np.max(tseries)
		t4_popt[4] -= xobar
		t4_popt[5] -= xobar
		t4_popt[6] -= xobar
		t4_popt[7] -= xobar
		ax_tseries.plot(
			times_ms-xobar,
			gauss4_model(
				times_ms-xobar,
				*t4_popt
			),
			'k--',
			alpha=0.8
		)

		### Summed Spectrum (summed over burst widths). Total and individual
		downfactor = 4
		bandpass_down = bandpass.reshape(-1, downfactor).mean(axis=1)
		axs['S'].stairs(
			bandpass_down/np.max(bandpass_down),
			np.linspace(*extent[2:], num=len(bandpass_down)+1),
			orientation='horizontal',
			# lw=2
		)
		for subbandmodel in subbandmodels:
			axs['S'].plot(
				subbandmodel,
				freqs,
				'k--',
				alpha=0.5,
				zorder=-1
			)

		axs['T'].sharex(axs['A'])
		axs['A'].sharey(axs['S'])
		axs['S'].set_xlabel('Intensity (arb.)')
		ax_wfall.set_xlim(extent[0], extent[1])
		ax_wfall.set_ylim(extent[2], extent[3])

		plt.setp(ax_tseries.get_xticklabels(), visible=False)
		# plt.setp(axs['S'].get_yticklabels(), visible=False)

		### Slope measurement plot
		ax_slope = axs['E']
		ax_slope.scatter(arrtimesdf['freqs'], arrtimesdf[tpoint]-xobar, c='k', s=20)
		ax_slope.errorbar(
			arrtimesdf['freqs'],
			arrtimesdf[tpoint]-xobar,
			yerr=arrtimesdf['tstart_err'],
			xerr=None,
			fmt='none',
			zorder=-1,
			color='#888888'
		)
		ax_slope.plot(freqs, dtdnu_start*freqs, 'k--')
		ax_slope.annotate(
			f"$dt/d\\nu =$ {scilabel(dtdnu, dtdnu_err)} ms/MHz",
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
		# plt.savefig(f"{bname}.png")
		plt.show()

	resultsdf = pd.DataFrame(data=results, columns=[
		'name',
		'DM',
		'dtdnu',
		'dtdnu_err'
	]).set_index('name')

	######## Delete this block if using for bursts other than Snelders
	# resultsdf['tau_w_ms'] = [ # copied manually from correlation spreadsheet
	# 	0.00646380671240664, # B43
	# 	0.00366414921653011, # B44
	# 	0.0116206417573761, # B06_a
	# 	0.00242298903767112, # B06_b
	# 	0.00303785178766883, # B30
	# 	0.0346724819083478, # B31_a
	# 	0.0346724819083478, # B31_b
	# 	0.00581992791903139, # B10
	# 	0.0081318875647676, # B38
	# 	0.00265068386683344, # B07
	# ]
	# resultsdf['center_f'] = [
	# 	4927.10624936093, # B43
	# 	5536.10665256537, # B44
	# 	7072.54790403014, # B06_a
	# 	7086.00586118911, # B06_b
	# 	5958.48369197416, # B30
	# 	5097.0214719115, # B31_a
	# 	5097.0214719115, # B31_b
	# 	5458.38960245244, # B10
	# 	5361.17461328926, # B38
	# 	4798.01572821007 # B07
	# ]
	# ######## end delete

	# resultsdf.to_csv('snelders/snelders2023_arrtimes_results.csv')

