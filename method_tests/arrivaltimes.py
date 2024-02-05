import numpy as np
import matplotlib.pyplot as plt
import your
from matplotlib.gridspec import GridSpec
import matplotlib.colors as colors
from matplotlib.patches import Rectangle, Ellipse
import scipy, glob
from itertools import zip_longest
from tqdm import tqdm
import pandas as pd
# from sklearn.mixture import GaussianMixture

import driftrate
from driftrate import scilabel, subburst_suffixes

# Based on https://github.com/mef51/subdriftlaw/blob/master/ArrivalTimes.ipynb

def line_model(nu, dtdnu):
	return dtdnu * nu

def gauss_model(x, a, xo, sigma):
	return a*np.exp(-(x-xo)**2/(2*(sigma**2)))

def gaussmix_model(x, *p):
	n = len(p)/3
	if n == 1:
		return gauss_model(x, *p)
	if n == 2:
		a1, a2, xo1, xo2, sigma1, sigma2 = p
		return (
			gauss_model(x, a1, xo1, sigma1) +
			gauss_model(x, a2, xo2, sigma2)
		)
	if n == 3:
		a1, a2, xo1, xo2, sigma1, sigma2 = p
		return (
			gauss_model(x, a1, xo1, sigma1) +
			gauss_model(x, a2, xo2, sigma2) +
			gauss_model(x, a3, xo3, sigma3)
		)
	if n == 4:
		a1, a2, a3, a4, xo1, xo2, xo3, xo4, sigma1, sigma2, sigma3, sigma4 = p
		return (
			gauss_model(x, a1, xo1, sigma1) +
			gauss_model(x, a2, xo2, sigma2) +
			gauss_model(x, a3, xo3, sigma3) +
			gauss_model(x, a4, xo4, sigma4)
		)
	if n == 5:
		a1, a2, a3, a4, xo1, xo2, xo3, xo4, sigma1, sigma2, sigma3, sigma4 = p
		return (
			gauss_model(x, a1, xo1, sigma1) +
			gauss_model(x, a2, xo2, sigma2) +
			gauss_model(x, a3, xo3, sigma3) +
			gauss_model(x, a4, xo4, sigma4) +
			gauss_model(x, a5, xo5, sigma5)
		)

	if n == 6:
		( # this is a lot
			a1, a2, a3, a4, a5, a6,
			xo1, xo2, xo3, xo4, xo5, xo6,
			sigma1, sigma2, sigma3, sigma4, sigma5, sigma6
		) = p
		return (
			gauss_model(x, a1, xo1, sigma1) +
			gauss_model(x, a2, xo2, sigma2) +
			gauss_model(x, a3, xo3, sigma3) +
			gauss_model(x, a4, xo4, sigma4) +
			gauss_model(x, a5, xo5, sigma5) +
			gauss_model(x, a6, xo6, sigma6)
		)


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
	x = np.linspace(0, duration, num=len(data))
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
	x = np.linspace(0, duration, num=len(data))
	sigma = np.sqrt(abs(sum(data*(x-np.mean(xos))**2)/sum(data)))/4
	guess = [*[np.max(data)]*4, *xos, *[sigma]*4]
	popt, pcov = scipy.optimize.curve_fit(gauss4_model, x, data, p0=guess)
	return popt, pcov

def fitgaussmix(data, duration, xos):
	n = len(xos) # Number of components
	if np.max(data) != 0:
		data = data / np.max(data) # normalize
	x = np.linspace(0, duration, num=len(data))
	sigmas = np.sqrt(abs(sum(data*(x-np.mean(xos))**2)/sum(data)))/4
	guess = [*[np.max(data)]*n, *xos, *[sigmas]*n]
	popt, pcov = scipy.optimize.curve_fit(gaussmix_model, x, data, p0=guess)
	return popt, pcov

def fitrows(wfall, dt, freqs):
	fitdata = np.zeros((wfall.shape[0], 10))
	for i, row in enumerate(wfall):
		try:
			popt, pcov = fitgauss(row, wfall.shape[1]*dt)
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

def measureburst(filename, xos=[]):
	""" Measure spectro-temporal properties of a burst

	Compute the inverse sub-burst slope (dt/dnu) using the per-row arrival time method
	Compute the duration and bandwidth by finding a 1-dimensional gaussian model
	to the integrated time series and spectrum, respectively. The duration and bandwidth are the
	1 sigma widths of either fit.
	Compute the center frequency as the center of the 1d spectrum model

	If multiple components are present, split them up and measure individually.
	"""
	results = []
	bname = filename.split('/')[-1].split('.')[0]
	data = np.load(filename, allow_pickle=True)
	wfall = data['wfall']
	wfall = driftrate.subtractbg(wfall, 0, int(wfall.shape[1]*0.1))
	wfall = driftrate.subsample(wfall, wfall.shape[0]//6, wfall.shape[1]//2) # for 4 component burst

	freqs_bin0 = min(data['dfs'])
	res_freq = data['bandwidth'] / wfall.shape[0] # MHz
	res_time_ms = 1000*data['duration'] / wfall.shape[1] # ms
	duration = wfall.shape[1]*res_time_ms

	freqs = np.linspace(freqs_bin0, max(data['dfs']), num=wfall.shape[0]) # channel width/2 is already added
	times_ms = np.linspace(0, duration, num=wfall.shape[1]) # array of timestamps
	tseries = np.nanmean(wfall, axis=0)

	tpoint = 'tstart' # 'tend', 'xo'
	pktime = np.nanargmax(tseries)*res_time_ms
	t_popt, _ = fitgauss(tseries, duration) # whether one or many components, for ref in plot
	print(f"{bname} {wfall.shape = }")

	if len(xos) == 0:
		xos.append(pktime)

	## Legacy 1 burst:
	window = 2*abs(t_popt[2]) # 2*stddev of a guassian fit to the integrated time series
	arrtimesdf = fitrows(wfall, res_time_ms, freqs)
	arrtimesdf = arrtimesdf[(arrtimesdf.amp > 0)]
	arrtimesdf = arrtimesdf[(pktime-window < arrtimesdf[tpoint]) & (arrtimesdf[tpoint] < pktime+window)]

	##### multi component model: use multiple 1d gaussians to make multiple windows in time,
	# then use the time windows to make frequency windows
	xos_chans = np.floor(np.array(xos)/res_time_ms)
	n_bursts = len(xos)

	tmix_popt, tmix_pcov = fitgaussmix(tseries, duration, xos=xos)
	tmix_perr = np.sqrt(np.diag(tmix_pcov))
	tmix_amps   = tmix_popt[:n_bursts]
	tmix_xos    = tmix_popt[n_bursts:n_bursts*2]
	tmix_sigmas = tmix_popt[n_bursts*2:n_bursts*3]
	tmix_sigma_errs = tmix_perr[n_bursts*2:n_bursts*3]

	subfalls = []
	subbands = []
	bandpass = np.zeros(wfall.shape[0])
	for xoi, s in zip(xos_chans, tmix_sigmas):
		s4 = np.floor(4*np.abs(s)/res_time_ms)
		s1 = np.floor(1*np.abs(s)/res_time_ms)
		subfall = wfall[..., int(xoi-s4):int(xoi+s4)] # 4sigma window around burst
		subband = wfall[..., int(xoi-s1):int(xoi+s1)].mean(axis=1) # sum only 1 sigma from burst peak
		bandpass += subband
		subfall = driftrate.subtractbg(subfall, 0, int(subfall.shape[1]*0.1)) # again
		subfalls.append(subfall)
		subbands.append(subband)

	#### Fitting (4 subburst model)
	dtdnus, subdfs = [], []
	subpeaks, subbandmodels = [], []
	for subfall, subband, xosi, sigma, sigma_err in zip(
		subfalls, subbands, xos, tmix_sigmas, tmix_sigma_errs
	):
		sigma = abs(sigma)
		subdf = fitrows(subfall, res_time_ms, freqs) # Fit a 1d gaussian in each row of the waterfall
		subpktime = np.nanargmax(np.nanmean(subfall, axis=0))*res_time_ms

		# Fit 1d gauss to burst spectrum
		fo = sum(freqs*subband)/sum(subband) # this is an estimate of center_f
		subband_popt, subband_pcov = scipy.optimize.curve_fit(
			gauss_model,
			freqs,
			subband/np.max(subband),
			p0=[
				1,
				fo,
				np.sqrt(abs(sum(subband*(freqs-fo)**2)/sum(subband))) # sigma
			]
		)
		subband_perr = np.sqrt(np.diag(subband_pcov))
		bwidth, bwidth_err = subband_popt[2], subband_perr[2] # sigma of spetrum fit
		pkfreq, pkfreq_err = subband_popt[1], subband_perr[1] # this is fitted center_f

		## Apply time and spectral filters to points
		subdf = subdf[(subdf.amp > 0)]
		subdf = subdf[(subpktime-2*sigma < subdf[tpoint]) & (subdf[tpoint] < subpktime+2*sigma)]
		subdf = subdf[
			(pkfreq-3*bwidth < subdf['freqs']) &
			(subdf['freqs'] < pkfreq+3*bwidth)
		]

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

		subdf[tpoint] = subdf[tpoint] + (xosi-4*sigma) # transform to full waterfall times

		dtdnus.append((dtdnu, dtdnu_err))
		subpeaks.append(subpktime)
		subdfs.append(subdf)
		print(f"{dtdnu = } +/- {dtdnu_err = }")

		results.append([
			f'{bname}_{subburst_suffixes[xos.index(xosi)]}',
			float(data['DM']),
			pkfreq,
			pkfreq_err,
			sigma,
			sigma_err,
			bwidth,
			bwidth_err,
			dtdnu,
			dtdnu_err
		])

	subdf = pd.concat(subdfs)

	##### Plotting
	extent = [
		-pktime,
		res_time_ms*wfall.shape[1]-pktime,
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
		arrtimesdf[tpoint]-pktime,
		arrtimesdf['freqs'],
		c=pt_colors,
		marker='o',
		s=25,
		alpha=0.0
	)
	subcolors = [(1, 1, 1, alpha) for alpha in np.clip(subdf['amp'], 0, 1)]
	ax_wfall.scatter( # component fit points
		subdf[tpoint]-pktime,
		subdf['freqs'],
		c=subcolors,
		marker='o',
		s=25
	)
	ax_wfall.set_xlabel("Time (ms)")
	ax_wfall.set_ylabel("Frequency (MHz)")

	# component lines
	for (dtdnu, dtdnu_err), xoi in zip(dtdnus, xos):
		ax_wfall.plot(
			times_ms-pktime,
			(1/dtdnu)*(times_ms-xoi),
			'w--',
			alpha=0.75,
			# label=f'$dt/d\\nu = $ {dtdnu:.2e} $\\pm$ {dtdnu_err:.2e}'
			label=f'{subburst_suffixes[xos.index(xoi)]}. $dt/d\\nu =$ {scilabel(dtdnu, dtdnu_err)} ms/MHz'
		)

	ax_wfall.set_title(f"{bname}")
	ax_wfall.legend(loc=1, handlelength=0)

	ax_tseries = axs['T']
	ax_tseries.plot(times_ms-pktime, tseries)
	ax_tseries.add_patch(Rectangle(
		(-window, ax_tseries.get_ylim()[0]),
		width=2*window,
		height=np.max(tseries)*0.075,
		color='tomato',
		alpha=0.0
	))

	for s, xoi in zip(tmix_sigmas, xos):
		w = 2*np.abs(s)
		ax_tseries.add_patch(Rectangle(
			(xoi-pktime-w, ax_tseries.get_ylim()[0]),
			width=2*w,
			height=np.max(tseries)*0.075,
			color='tomato',
			alpha=0.5
		))

	ax_tseries.plot(
		times_ms-pktime,
		gauss_model(
			times_ms-pktime,
			np.max(tseries)*t_popt[0],
			t_popt[1]-pktime,
			t_popt[2]
		),
		'k--',
		alpha=0.1
	)

	# 4 Gaussian model
	tmix_popt[:n_bursts] = [a*np.max(tseries) for a in tmix_amps]
	tmix_popt[n_bursts:n_bursts*2] = [x-pktime for x in tmix_xos]
	ax_tseries.plot(
		times_ms-pktime,
		gaussmix_model(
			times_ms-pktime,
			*tmix_popt
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

	### Slope measurement plot. Plot last component
	arrtimesdf = subdfs[-1]
	ax_slope = axs['E']
	ax_slope.scatter(arrtimesdf['freqs'], arrtimesdf[tpoint]-pktime, c='k', s=20)
	ax_slope.errorbar(
		arrtimesdf['freqs'],
		arrtimesdf[tpoint]-pktime,
		yerr=arrtimesdf['tstart_err'],
		xerr=None,
		fmt='none',
		zorder=-1,
		color='#888888'
	)
	ax_slope.plot(freqs, dtdnu*freqs, 'k--')
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
	ax_slope.set_xlim(min(freqs), max(freqs))
	plt.tight_layout()
	plt.savefig(f"{bname}.png")
	plt.show()
	return results

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

	results = []
	for filename in files:
		burst_results = measureburst(filename, xos=[9.26, 15.75, 22.61, 27.94]) # ms, burst B10 of Sheikh2023
		for row in burst_results:
			results.append(row)

	resultsdf = pd.DataFrame(data=results, columns=[
		'name',
		'DM',
		'center_f',
		'center_f_err',
		'duration',
		'duration_err',
		'bandwidth',
		'bandwidth_err',
		'dtdnu',
		'dtdnu_err'
	]).set_index('name')
	resultsdf.to_csv('results.csv')

