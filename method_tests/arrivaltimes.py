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
		a1, a2, a3, xo1, xo2, xo3, sigma1, sigma2, sigma3 = p
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
		( # this is a lot
			a1, a2, a3, a4, a5,
			xo1, xo2, xo3, xo4, xo5,
			sigma1, sigma2, sigma3, sigma4, sigma5
		) = p
		return (
			gauss_model(x, a1, xo1, sigma1) +
			gauss_model(x, a2, xo2, sigma2) +
			gauss_model(x, a3, xo3, sigma3) +
			gauss_model(x, a4, xo4, sigma4) +
			gauss_model(x, a5, xo5, sigma5)
		)

	if n == 6:
		(
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
	try:
		popt, pcov = scipy.optimize.curve_fit(
			gauss_model,
			x,
			data,
			p0=[
				np.max(data),
				xo,
				np.sqrt(abs(sum(data*(x-xo)**2)/sum(data))) # sigma
			],
		)
	except RuntimeError as e:
		popt = [np.nan, np.nan, np.nan]
		pcov = [np.nan, np.nan, np.nan]
	finally:
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

def fitrows(wfall, dt, freqs, plot=False):
	fitdata = np.zeros((wfall.shape[0], 10))
	for i, row in enumerate(wfall):
		popt, pcov = fitgauss(row, wfall.shape[1]*dt)
		perr = np.sqrt(np.diag(pcov))
		if len(perr.shape) == 2: perr = np.diag(perr) # handles when pcov is nans
		sigma = abs(popt[2])
		tstart = (popt[1]-np.sqrt(2)*sigma)
		tstart_err = np.sqrt(perr[1]**2 + 2*perr[2]**2)
		tend   = (popt[1]+np.sqrt(2)*sigma)
		fitdata[i,:] = [freqs[i], tstart, tend, popt[0], popt[1], tstart_err, sigma, *perr]

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

def measureburst(filename, xos=[], outdir='', show=True):
	""" Measure spectro-temporal properties of a burst, and output a figure

	Compute the inverse sub-burst slope (dt/dnu) using the per-row arrival time method
	Compute the duration and bandwidth by finding a 1-dimensional gaussian model
	to the integrated time series and spectrum, respectively. The duration and bandwidth are the
	1 sigma widths of either fit.
	Compute the center frequency as the center of the 1d spectrum model

	If multiple components are present, split them up and measure individually.

	Args:
		filename (str): filename to npz of burst waterfall
		xos (List[float]): list of times in ms of burst centers. Can be approximate.
		outdir (str): string of output folder for figures
		show (bool): if True show interactive figure window for each file
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

	## Assuming 1 burst:
	window = 2*abs(t_popt[2]) # 2*stddev of a guassian fit to the integrated time series

	##### multi component model: use multiple 1d gaussians to make multiple windows in time,
	# then use the time windows to make frequency windows
	n_bursts = len(xos)

	tmix_popt, tmix_pcov = fitgaussmix(tseries, duration, xos=xos)
	tmix_perr = np.sqrt(np.diag(tmix_pcov))
	tmix_amps   = tmix_popt[:n_bursts]
	tmix_xos    = tmix_popt[n_bursts:n_bursts*2]
	tmix_sigmas = tmix_popt[n_bursts*2:n_bursts*3]
	tmix_sigma_errs = tmix_perr[n_bursts*2:n_bursts*3]

	xos = tmix_xos.tolist() # align to fit component centers

	subfalls = []
	subbands = []
	bandpass = np.zeros(wfall.shape[0])
	xos_chans = np.floor(np.array(xos)/res_time_ms)
	for xoi, s in zip(xos_chans, tmix_sigmas):
		s4 = np.floor(4*np.abs(s)/res_time_ms)
		s1 = np.floor(1*np.abs(s)/res_time_ms)

		# account for when the edge is outside of wfall
		if xoi-s4 < 0:
			subfall = wfall[..., :int(xoi+s4)]
		else:
			subfall = wfall[..., int(xoi-s4):int(xoi+s4)] # 4sigma window around burst
		if xoi-s1 < 0:
			subband = wfall[..., :int(xoi+s1)].mean(axis=1)
		else:
			subband = wfall[..., int(xoi-s1):int(xoi+s1)].mean(axis=1) # sum only 1 sigma from burst peak

		bandpass += subband
		subfall = driftrate.subtractbg(subfall, 0, 10) # subtract bg again, left
		subfall = driftrate.subtractbg(subfall, subfall.shape[1]-1-10, None) # right
		subfalls.append(subfall)
		subbands.append(subband)

	#### Fitting
	dtdnus, subdfs = [], []
	subpeaks, subbandmodels = [], []
	for subfall, subband, xosi, sigma, sigma_err in zip(
		subfalls, subbands, xos, tmix_sigmas, tmix_sigma_errs
	):
		sigma = abs(sigma)
		subdf = fitrows(subfall, res_time_ms, freqs) # Fit a 1d gaussian in each row of the waterfall
		subpktime = 4*sigma # since we made a 4 sigma window

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
		subdf = subdf[subdf.tstart_err/subdf.tstart < 10]
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
		subaxs['W'].plot(
			subtimes,
			(1/dtdnu)*(subtimes-subpktime),
			'w--',
			label=f'{dtdnu=:.2e} ms/MHz'
		)
		subaxs['W'].legend()

		subaxs['B'].plot(
			gauss_model(freqs, *subband_popt),
			freqs
		)
		subbandmodels.append(gauss_model(freqs, *subband_popt))
		# plt.show() # show figure of the cutout sub-burst
		plt.close()

		subdf[tpoint] = subdf[tpoint] + (xosi-4*sigma) # transform to full waterfall times

		dtdnus.append((dtdnu, dtdnu_err))
		subpeaks.append(subpktime)
		subdfs.append(subdf)
		# print(f"{dtdnu = } +/- {dtdnu_err = }")

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

	# Gaussian mix model
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

	### Slope measurement plot. Plot last component
	lastdf = subdfs[-1]
	ax_slope = axs['E']
	ax_slope.scatter(lastdf['freqs'], lastdf[tpoint]-xos[-1], c='k', s=20)
	ax_slope.errorbar(
		lastdf['freqs'],
		lastdf[tpoint]-xos[-1],
		yerr=lastdf[f'{tpoint}_err'],
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

	[ax_tseries.axvline(x=t-pktime, ls='--') for t in xos]
	def printtime(event):
		if event.dblclick:
			x, y = event.xdata, event.ydata
			xos.append(x+pktime)
			# print(f"{x+pktime} ms")
			print(f"xos = {[round(xi,2) for xi in xos]}")
			[ax_tseries.axvline(x=t-pktime, ls='--') for t in xos]
			fig.canvas.draw()
			# print("Adding component...")
			# plt.close()
			# return measureburst(
			# 	filename,
			# 	xos=xos,
			# 	outdir='sheikh/',
			# 	show=True
			# )
	cid = fig.canvas.mpl_connect('button_press_event', printtime)

	if show: plt.show()
	if '/' in outdir or outdir == '':
		outname = f"{outdir}{bname}.png"
	else:
		outname = f"{outdir}/{bname}.png"
	plt.savefig(outname)
	print(f"Saved {outname}.")

	plt.close()
	return results

if __name__ == '__main__':
	# files = glob.glob('/Users/mchamma/dev/SurveyFRB20121102A/data/snelders2023/figure_1/*.npz')
	# files = glob.glob('/Users/mchamma/dev/frbdata/FRB20220912A/sheikh2023/npzs/*.npz')
	# [print(f) for f in sorted(files)]
	# exit()

	prefix = '/Users/mchamma/dev/frbdata/FRB20220912A/sheikh2023/npzs/'
	files = {
		'fil_59872_37843_186776977_FRB20220912a_0001_lob_10sec_crop.npz' : [7, 12.2],
		'fil_59873_29353_191532226_FRB20220912a_0001_lob_10sec_crop.npz' : [7.9, 9.6],
		'fil_59875_21343_201590209_FRB20220912a_0001_lob_10sec_crop.npz' : [6.95, 7.41],
		'fil_59880_21831_227987182_FRB20220912a_0001_lob_10sec_crop.npz' : [8.1, 29.6, 31.7, 35.1, 41.7],
		'fil_59880_32784_228655700_FRB20220912a_0001_lob_10sec_crop.npz' : [7.86, 12.5, 13.84],
		'fil_59880_34610_228767150_FRB20220912a_0001_lob_10sec_crop.npz' : [9.34, 10.95, 17.5],
		'fil_59882_22768_238591247_FRB20220912a_0001_lob_10sec_crop.npz' : [],
		'fil_59882_39198_239594055_FRB20220912a_0001_lob_10sec_crop.npz' : [],#[10.45, 7.35],
		'fil_59883_08510_242994445_FRB20220912a_0001_lob_10sec_crop.npz' : [13.22, 19.08],
		'fil_59883_37405_244758056_FRB20220912a_0001_lob_10sec_crop.npz' : [9.26, 15.75, 22.61, 27.94], # ms, burst B10 of Sheikh2023
		'fil_59887_29267_265355102_FRB20220912a_0001_lob_10sec_crop.npz' : [],#[13.89, 22.89, 30.1, 23.93],
		'fil_59887_36569_265800781_FRB20220912a_0001_lob_10sec_crop.npz' : [],
		'fil_59889_10344_274747009_FRB20220912a_0001_lob_10sec_crop.npz' : [],
		'fil_59889_32251_276084106_FRB20220912a_0001_lob_10sec_crop.npz' : [],
		'fil_59890_26883_281029907_FRB20220912a_0001_lob_10sec_crop.npz' : [],
		'fil_59890_30534_281252746_FRB20220912a_0001_lob_10sec_crop.npz' : [],
		'fil_59891_25899_286243286_FRB20220912a_0001_lob_10sec_crop_1443.npz' : [],
		'fil_59891_25899_286243286_FRB20220912a_0001_lob_10sec_crop_967.npz' : [],
		'fil_59897_05894_316662902_FRB20220912a_0001_lob_10sec_crop.npz' : [],
		'fil_59898_37955_323893188_FRB20220912a_0001_lob_10sec_crop.npz' : [],
		'fil_59900_11856_332847106_FRB20220912a_0001_lob_10sec_crop.npz' : [],
		'fil_59900_22810_333515686_FRB20220912a_0001_lob_10sec_crop.npz' : [],
		'fil_59901_37464_339683532_FRB20220912a_0001_lob_10sec_crop.npz' : [],
		'fil_59903_09465_348521484_FRB20220912a_0001_lob_10sec_crop.npz' : [],
		'fil_59903_20417_349189941_FRB20220912a_0001_lob_10sec_crop.npz' : [],
		'fil_59904_22142_354568664_FRB20220912a_0001_lob_10sec_crop.npz' : [],
		'fil_59907_21795_370367797_FRB20220912a_0001_lob_10sec_crop.npz' : [],
		'fil_59912_25267_396946899_FRB20220912a_0001_lob_10sec_crop.npz' : [],
		'fil_59914_20362_407194396_FRB20220912a_0001_lob_10sec_crop.npz' : [],
		'fil_59915_24158_412699523_FRB20220912a_0001_lob_10sec_crop.npz' : [],
		'fil_59916_25852_418076354_FRB20220912a_0001_lob_10sec_crop.npz' : [],
		'fil_59920_19510_438783020_FRB20220912a_0001_lob_10sec_crop.npz' : [],
		'fil_59930_10525_490968994_FRB20220912a_0001_lob_10sec_crop.npz' : [],
		'fil_59931_23281_497020996_FRB20220912a_0001_lob_10sec_crop.npz' : [],
		'fil_59933_14701_507044189_FRB20220912a_0001_lob_10sec_crop.npz' : []
	}

	results = []
	for filename, xos in files.items():
		filename = f'{prefix}{filename}'
		burst_results = measureburst(
			filename,
			xos=xos,
			outdir='sheikh/',
			show=False
		)
		for row in burst_results:
			results.append(row)

	resultsdf = pd.DataFrame(data=results, columns=[
		'name',
		'DM',
		'center_f (MHz)',
		'center_f_err',
		'duration (ms)',
		'duration_err',
		'bandwidth (MHz)',
		'bandwidth_err',
		'dtdnu (ms/MHz)',
		'dtdnu_err'
	]).set_index('name')
	resultsdf.to_csv('sheikh/results_Feb_9_2024.csv')

