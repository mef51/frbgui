import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import your
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle, Ellipse
import scipy, glob
from itertools import zip_longest, cycle
from tqdm import tqdm
import pandas as pd
# from sklearn.mixture import GaussianMixture
import driftrate
from driftrate import scilabel, subburst_suffixes

# Based on https://github.com/mef51/subdriftlaw/blob/master/ArrivalTimes.ipynb

def zero_line_model(nu, dtdnu):
	return dtdnu * nu

def line_model(nu, dtdnu, nu0dtaudnu):
	return dtdnu * nu + nu0dtaudnu

def gauss_model(x, a, xo, sigma):
	return a*np.exp(-(x-xo)**2/(2*(sigma**2)))

def smallestdivisor(n):
	for i in range(2, n):
		if n % i == 0:
			return i

###### N component models. Can these be generated?
# See GaussianMixture models (sklearn, lmfit, general) maybe
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
	if n == 7:
		(
			a1, a2, a3, a4, a5, a6, a7,
			xo1, xo2, xo3, xo4, xo5, xo6, xo7,
			sigma1, sigma2, sigma3, sigma4, sigma5, sigma6, sigma7
		) = p
		return (
			gauss_model(x, a1, xo1, sigma1) +
			gauss_model(x, a2, xo2, sigma2) +
			gauss_model(x, a3, xo3, sigma3) +
			gauss_model(x, a4, xo4, sigma4) +
			gauss_model(x, a5, xo5, sigma5) +
			gauss_model(x, a6, xo6, sigma6) +
			gauss_model(x, a7, xo7, sigma7)
		)
	if n == 8:
		(
			a1, a2, a3, a4, a5, a6, a7, a8,
			xo1, xo2, xo3, xo4, xo5, xo6, xo7, xo8,
			sigma1, sigma2, sigma3, sigma4, sigma5, sigma6, sigma7, sigma8
		) = p
		return (
			gauss_model(x, a1, xo1, sigma1) +
			gauss_model(x, a2, xo2, sigma2) +
			gauss_model(x, a3, xo3, sigma3) +
			gauss_model(x, a4, xo4, sigma4) +
			gauss_model(x, a5, xo5, sigma5) +
			gauss_model(x, a6, xo6, sigma6) +
			gauss_model(x, a7, xo7, sigma7) +
			gauss_model(x, a8, xo8, sigma8)
		)
	if n == 9:
		(
			a1, a2, a3, a4, a5, a6, a7, a8, a9,
			xo1, xo2, xo3, xo4, xo5, xo6, xo7, xo8, xo9,
			sigma1, sigma2, sigma3, sigma4, sigma5, sigma6, sigma7, sigma8, sigma9
		) = p
		return (
			gauss_model(x, a1, xo1, sigma1) +
			gauss_model(x, a2, xo2, sigma2) +
			gauss_model(x, a3, xo3, sigma3) +
			gauss_model(x, a4, xo4, sigma4) +
			gauss_model(x, a5, xo5, sigma5) +
			gauss_model(x, a6, xo6, sigma6) +
			gauss_model(x, a7, xo7, sigma7) +
			gauss_model(x, a8, xo8, sigma8) +
			gauss_model(x, a9, xo9, sigma9)
		)

########### End n-component models

def fitgauss(data, duration):
	# use curve-fit (non-linear leastsq)
	if len(data) == 0:
		popt = [np.nan, np.nan, np.nan]
		pcov = [np.nan, np.nan, np.nan]
		return popt, pcov
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

def fitgaussmix(data, duration, xos):
	n = len(xos) # Number of components
	if np.max(data) != 0:
		data = data / np.max(data) # normalize
	x = np.linspace(0, duration, num=len(data))
	sigmas = np.sqrt(abs(sum(data*(x-np.mean(xos))**2)/sum(data)))/4
	guess = [*[np.max(data)]*n, *xos, *[sigmas]*n]
	try:
		popt, pcov = scipy.optimize.curve_fit(gaussmix_model, x, data, p0=guess)
	except RuntimeError as e:
		popt = [0]*3*n
		pcov = [0]*3*n
	finally:
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
		extent=extent,
		norm='linear',
		vmax=np.quantile(data, 0.999),
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

def measureburst(
	filename,
	xos=[],
	targetDM=None,
	downfactors=(1,1),
	subtractbg=False,
	crop=None,
	masks=[],
	outdir='',
	show=True,
	show_components=False,
	save=True
):
	""" Measure spectro-temporal properties of a burst, and output a figure

	Compute the inverse sub-burst slope (dt/dnu) using the per-row arrival time method
	Compute the duration and bandwidth by finding a 1-dimensional gaussian model
	to the integrated time series and spectrum, respectively. The duration and bandwidth are the
	1 sigma widths of either fit.
	Compute the center frequency as the center of the 1d spectrum model

	If multiple components are present, split them up and measure individually. The number
	of components to fit for is equal to ``len(xos)``

	Args:
		filename (str): filename to npz of a *dedispersed* burst waterfall
		xos (List[float] or 2-tuple of List[float], optional): List of times in ms of sub-burst centers.
			Can be approximate. If a 2-tuple, the second list is used as the location(s) to cut the waterfall.
			Only supports one cut at the moment.
		targetDM (float, optional): the DM (pc/cm^3) to perform the measurement at.
			Default is to perform the measurement at the DM of the npz file.
		downfactors (tuple[int], optional): 2-tuple of factors to downsample by in frequency and time (respectively)
		subtractbg (bool, tuple[bool], optional): Perform a second background subtraction on subbursts.
			By default will do a background subtraction using 10% of channels on the whole waterfall.
			Pass `(False, False)` to skip both rounds of background subtraction.
		outdir (str, optional): string of output folder for figures
		crop (tuple[int], optional): pair of indices to crop the waterfall in time
		masks (List[int], optional): frequency indices to mask
		show (bool, optional): if True show interactive figure window for each file
		show_components (bool, optional): if True show figure window for each sub-burst
		save (bool, optional): if True save a figure displaying the measurements.
	"""
	cuts = []
	if type(xos) == tuple:
		if len(xos) != 2:
			raise "Error: xos must be a list or tuple of two lists"
		cuts = xos[1]
		if len(cuts) > 1:
			raise "Error: more than 1 cut is not implemented"
		xos = xos[0]
	if type(xos) != list:
		raise "Error: xos must be a list"

	xos = sorted(xos)

	presubtractbg = True
	if type(subtractbg) == tuple:
		presubtractbg = subtractbg[0]
		subtractbg = subtractbg[1]

	results = []
	bname = filename.split('/')[-1].split('.')[0]
	data = np.load(filename, allow_pickle=True)
	wfall = np.copy(data['wfall'])

	if targetDM:
		ddm = targetDM - data['DM']
		wfall = driftrate.dedisperse(
			wfall,
			ddm,
			min(data['dfs']),
			data['bandwidth']/wfall.shape[0],
			1000*data['duration']/wfall.shape[1]
		)
	else:
		targetDM = data['DM']

	for mask in masks:
		wfall[mask] = 0
	if presubtractbg:
		wfall = driftrate.subtractbg(wfall, 0, int(wfall.shape[1]*0.1))
	if type(crop) == tuple or type(crop) == list and len(crop) == 2:
		wfall = wfall[..., crop[0]:crop[1]]
	wfall = driftrate.subsample(
		wfall,
		wfall.shape[0]//downfactors[0],
		wfall.shape[1]//downfactors[1]
	)

	# determine resolutions accounting for downsampling and cropping
	freqs_bin0  = min(data['dfs'])
	res_freq    = downfactors[0] * data['bandwidth'] / data['wfall'].shape[0] # MHz
	res_time_ms = downfactors[1] * 1000*data['duration'] / data['wfall'].shape[1] # ms
	duration    = wfall.shape[1] * res_time_ms # duration of potentially cropped waterfall

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
	if len(tmix_perr.shape) == 2: tmix_perr = np.diag(tmix_perr) # handles when pcov is nans
	tmix_amps   = tmix_popt[:n_bursts]
	tmix_xos    = tmix_popt[n_bursts:n_bursts*2]
	tmix_sigmas = tmix_popt[n_bursts*2:n_bursts*3]
	tmix_sigma_errs = tmix_perr[n_bursts*2:n_bursts*3]

	xos = tmix_xos if type(tmix_xos) == list else tmix_xos.tolist() # align to fit component centers

	subfalls = []
	subbands = []
	bandpass = np.zeros(wfall.shape[0])
	xos_chans = np.floor(np.array(xos)/res_time_ms)
	for xoi, s in zip(xos_chans, tmix_sigmas):
		s4 = np.floor(4*np.abs(s)/res_time_ms)
		s1 = np.floor(1*np.abs(s)/res_time_ms)
		if s4 == 0 or s1 == 0:
			s4 = 10 # hack
			s1 = 10 # hack
			lbl = subburst_suffixes[np.where(xos_chans == xoi)[0][0]]
			print(
				f"Warning: Component ({lbl}) has 0 width, likely due to poor 1D fit"
			)

		if len(cuts) == 0:
			# account for when the edge is outside of wfall
			if xoi-s4 < 0:
				subfall = wfall[..., :int(xoi+s4)]
			else:
				subfall = wfall[..., int(xoi-s4):int(xoi+s4)] # 4sigma window around burst
		elif len(cuts) == 1:
			cutchan = int(np.floor(np.array(cuts)/res_time_ms)[0])
			if xoi < cutchan:
				subfall = wfall[..., :cutchan]
			elif xoi > cutchan:
				subfall = wfall[..., cutchan:]

		if xoi-s1 < 0:
			subband = wfall[..., :int(xoi+s1)].mean(axis=1)
		elif xoi-s1 > wfall.shape[1]:
			subband = wfall.mean(axis=1) # probably bad fit, take it all as a fallback
		else:
			subband = wfall[..., int(xoi-s1):int(xoi+s1)].mean(axis=1) # sum only 1 sigma from burst peak

		bandpass += subband
		if len(cuts) == 0 and subtractbg: # need to be careful about bg subtraction when cutting
			subfall = driftrate.subtractbg(subfall, 0, 10) # subtract bg again, left
			subfall = driftrate.subtractbg(subfall, subfall.shape[1]-1-10, None) # right
		subfalls.append(subfall)
		subbands.append(subband)

	#### Fitting
	dtdnus, intercepts, subdfs = [], [], []
	subbandpopts, subbandmodels = [], []
	colors = cycle([
		'white',
		'black',
		'red',
		'green',
		'blue',
		'yellow',
		'darkgreen',
		'brown'
	])
	for subfall, subband, xosi, sigma, sigma_err in zip(
		subfalls, subbands, xos, tmix_sigmas, tmix_sigma_errs
	):
		sigma = abs(sigma)
		subdf = fitrows(subfall, res_time_ms, freqs) # Fit a 1d gaussian in each row of the waterfall
		if len(cuts) == 0:
			subpktime = 4*sigma # since we made a 4 sigma window
		else:
			subpktime = np.nanargmax(subfall.mean(axis=0))*res_time_ms

		# Fit 1d gauss to burst spectrum
		fo = sum(freqs*subband)/sum(subband) # this is an estimate of center_f
		try:
			subband_popt, subband_pcov = scipy.optimize.curve_fit(
				gauss_model,
				freqs,
				subband/np.max(subband),
				p0=[
					1,
					fo,
					np.sqrt(abs(sum(subband*(freqs-fo)**2)/sum(subband))) # sigma
				],
			)
			subband_perr = np.sqrt(np.diag(subband_pcov))
		except RuntimeError as e:
			print(f"Warning: Spectrum fit failed. Skipping burst bandwidth filter.")
			subband_popt, subband_perr = [0, 1, 1], [0, 0, 0]

		bwidth, bwidth_err = subband_popt[2], subband_perr[2] # sigma of spetrum fit
		pkfreq, pkfreq_err = subband_popt[1], subband_perr[1] # this is fitted center_f

		## Apply time and spectral filters to points
		subdf = subdf[(subdf.amp > 0)]
		subdf = subdf[subdf.tstart_err/subdf.tstart < 10]
		subdf = subdf[(subpktime-2*sigma < subdf[tpoint]) & (subdf[tpoint] < subpktime+2*sigma)]
		if bwidth != 1:
			subdf = subdf[
				(pkfreq-3*bwidth < subdf['freqs']) &
				(subdf['freqs'] < pkfreq+3*bwidth)
			]

		# Measure dt/dnu
		if len(subdf) > 1: # only fit a line if more than 1 point
			popt, pcov = scipy.optimize.curve_fit(
				line_model,
				subdf['freqs'],
				subdf[tpoint] - subpktime,
				sigma=subdf[f'{tpoint}_err'],
				absolute_sigma=True,
			)
			perr = np.sqrt(np.diag(pcov))
			dtdnu, dtdnu_err = popt[0], perr[0]
			nu0dtaudnu, nu0dtaudnu_err = popt[1], perr[1]

			# print(f"{dtdnu = :.5f} +/- {dtdnu_err:.5f} ms/MHz")
			# print(f"{nu0dtaudnu = :.5f} +/- {nu0dtaudnu_err:.5f} ms")
		else: # no measurement
			dtdnu, dtdnu_err = 0, 0
			nu0dtaudnu, nu0dtaudnu_err = 0, 0


		# Sub-burst plot
		if show_components:
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
			if dtdnu != 0:
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
			plt.show()
			plt.close()
		subbandmodels.append(gauss_model(freqs, *subband_popt))

		# transform times to full waterfall times
		if len(cuts) == 1:
			if xosi < cuts[0]:
				pass # times are already good
			elif xosi > cuts[0]:
				subdf[tpoint] = subdf[tpoint] + cuts[0]
		else:
			subdf[tpoint] = subdf[tpoint] + (xosi-4*sigma)
		subdf['color'] = next(colors) # assign color to points

		dtdnus.append((dtdnu, dtdnu_err))
		intercepts.append((nu0dtaudnu, nu0dtaudnu_err))
		subbandpopts.append(subband_popt)
		subdfs.append(subdf)
		# print(f"{dtdnu = } +/- {dtdnu_err = }")

		rowname = bname if len(xos) == 1 else f'{bname}_{subburst_suffixes[xos.index(xosi)]}'
		results.append([
			rowname,
			targetDM,
			pkfreq,
			pkfreq_err,
			sigma,
			sigma_err,
			bwidth,
			bwidth_err,
			dtdnu,
			dtdnu_err,
			nu0dtaudnu,
			nu0dtaudnu_err
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
	ax_wfall.annotate(
		f"$DM =$ {targetDM:.3f} pc/cm$^3$",
		xy=(0.05, 0.925),
		xycoords='axes fraction',
		color='white',
		weight='black',
		size=10,
		bbox={"boxstyle":"round"}
	)
	if len(subdf) > 0:
		ax_wfall.scatter( # component fit points
			subdf[tpoint]-pktime,
			subdf['freqs'],
			c='w',
			edgecolors=subdf['color'],
			marker='o',
			s=25,
			alpha=np.clip(subdf['amp'], 0, 1)
		)
	ax_wfall.set_xlabel("Time (ms)")
	ax_wfall.set_ylabel("Frequency (MHz)")

	# component lines
	for (dtdnu, dtdnu_err), (tb, tb_err), xoi in zip(dtdnus, intercepts, xos):
		if dtdnu != 0:
			ax_wfall.plot(
				times_ms-pktime,
				(1/dtdnu)*(times_ms-xoi) - tb/dtdnu,
				'w--',
				alpha=0.75,
				# label=f'$dt/d\\nu = $ {dtdnu:.2e} $\\pm$ {dtdnu_err:.2e}'
				label=f'{subburst_suffixes[xos.index(xoi)]}. $dt/d\\nu =$ {scilabel(dtdnu, dtdnu_err)} ms/MHz'
			)

	ax_wfall.set_title(f"{bname}")
	ax_wfall.legend(loc=1, handlelength=0)

	ax_tseries = axs['T']
	ax_tseries.plot(times_ms-pktime, tseries)

	# plot filter windows (time)
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
	downband = 4
	if len(bandpass) % downband != 0:
		downband = smallestdivisor(len(bandpass))
	bandpass_down = bandpass.reshape(-1, downband).mean(axis=1)
	bandpass_down = bandpass_down / np.max(bandpass_down)
	axs['S'].stairs(
		bandpass_down,
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

	# plot filter windows (frequency)
	for i, subbandpopt in enumerate(subbandpopts):
		pkfreq, bw = subbandpopt[1], 3*subbandpopt[2]
		rectwidth = np.max(bandpass_down)*0.035
		axs['S'].add_patch(Rectangle(
			(0.025+axs['S'].get_xlim()[0] + i*(rectwidth+0.025), pkfreq-bw),
			height=2*bw,
			width=rectwidth,
			color='tomato',
			alpha=0.5
		))

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
	ax_slope.plot(freqs, dtdnu*freqs+nu0dtaudnu, 'k--')
	ax_slope.annotate(
		f"$dt/d\\nu =$ {scilabel(dtdnu, dtdnu_err)} ms/MHz \n$t_b =$ {scilabel(nu0dtaudnu, nu0dtaudnu_err)} ms",
		xy=(0.7, 0.6),
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


	for t in xos:
		ax_tseries.axvline(x=t-pktime, ls='--')
	for cut in cuts:
		ax_tseries.axvline(
			x=cut-pktime,
			ls='-.',
			c='k',
			zorder=-1,
		)

	def printtime(event):
		if event.dblclick:
			x, y = event.xdata, event.ydata
			xos.append(x+pktime)
			# print(f"{x+pktime} ms")
			print(f"'xos' : {[round(xi,2) for xi in xos]}")
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
	if save:
		if outdir[-1] == '/' or outdir == '':
			outname = f"{outdir}{bname}.png"
		else:
			outname = f"{outdir}/{bname}.png"
		fig.savefig(outname)
		print(f"Saved {outname}.")

	plt.close()
	return results

def listnpzs(path):
	""" List all npz files in path """
	files = glob.glob(path+'*.npz')
	[print(f) for f in sorted(files)]
	exit()

results_columns = [
	'name',
	'DM',
	'center_f (MHz)',
	'center_f_err',
	'duration (ms)',
	'duration_err',
	'bandwidth (MHz)',
	'bandwidth_err',
	'dtdnu (ms/MHz)',
	'dtdnu_err',
	'tb (ms)', # nu0dtaudnu
	'tb_err'
]

# if __name__ == '__main__':
# 	files = glob.glob('/Users/mchamma/dev/SurveyFRB20121102A/data/snelders2023/figure_1/*.npz')
# 	files = glob.glob('/Users/mchamma/dev/frbdata/FRB20220912A/sheikh2023/npzs/*.npz')
#	files = glob.glob(prefix+'*.npz')
# 	[print(f) for f in sorted(files)]
# 	exit()



