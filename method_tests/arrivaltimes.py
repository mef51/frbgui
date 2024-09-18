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

def line_model(nu, dtdnu, t_b):
	return dtdnu * nu + t_b

def gauss_model(x, a, xo, sigma):
	return a*np.exp(-(x-xo)**2/(2*(sigma**2)))

def smallestdivisor(n):
	for i in range(2, n):
		if n % i == 0:
			return i

def listnpzs(path):
	""" List all npz files in path """
	files = glob.glob(path+'*.npz')
	[print(f) for f in sorted(files)]
	exit()

# N component model
def gaussmix_model(x, *p):
	n = len(p)//3
	model = 0
	for i in range(0, n): # stops at n-1
		model += gauss_model(x, p[0*n+i], p[1*n+i], p[2*n+i])
	return model

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

def fitgaussmix(data, duration, xos, sigmas=None, fix_xos=False, tol=0.01):
	n = len(xos) # Number of components
	if np.max(data) != 0:
		data = data / np.max(data) # normalize
	x = np.linspace(0, duration, num=len(data))
	if not sigmas:
		sigmas = np.sqrt(abs(sum(data*(x-np.mean(xos))**2)/sum(data)))/4
		guess = [*[np.max(data)]*n, *xos, *[sigmas]*n]
	else:
		guess = [*[np.max(data)]*n, *xos, *sigmas]

	bounds = (-np.inf, np.inf)
	if fix_xos:
		bounds = (  # fix xos
			[*[-np.inf]*n, *[xoi - tol for xoi in xos], *[-np.inf]*n],
			[*[np.inf]*n, *[xoi + tol for xoi in xos], *[np.inf]*n]
		)

	try:
		popt, pcov = scipy.optimize.curve_fit(
			gaussmix_model,
			x,
			data,
			p0=guess,
			bounds=bounds
		)
	except RuntimeError as e:
		popt = [0]*3*n
		pcov = [0]*3*n
	finally:
		return popt, pcov

def fitrows(wfall, dt, freqs, plot=False):
	fitdata = np.zeros((wfall.shape[0], 10))
	for i, row in enumerate(wfall):
		popt, pcov = fitgauss(row, wfall.shape[1]*dt)
		# print(f'row {i}: {popt = } {np.mean(row) = } {freqs[i] = }')
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

logdebug = False
def printd(*args):
	if logdebug: print(*args)

results_columns = [
	'name',
	'DM',
	't0 (ms)',
	't0_err',
	'center_f (MHz)',
	'center_f_err',
	'duration (ms)',
	'duration_err',
	'bandwidth (MHz)',
	'bandwidth_err',
	'dtdnu (ms/MHz)',
	'dtdnu_err',
	'tb (ms)',
	'tb_err'
]
def measureburst(
	filename,
	xos=[],
	cuts=[],
	sigmas=None,
	fix_xos=False,
	tolms=0.01,
	targetDM=None,
	correctTimes=False,
	downfactors=(1,1),
	subtractbg=False,
	bw_filter='data_cutoff',
	bw_width_factor=3, # burst based filter factors are likely to help a lot with complex and blended components, and avoids having to do submasks
	snr_cutoff=3,
	t_filter_factor=2,
	crop=None,
	masks=[],
	submasks=None,
	measure_drift=True,
	show=True,
	show_components=False,
	cmap_norm='linear',
	cmap='viridis',
	save=True,
	outdir='',
	outfmt='.png',
	return_arrivaltimes=False,
	return_fig=False,
	loadonly=False,
	save_solutions=False,
	load_solutions=None,
	hide_legend=False,
	legendloc=1,
	label_components=False,
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
		filename (str): filename to npz of a *dedispersed* burst waterfall. Conforms to frbgui's burst format.
		xos (List[float] or 2-tuple of List[float], optional): List of times in ms of sub-burst centers.
			Can be approximate. If a 2-tuple, the second list is used as the location(s) to cut the waterfall.
			Using the ``cuts`` option instead is equivalent to using the 2 tuple option.
		cuts (List[float], optional): List of times in ms to cut the waterfall. Useful for blended bursts.
			User must make sure their cuts make sense (i.e. in between burst centers).
			Typically, if one cut is needed, then all bursts should be cut as well even if well separated.
		sigmas (List[float], optional): initial guesses for the width sigma when finding the 1-dimensional gaussian model
			to the time series. Must be the same length as ``xos``
		fix_xos (bool, optional): Default False. Whether or not to fix the passed ``xos`` when fitting the 1d model.
			Useful when bursts are blended and one can visually distinguish where a burst should be from the waterfall
			even if it appears completely absorbed in the integrated time series.
		tolms (float, optional): Tolerance in milliseconds to use when ``fix_xos`` is True. Default is 0.01 ms.
		targetDM (float, optional): the DM (pc/cm^3) to perform the measurement at.
			Default is to perform the measurement at the DM of the npz file.
		correctTimes (bool, optional): Shift xos and cuts to account for dispersive shift
			when applying a targetDM other than the burst DM. Note that this shift will occur even when ``fix_xos`` is True.
		downfactors (tuple[int], optional): 2-tuple of factors to downsample by in frequency and time (respectively)
		subtractbg (bool, tuple[bool], optional): Perform a second background subtraction on subbursts.
			By default will do a background subtraction using 10% of channels on the whole waterfall.
			Pass `(False, False)` to skip both rounds of background subtraction.
		bw_filter (str, optional): The type of spectral/bandwidth filter to apply on arrival times. Default is `'model_width'`. Options are
			1. `'data_cutoff'`: filter out arrival times in channels where the 1 \\(\\sigma\\) on-pulse mean amplitude is < 10 times the noise amplitude
			2. `'model_cutoff'`: filter out arrival times in channels where the 1d spectral model amplitude is < 10 times the noise amplitude
			3. `'model_width'`: filter out arrival times that lie beyond a multiple of the 1d spectral model width \\(\\sigma\\). See ``bw_width_factor``.
		bw_width_factor (int, optional): By default 3 \\(\\sigma\\) of the burst bandwidth is applied as a spectral filter.
			For bursts with lots of frequency structure this may be inadequate,
			and this parameter can be used to override the filter width. It's recommended to try downsampling first. Note that a
			high `bw_width_factor` such as 10-15 likely indicates the bandwidth measurement is being understimated.
		snr_cutoff (int, optional): The S/N cutoff to use when `bw_filter='data_cutoff'` or `bw_filter='model_cutoff'`.
			By default equals 3.
		t_filter_factor (int, optional): By default 2 \\(\\sigma\\) of the burst duration is applied as a temporal filter.
		outdir (str, optional): string of output folder for figures. Defaults to ''.
		crop (tuple[int], optional): pair of indices to crop the waterfall in time
		masks (List[int], optional): frequency indices to mask. Masks are applied before downsampling
		submasks (tuple[List[int]], optional): tuple of length `xos` of lists of indices to mask on a subcomponent's waterfall.
			Note that contrary to `masks`, these are applied after downsampling.
			Indices are scaled from the original size to the downsampled size and so can cover more than one channel.
			The length of `submasks` must match the length of `xos`.
			Example: To specify a mask on the 4th component of a waterfall with 4 components, pass
			``submask=([],[],[],[22])``.
			This is also a good way to filter out misbehaving components in an otherwise well-measured waterfall
			and is useful for complicated bursts.
		measure_drift (bool, optional): When True (default), and if `len(xos) > 1` (i.e. there are multiple burst components), will measure
			the drift rate using the times and center frequencies of the bursts to fit a line. Will also plot a corresponding
			line showing the drift rate measurement. Set to False to disable this behaviour.
		show (bool, optional): if True show interactive figure window for each file
		show_components (bool, optional): if True show figure window for each sub-burst
		cmap_norm (str, optional) The colormap normalization `norm` parameter passed to matplotlib's imshow command
			when plotting the waterfall. Default is 'linear', other options are 'log', 'symlog', 'logit',
			or matplotlib's Normalize class.
		cmap (str, optional): matplotlib colormap to use for waterfall
		return_arrivaltimes (bool, optional): If True, will a dataframe of the arrival times per channel
		return_fig (bool, optional): if True, return the matplotlib figure. The figure will not be closed.
		save (bool, optional): if True save a figure displaying the measurements.
		loadonly (bool, optional): if True will perform loading steps such as masking, dedispersing,
			and downsampling, then return a tuple of wfall, freqs, times_ms, t_popt, DM, etc.
		outfmt (str, optional): string of file format to save figure as. Default is '.png'. Include the '.' character.
		save_solutions (bool, optional): setting to True will save a file inside of the folder specified by outdir
			that contains the fit solution data for the 1d time series and the spectrum of each component. Useful for reviewing
			measurements of bursts with many components that take a long time to analyse. Output filename will be of the form
			`f'{bname}.sols.npz'. Default False
		load_solutions (str, optional): Filename of solutions file generated by `save_solutions` option. Default is None
		hide_legend (bool, optional): Hides the legend in the output if True.
		legendloc (int or str, optional): Set the location of the legend. Passed to matplotlib's loc argument when the legend is called.
		label_components (bool, optional): If True, label components filters in the time series plot. Useful for complicated waterfalls.


	Returns:
		results (list): list of lists where each list is the result of the measurement.
			This array can be used to make a pandas dataframe in the following way:
			```
			resultsdf = pd.DataFrame(
				data=results,
				columns=arrivaltimes.results_columns
			).set_index('name')
			```

			where the columns of the dataframe are
			```
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
			'tb (ms)', # t_b
			'tb_err'
			```

		arrtimesdf (pd.DataFrame): Only returned when `return_arrivaltimes` is True.
		fig (matplotlib.fig): Matplotlib figure. Only returned when `return_fig` is True.
	"""
	if type(xos) == tuple:
		if len(xos) != 2:
			raise "Error: xos must be a list or tuple of two lists"
		cuts = xos[1]
		xos = xos[0]
	if type(xos) != list:
		raise "Error: xos must be a list"

	xos = sorted(xos)
	cuts = sorted(cuts)

	presubtractbg = True
	if type(subtractbg) == tuple:
		presubtractbg = subtractbg[0]
		subtractbg = subtractbg[1]

	results = []
	bname = filename.split('/')[-1].split('.npz')[0]
	data = np.load(filename, allow_pickle=True)
	wfall = np.copy(data['wfall'])

	if targetDM:
		print(f"Info: Dedispersing from {data['DM']} to {targetDM} pc/cm3")
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
		print(f"Info: {bname}: cropped to {wfall.shape = }")
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
	print(f"Info: {res_freq = :.3f} MHz {res_time_ms = :.5f} ms {min(data['dfs']) = } -- {max(data['dfs']) = } MHz")

	if targetDM and correctTimes:
		ddm = targetDM - data['DM']

		a_dm = 4.14937759336e6
		center_i, errorsumi = driftrate.findCenter(wfall)
		center_f = center_i*res_freq + freqs_bin0
		high_ref_freq = max(data['dfs'])

		deltat = - a_dm * (center_f**-2 - high_ref_freq**-2) * ddm
		xos = [x+deltat for x in xos]
		cuts = [c+deltat for c in cuts]
		print(f'Info: shifting xos and cuts by {deltat} ms for {targetDM = } pc/cm3')

	freqs = np.linspace(freqs_bin0, max(data['dfs']), num=wfall.shape[0]) # channel width/2 is already added
	times_ms = np.linspace(0, duration, num=wfall.shape[1]) # array of timestamps
	tseries = np.nanmean(wfall, axis=0)

	tpoint = 'tstart' # 'tend', 'xo'
	pktime = np.nanargmax(tseries)*res_time_ms
	t_popt, _ = fitgauss(tseries, duration) # whether one or many components, for ref in plot
	print(f"Info: {bname}: {data['wfall'].shape = }, {wfall.shape = }.")
	print(f"Info: Using {bw_filter = } and {snr_cutoff = }")

	if loadonly:
		return (
			wfall,
			freqs,
			times_ms,
			res_freq,
			res_time_ms,
			targetDM,
			t_popt
		)

	if len(xos) == 0:
		xos.append(pktime)

	## Assuming 1 burst:
	# window = t_filter_factor*abs(t_popt[2]) # 2*stddev of a guassian fit to the integrated time series

	##### multi component model: use multiple 1d gaussians to make multiple windows in time,
	# then use the time windows to make frequency windows
	n_bursts = len(xos)

	if load_solutions:
		solsdata = np.load(load_solutions, allow_pickle=True)
		tmix_popt, tmix_perr = solsdata['tmix_popt'], solsdata['tmix_perr']
		subbandpopts, subbandperrs = list(solsdata['subbandpopts']), list(solsdata['subbandperrs'])
		subbandmodels = []
	else:
		tmix_popt, tmix_pcov = fitgaussmix(
			tseries,
			duration,
			xos=xos,
			sigmas=sigmas,
			fix_xos=fix_xos,
			tol=tolms
		)
		tmix_perr = np.sqrt(np.diag(tmix_pcov))
		subbandpopts, subbandmodels, subbandperrs = [], [], []

	if len(tmix_perr.shape) == 2: tmix_perr = np.diag(tmix_perr) # handles when pcov is nans
	tmix_amps       = tmix_popt[:n_bursts]
	tmix_xos        = tmix_popt[n_bursts:n_bursts*2]
	tmix_xos_errs   = tmix_perr[n_bursts:n_bursts*2]
	tmix_sigmas     = tmix_popt[n_bursts*2:n_bursts*3]
	tmix_sigma_errs = tmix_perr[n_bursts*2:n_bursts*3]
	printd(f"'sigmas': {[f for f in tmix_sigmas]}")

	xos = tmix_xos if type(tmix_xos) == list else tmix_xos.tolist() # align to fit component centers
	xos_errs = tmix_xos_errs if type(tmix_xos_errs) == list else tmix_xos_errs.tolist()

	if not submasks:
		submasks = ([],)*len(xos)
	else:
		if len(submasks) != len(xos):
			raise ValueError(f"Please ensure the length of xos and submasks match. {len(xos) = } {len(submasks) = }")

	subfalls = []
	subbands = []
	# sample of noise levels matching
	# number of channels used in corresponding subband integration
	# from beginning of waterfall
	noisesmpls = []
	bandpass = np.zeros(wfall.shape[0])
	xos_chans = np.floor(np.array(xos)/res_time_ms)
	noise_edges = []
	for xoi, s in zip(xos_chans, tmix_sigmas):
		s4 = np.floor(4*np.abs(s)/res_time_ms)
		s1 = np.floor(1*np.abs(s)/res_time_ms)
		if s4 == 0 or s1 == 0:
			s4 = 4 # hack
			s1 = 1 # hack
			lbl = subburst_suffixes[np.where(xos_chans == xoi)[0][0]]
			print(
				f"Warning: Component ({lbl}) has width below the time resolution, possibly due to poor 1D fit. Using 1σ width as 1 channel."
			)

		if len(cuts) == 0:
			# account for when the edge is outside of wfall
			if xoi-s4 < 0:
				subfall = wfall[..., :int(xoi+s4)+1]
			else:
				subfall = wfall[..., int(xoi-s4):int(xoi+s4)+1] # 4sigma window around burst
		else:
			cutchans = np.floor(np.array(cuts)/res_time_ms).astype(int)
			if xoi < cutchans[0]: # left edge
				subfall = wfall[..., :cutchans[0]]
			elif xoi > cutchans[-1]: # right edge
				subfall = wfall[..., cutchans[-1]:]
			else: # middle
				ci = 0
				while xoi > cutchans[ci]: ci += 1
				prev_ci = ci-1

				ci = -1
				while xoi < cutchans[ci]: ci -= 1
				next_ci = ci+1

				subfall = wfall[..., cutchans[prev_ci]:cutchans[next_ci]]

		# Compute component and sample noise
		# Slicing syntax: a[start:stop]  means items start through stop-1
		# Therefore When slicing we add 1 to the end to include the ending channel.
		ddof = 0 # Bessel's correction = 1
		if xoi-s1 < 0: # left edge
			print("Info: Spectral filter noise level sampled from end of waterfall.")
			subband = wfall[..., :int(xoi+s1)+1].mean(axis=1)
			noisesmpls.append(
				wfall[..., -int(xoi+s1):].std(axis=1, ddof=ddof)
			)
			noise_edges.append((len(wfall)-int(xoi+s1), len(wfall)-1))
		elif xoi-s1 > wfall.shape[1]:
			subband = wfall.mean(axis=1) # probably bad fit, take it all as a fallback
			noisesmpls.append(wfall.std(axis=1, ddof=ddof))
			noise_edges.append((0, len(wfall)-1))
			if bw_filter != 'model_width':
				print("Warning: Noise sample taken from whole waterfall. Spectral filter may be overly aggressive.")
		else: # Compute spectrum by summing only 1 sigma from burst peak
			if int(xoi-s1) <= int(xoi+s1)-int(xoi-s1):
				print("Warning: Noise sample overlaps with pulse region.")

			subband = wfall[..., int(xoi-s1):int(xoi+s1)+1].mean(axis=1)
			noisesmpls.append(
				wfall[..., :int(xoi+s1)-int(xoi-s1)+1].std(axis=1, ddof=ddof)
			)
			noise_edges.append((0, int(xoi+s1)-int(xoi-s1)+1))
			printd(
				f"{int(xoi-s1) = }, {int(xoi+s1) = }",
				wfall[..., :int(xoi+s1)-int(xoi-s1)+1].shape,
				wfall[..., int(xoi-s1):int(xoi+s1)+1].shape,
				int(xoi-s1),
				int(xoi+s1)
			)
			if wfall[..., :int(xoi+s1)-int(xoi-s1)+1].shape != wfall[..., int(xoi-s1):int(xoi+s1)+1].shape:
				print("Warning!!!: Subband and noise sample regions differ in size. Check xos.")

		bandpass += subband
		if len(cuts) == 0 and subtractbg: # need to be careful about bg subtraction when cutting
			subfall = driftrate.subtractbg(subfall, 0, 10) # subtract bg again, left
			subfall = driftrate.subtractbg(subfall, subfall.shape[1]-1-10, None) # right
		subfalls.append(subfall)
		subbands.append(subband)

	#### Measurements
	dtdnus, intercepts, subdfs = [], [], []
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

	for subfall, subband, xosi, xosi_err, sigma, sigma_err, submask, noisesmpl in zip(
		subfalls, subbands, xos, xos_errs, tmix_sigmas, tmix_sigma_errs, submasks, noisesmpls
	):
		for m in submask:
			if type(m) == range:
				m = np.array(m)
			subfall[m//downfactors[0]] = 0

		sigma = abs(sigma)
		subdf = fitrows(subfall, res_time_ms, freqs) # Fit a 1d gaussian in each row of the waterfall
		if len(cuts) == 0:
			subpktime = 4*sigma # since we made a 4 sigma window
		else:
			ci, edge = 0, 0
			while ci < len(cuts) and xosi >= cuts[ci]:
				edge = cuts[ci]
				ci += 1
			subpktime = xosi - edge

		# Fit 1d gauss to burst spectrum
		fo = sum(freqs*subband)/sum(subband) # this is an estimate of center_f
		if load_solutions:
			subband_popt, subband_perr = subbandpopts[0], subbandperrs[0]
			subbandpopts = subbandpopts[1:]
			subbandperrs = subbandperrs[1:]
		else:
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
			except (RuntimeError,ValueError) as e:
				print(f"Warning: Spectrum fit failed.")
				subband_popt, subband_perr = [0, 1, 1], [0, 0, 0]

		bwidth, bwidth_err = subband_popt[2], subband_perr[2] # sigma of spetrum fit
		pkfreq, pkfreq_err = subband_popt[1], subband_perr[1] # this is fitted center_f and center_f_err

		## Apply time and spectral filters to points
		printd(f"Debug: pre-filters {len(subdf) = }")
		printd(f"Info: Applying '{bw_filter}' spectral filter ")
		if bw_filter not in ['data_cutoff', 'model_cutoff', 'model_width']:
			print(f"Warning: unrecognized {bw_filter = }. Reverting to 'data_cutoff'")
			bw_filter = 'data_cutoff'
		if bw_filter == 'data_cutoff':
			if logdebug:
				for f,n,s in zip(freqs,noisesmpl,subband):
					printd(f"{f:.3f} MHz: {n = :.8f}\t{s = :.8f}\t{s/n = :.8f}")

			subdf = subdf[ # freqs is the implied axis
				subband/noisesmpl > snr_cutoff
			]
		elif bw_filter == 'model_cutoff' and bwidth != 1: # there must be a fit
			model = gauss_model(subdf['freqs'], *subband_popt)
			subdf = subdf[
				model/noisesmpl > snr_cutoff
			]
		elif bw_filter == 'model_width' and bwidth != 1: # there must be a fit
			subdf = subdf[
				(pkfreq-bw_width_factor*bwidth < subdf['freqs']) &
				(subdf['freqs'] < pkfreq+bw_width_factor*bwidth)
			]

		subdf = subdf[(subdf.amp > 0)]
		subdf = subdf[subdf.tstart_err/subdf.tstart < 10]
		subdf = subdf[
			(subpktime-t_filter_factor*sigma < subdf[tpoint]) &
			(subdf[tpoint] < subpktime+t_filter_factor*sigma) # full witdh
			# (subdf[tpoint] < subpktime) # arrival time must be before pktime
		]
		printd(f"Debug: post-filters {len(subdf) = }")

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
			t_b, tb_err = popt[1], perr[1]

			# print(f"{dtdnu = :.5e} +/- {dtdnu_err:.5e} ms/MHz")
			# print(f"{dtdnu2 = :.5e} +/- {dtdnu_err2:.5e} ms/MHz")
			# print(f"{nu0fit = } +/- {nu0fit_err = }")
			# print(f"{t_b = :.5f} +/- {tb_err:.5f} ms")
		else: # no measurement
			dtdnu, dtdnu_err = 0, 0
			t_b, tb_err = 0, 0

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
				edgecolor='r',
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
		if len(cuts) == 0:
			subdf[tpoint] = subdf[tpoint] + (xosi-4*sigma)
		elif len(cuts) > 0:
			if xosi < cuts[0]: # left edge
				pass # times are already good
			elif xosi > cuts[-1]: # right edge
				subdf[tpoint] = subdf[tpoint] + cuts[-1]
			else: # middle
				ci = 0
				while xosi > cuts[ci]: ci += 1
				prev_ci = ci-1
				subdf[tpoint] = subdf[tpoint] + cuts[prev_ci]

		subdf['color'] = next(colors) # assign color to points

		dtdnus.append((dtdnu, dtdnu_err))
		intercepts.append((t_b, tb_err))
		subbandpopts.append(subband_popt)
		subbandperrs.append(subband_perr)
		subdfs.append(subdf)
		# print(f"{dtdnu = } +/- {dtdnu_err = }")

		rowname = bname if len(xos) == 1 else f'{bname}_{subburst_suffixes[xos.index(xosi)]}'
		results.append([ # see `results_columns`
			rowname,	# 'name',
			targetDM,	# 'DM',
			xosi,		# 't0 (ms)',
			xosi_err,	# 't0_err'
			pkfreq,		# 'center_f (MHz)',
			pkfreq_err,	# 'center_f_err',
			sigma,		# 'duration (ms)',
			sigma_err,	# 'duration_err',
			bwidth,		# 'bandwidth (MHz)',
			bwidth_err,	# 'bandwidth_err',
			dtdnu,		# 'dtdnu (ms/MHz)',
			dtdnu_err,	# 'dtdnu_err',
			t_b,		# 'tb (ms)',
			tb_err		# 'tb_err'
		])

	subdf = pd.concat(subdfs)
	if save_solutions:
		if outdir == '' or outdir[-1] == '/':
			solname = f"{outdir}{bname}-{datetime.now().strftime('%b-%d-%Y')}.sols.npz"
		else:
			solname = f"{outdir}/{bname}-{datetime.now().strftime('%b-%d-%Y')}.sols.npz"
		np.savez(
			solname,
			tmix_popt=tmix_popt,
			tmix_perr=tmix_perr,
			subbandpopts=subbandpopts,
			subbandperrs=subbandperrs
		)
		print(f'Info: Saved {solname} solutions file')

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
		figsize=(10, 9),
		width_ratios=[3,1],
		# gridspec_kw={'hspace':0.464}
	)

	### Waterfall
	ax_wfall = axs['A']
	ax_wfall.imshow(
		wfall,
		aspect='auto',
		origin='lower',
		interpolation='none',
		cmap=cmap,
		extent=extent,
		norm=cmap_norm,
		vmax=np.quantile(wfall, 0.999),
	)
	ax_wfall.annotate(
		f"DM = {targetDM:.3f} pc/cm$^3$",
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

	# Component lines
	for (dtdnu, dtdnu_err), (tb, tb_err), xoi in zip(dtdnus, intercepts, xos):
		if dtdnu != 0:
			ax_wfall.plot(
				times_ms-pktime,
				(1/dtdnu)*(times_ms-xoi) - tb/dtdnu,
				'w-.',
				alpha=0.75,
				# label=f'$dt/d\\nu = $ {dtdnu:.2e} $\\pm$ {dtdnu_err:.2e}'
				label=f'{subburst_suffixes[xos.index(xoi)]}. $dt/d\\nu =$ {scilabel(dtdnu, dtdnu_err)} ms/MHz'
			)

	# Noise sample lines
	# s1 = np.floor(1*np.abs(s)/res_time_ms)
	for ns in noise_edges:
		for n in ns:
			ax_wfall.axvline(
				x=n*res_time_ms-pktime,
				c='r',
				ls='--'
			)

	if len(xos) > 1 and measure_drift:
		targetdf = pd.DataFrame(
			data=results,
			columns=results_columns
		).set_index('name')
		odrjob = scipy.odr.ODR(
			scipy.odr.RealData(
				targetdf['center_f (MHz)'],
				targetdf['t0 (ms)']-pktime,
				sx=targetdf['center_f_err'],
				sy=targetdf['t0_err'],
			),
			scipy.odr.Model(lambda B, x: B[0]*x + B[1]),
			beta0=[-1, 0]
		)
		odrjob.set_job(fit_type=0)
		odrfit = odrjob.run()
		drift, drift_err = odrfit.beta[0], np.sqrt(np.diag(odrfit.cov_beta))[0]
		ax_wfall.plot(
			times_ms-pktime,
			(1/drift)*(times_ms-pktime)+(-odrfit.beta[1]/drift),
			'r-.',
			label=f'Drift: $\\Delta t / \\Delta \\nu = ${scilabel(drift, drift_err)} ms/MHz',
		)
		ax_wfall.errorbar(
			targetdf['t0 (ms)']-pktime,
			targetdf['center_f (MHz)'],
			xerr=targetdf['t0_err'],
			yerr=targetdf['center_f_err'],
			fmt='rX',
			markeredgecolor='k'
		)

	# Test line for 1 parameter line model:
	# ax_wfall.plot(
	# 	times_ms-pktime,
	# 	(1/dtdnu2)*(times_ms-xoi)+pkfreq,
	# 	'y-.',
	# 	alpha=0.75,
	# 	# label=f'$dt/d\\nu = $ {dtdnu2:.2e} $\\pm$ {dtdnu_err:.2e}'
	# 	label=f'{subburst_suffixes[xos.index(xoi)]}. $dt/d\\nu =$ {scilabel(dtdnu2, dtdnu_err2)} ms/MHz'
	# )

	ax_wfall.set_title(f"{bname}")
	if not hide_legend: ax_wfall.legend(loc=legendloc, handlelength=0)

	ax_tseries = axs['T']
	ax_tseries.plot(times_ms-pktime, tseries)

	# plot filter windows (time)
	sp = 0
	for s, xoi in zip(tmix_sigmas, xos):
		w = t_filter_factor*np.abs(s)
		ax_tseries.add_patch(Rectangle(
			(xoi-pktime-w, ax_tseries.get_ylim()[0] + sp*(np.max(tseries)*0.075)),
			width=2*w,
			height=np.max(tseries)*0.075,
			color='tomato',
			alpha=0.5
		))
		if label_components:
			ax_tseries.annotate(
				f"{subburst_suffixes[xos.index(xoi)]}",
				(xoi-pktime, ax_tseries.get_ylim()[0] + sp*(np.max(tseries)*0.075)),
			)
		sp += 1

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
		pkfreq, bw = subbandpopt[1], bw_width_factor*subbandpopt[2]
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
	plt.setp(axs['S'].get_yticklabels(), visible=False) # For paper

	### Slope measurement plot. Plot last component
	plotdf = subdfs[-1]
	ax_slope = axs['E']
	ax_slope.scatter(plotdf['freqs'], plotdf[tpoint]-xos[-1], c='k', s=20)
	ax_slope.errorbar(
		plotdf['freqs'],
		plotdf[tpoint]-xos[-1],
		yerr=plotdf[f'{tpoint}_err'],
		xerr=None,
		fmt='none',
		zorder=-1,
		color='#888888'
	)
	ax_slope.plot(freqs, dtdnu*freqs+t_b, 'k--')
	ax_slope.annotate(
		f"$dt/d\\nu =$ {scilabel(dtdnu, dtdnu_err)} ms/MHz \n$\\sigma_t =$ {scilabel(sigma, sigma_err)} ms",
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

	for s, xoi in zip(tmix_sigmas, xos):
		# 1 sigma regions used for spectrum
		# for xsig in [xoi-pktime-np.abs(s), xoi-pktime+np.abs(s)]:
		# 	ax_tseries.axvline(x=xsig, alpha=1, color='beige', zorder=-20)
		tmplims = ax_tseries.get_ylim() # hack for region not extending all the way
		ax_tseries.fill_between(
			[xoi-pktime-np.abs(s), xoi-pktime+np.abs(s)],
			# 0, 0.1,
			*(10*np.array(ax_tseries.get_ylim())),
			zorder=-20,
			alpha=0.75,
			color='cyan'
		)
		ax_tseries.set_ylim(tmplims)

	def printinfo(event):
		x, y = event.xdata, event.ydata
		if event.dblclick:
			if event.button == 1: # Left click: Select component
				xos.append(x+pktime)
				# print(f"{x+pktime} ms")
				print(f"'xos' : {[round(xi,2) for xi in xos]}")
				[ax_tseries.axvline(x=t-pktime, ls='--') for t in xos]
				fig.canvas.draw()
			elif event.button == 3: # dbl Right click
				cuts.append(x+pktime)
				# print(f"{x+pktime} ms")
				print(f"'cuts' : {[round(xi,2) for xi in cuts]}")
				[ax_tseries.axvline(x=t-pktime, ls='-.', c='k') for t in cuts]
				fig.canvas.draw()
		if event.button == 3: # Right click: select mask channel
			ychan = int(downfactors[0]*(y - freqs_bin0)/res_freq)
			xchan = int(downfactors[1]*(x+pktime)/res_time_ms)
			print(f"freq chan: {ychan} time chan: {xchan}")
			# print(f"freq chan: {ychan} time chan: {xchan} {np.mean(wfall[ychan//downfactors[0]]) = }")
	cid = fig.canvas.mpl_connect('button_press_event', printinfo)

	if show:
		plt.show()
		# Monitor hack
		# mngr = plt.get_current_fig_manager()
		# # to put it into the upper left corner for example:
		# # mngr.window.setGeometry(50,100,640, 545) # newX, newY, dx, dy
		# geom = mngr.window.geometry()
		# x,y,dx,dy = geom.getRect()
		# print(x, y, dx, dy)
	if save:
		if outdir == '' or outdir[-1] == '/':
			outname = f"{outdir}{bname}{outfmt}"
		else:
			outname = f"{outdir}/{bname}{outfmt}"
		fig.savefig(outname)
		print(f"Info: Saved {outname}.")

	if return_arrivaltimes and not return_fig:
		plt.close()
		return results, subdf
	elif not return_arrivaltimes and return_fig:
		return results, fig
	elif return_arrivaltimes and return_fig:
		return results, subdf, fig
	else:
		plt.close()
		return results

drift_columns = [
	'name',
	'DM',
	'center_f (MHz)',
	'center_f_err',
	'duration (ms)',
	'duration_err',
	'drift (ms/MHz)',
	'drift_err',
	'source'
]
def measuredrifts(
	df,
	show_plot=True,
	verbose=False,
	search_paths=[]
):
	"""
	Using a results spreadsheet created from `measureburst`, find bursts with multiple components and compute their drift rate
	using the t0 and center_f measurement data.

	Args:
		df
		show_plot (bool, optional): If True, show a plot of the drift rate measurement. False by default.
		search_paths (list[str], optional): If provided, and if show_plot is True,
			will search for a corresponding .npz file of the burst waterfall and load
			it in order to display a plot of the waterfall with the measured drift rate overlaid.
		verbose (bool, optional): If True, print a lot of information about the drift measurement, including the components used.
	"""
	df = df[df['t0 (ms)'].notna()] # ignore rows that don't have t0 data
	seldf = df.loc[df.index.str.endswith('_a')]
	# bursts with components as identified by the '_a' suffix.
	# If source is a column, then the spreadsheet is multi-source data

	if 'source' in seldf.columns:
		targetbursts = zip(
			[s[:-2] for s in seldf.index],
			seldf['source']
		)
	else:
		targetbursts = [s[:-2] for s in seldf.index]

	if verbose: print(f"Info: {targetbursts = }")
	rows = [] # name, drift (ms/MHz), drift_err, duration (ms), duration_err
	for burst in targetbursts:
		if 'source' in seldf.columns:
			burst, source = burst
			targetdf = df.loc[(df.index.str.startswith(burst+'_')) & (df.source == source)]
		else:
			targetdf = df.loc[df.index.str.startswith(burst+'_')]
			source = ''
		if verbose: print(targetdf)

		duration = targetdf['t0 (ms)'][-1] - targetdf['t0 (ms)'][0] # Δt between first and last components
		duration_err = np.sqrt(targetdf['t0_err'][-1]**2 + targetdf['t0_err'][0]**2)
		mean_freq = targetdf['center_f (MHz)'].mean()
		mean_freq_err = np.sqrt((targetdf['center_f_err']**2).sum())/len(targetdf)
		odrjob = scipy.odr.ODR(
			scipy.odr.RealData(
				targetdf['center_f (MHz)'],
				targetdf['t0 (ms)'],
				sx=targetdf['center_f_err'],
				sy=targetdf['t0_err'],
			),
			scipy.odr.Model(lambda B, x: B[0]*x + B[1]),
			beta0=[-1, 0]
		)
		odrjob.set_job(fit_type=0)
		odrfit = odrjob.run()
		# odrfit.pprint()

		drift, drift_err = odrfit.beta[0], np.sqrt(np.diag(odrfit.cov_beta))[0]
		rows.append([
			burst,
			targetdf['DM'][0],
			mean_freq,
			mean_freq_err,
			duration,
			duration_err,
			drift,
			drift_err,
			source
		])

		if show_plot:
			fig, ax = plt.subplots(1,1, figsize=(6,5))
			times = np.linspace(targetdf['t0 (ms)'].min()*0.9, targetdf['t0 (ms)'].max()*1.1, num=100)
			ax.errorbar(
				targetdf['t0 (ms)'],
				targetdf['center_f (MHz)'],
				xerr=targetdf['t0_err'],
				yerr=targetdf['center_f_err'],
				fmt='ko',
			)
			ax.plot(
				times,
				(1/drift)*times+(-odrfit.beta[1]/drift),
				label=f'ODR: $\\Delta t / \\Delta \\nu = ${scilabel(drift, drift_err)} ms/MHz',
			)
			ax.set_title(f"{burst} $\\sigma_T = $ {duration:.2f} ms")
			plt.legend()
			plt.show()

	return rows

if __name__ == 'n__main__':
	df = pd.read_csv('measurements/frb20121102A/gajjar2018/results_Aug-16-2024-560.105.csv').set_index('name')
	measuredrifts(
		df,
		search_paths=['/Users/mchamma/dev/SurveyFRB20121102A/data/gajjar2018/']
	)


def measure_allmethods(filename, show=True, p0tw=0.01, p0bw=100, **kwargs):
	"""
	Collect spectro-temporal measurements of a burst using multiple techniques.
	Utility for comparing the result of spectro-temporal measurements of a burst obtained from the following
	techniques:
	1. arrival times method
	2. ACF method
	3. direct 2d gaussian with three gaussian forms

	Args:
		filename (str): path to .npz of burst in frbgui burst format
		kwargs: keywords to pass to `measureburst` for loading/processing the waterfall.
			`measureburst` with `loadonly=True` is used to load the burst.
	"""
	model_results = []
	precalc_results = []

	(
		wfall,
		freqs,
		times_ms,
		res_freq,
		res_time_ms,
		targetDM,
		t_popt
	) = measureburst(
		filename,
		loadonly=True,
		**kwargs
	)

	freqs_bin0 = min(freqs)
	tseries = np.nanmean(wfall, axis=0)
	pktime = np.nanargmax(tseries)*res_time_ms
	extent, corrext = driftrate.getExtents(wfall, res_freq, res_time_ms, freqs_bin0)

	fig, axs = plt.subplots(1,2, figsize=(10,5))
	axs[0].imshow(
		wfall,
		aspect='auto',
		origin='lower',
		extent=extent,
		interpolation='none',
		norm='linear',
		vmax=np.quantile(wfall, 0.999),
	)
	axs[0].annotate(
		f"$DM =$ {targetDM:.3f} pc/cm$^3$",
		xy=(0.05, 0.925),
		xycoords='axes fraction',
		color='white',
		weight='black',
		size=10,
		bbox={"boxstyle":"round"}
	)
	axs[0].set_xlabel("Time (ms)")
	axs[0].set_ylabel("Frequency (MHz)")

	## Arrival times measurement
	arr_result, arrtimesdf = measureburst(
		filename,
		save=False,
		show=False,
		return_arrivaltimes=True,
		**kwargs
	)

	arrdf = pd.DataFrame(
		data=arr_result,
		columns=results_columns
	).set_index('name')
	# print(arrdf[['dtdnu (ms/MHz)']], 'ms/MHz', arrdf[['dtdnu_err']], 'ms/MHz')

	print("Arrival Times method:")
	print(f"\t{arrdf.iloc[0]['dtdnu (ms/MHz)']:.4e} +/- {arrdf.iloc[0]['dtdnu_err']:.4e} ms/MHz")

	if len(arrtimesdf) > 0:
		axs[0].scatter( # component fit points
			arrtimesdf['tstart'],
			arrtimesdf['freqs'],
			c='w',
			edgecolors=arrtimesdf['color'],
			marker='o',
			s=25,
			alpha=np.clip(arrtimesdf['amp'], 0, 1),
			label=f"$dt/d\\nu =$ {scilabel(arrdf.iloc[0]['dtdnu (ms/MHz)'], arrdf.iloc[0]['dtdnu_err'])} ms/MHz"
		)

		axs[0].plot( # <--- This line does not appear in the correct spot
			times_ms-pktime,
			(1/arrdf.iloc[0]['dtdnu (ms/MHz)'])*(times_ms-pktime) - arrdf.iloc[0]['tb (ms)']/arrdf.iloc[0]['dtdnu (ms/MHz)'],
			'w--',
			alpha=0.75,
			# label=f'$dt/d\\nu = $ {dtdnu:.2e} $\\pm$ {dtdnu_err:.2e}'
			label='arrival times fit'
		)
	axs[0].set_xlim(extent[0], extent[1])
	axs[0].set_ylim(extent[2], extent[3])


	## ACF measurement (-3470 mhz/ms in frbgui)
	p0s = [
		[1, pktime, freqs[len(freqs)//2], p0tw, p0bw, 0],  # amp, xo, yo, sigma_x, sigma_y, theta
		[1, pktime, 0., freqs[len(freqs)//2], p0tw, p0bw], # amp, t0, dt, nu0, sigma_t, sigma_nu
		[1, pktime, 0., freqs[len(freqs)//2], p0tw, p0bw]  # amp, t0, dnu, nu0, w_t, w_nu
	]

	(
		slope,
		slope_error,
		cpopt,
		cperr,
		_, # theta,
		_, # red_chisq,
		_, # center_f,
		_, # center_f_err,
		fitmap,
	) = driftrate.processBurst(
		wfall,
		res_freq,
		res_time_ms,
		freqs_bin0,
		p0=[],
		verbose=False,
		plot=False
	)
	print(f"ACF method: \n\t{slope = } +/- {slope_error} MHz/ms")
	print(f"\t{1/slope:.4e} +/- {slope_error/slope**2:.4e} ms/MHz")

	precalc_results += [1/slope, slope_error/slope**2]
	model_results += list(cpopt)+list(cperr)

	corr = driftrate.autocorr2d(wfall)

	corrplot = axs[1].imshow(
		corr,
		interpolation='none',
		aspect='auto',
		cmap='gray',
		extent=corrext,
		origin='lower',
	)
	corrplot.set_clim(0, np.max(corr))
	if cpopt[0] > 0:
		c = axs[1].contour(
			fitmap,
			[cpopt[0]/4, cpopt[0]*0.9],
			colors='b',
			alpha=0.75,
			extent=corrext,
		)
		c.collections[0].set_label(
			f"ACF: $dt/d\\nu =$ {scilabel(1/slope, slope_error/slope**2)} ms/MHz"
		)
	axs[1].set_xlabel("Time lag (ms)")
	axs[1].set_ylabel("Frequency lag (MHz)")

	## ACF with constant floor
	print(f"{cpopt = }")
	print("ACF with floor:")
	(
		slope2,
		slope_error2,
		cpopt2,
		cperr2,
		_, # theta,
		_, # red_chisq,
		_, # center_f,
		_, # center_f_err,
		fitmap2,
	) = driftrate.processBurst(
		wfall,
		res_freq,
		res_time_ms,
		freqs_bin0,
		p0=[],
		verbose=False,
		plot=False,
		usefloor=True
	)
	print(f"ACF with floor method: \n\t{slope2 = } +/- {slope_error2} MHz/ms, floor = {cpopt2[6]}")
	print(f"\t{1/slope2:.4e} +/- {slope_error2/slope2**2:.4e} ms/MHz")

	if cpopt2[0] > 0:
		c = axs[1].contour(
			fitmap2,
			[cpopt2[0]/4, cpopt2[0]*0.9],
			colors='g',
			alpha=0.75,
			extent=corrext,
		)
		c.collections[0].set_label(
			f"ACF floored: $dt/d\\nu =$ {scilabel(1/slope2, slope_error2/slope2**2)} ms/MHz \nfloor = {cpopt2[6]:.2e}"
		)
	axs[0].legend()

	## Gaussian models measurements
	models = [driftrate.twoD_Gaussian, driftrate.gaussian_dt, driftrate.gaussian_dnu] # preserve order
	sigma = np.std( wfall[:, 0:50] )
	for model, p0 in zip(models, p0s):
		popt, pcov = driftrate.fitdatagaussiannlsq(
			wfall,
			extent,
			p0=p0,
			sigma=sigma,
			model=model
		)
		perr = np.sqrt(np.diag(pcov))
		# print(popt)
		print(f"{model.__name__}:")

		if model == driftrate.twoD_Gaussian: # amp, xo, yo, sigma_x, sigma_y, theta
			units = ['']*len(p0)
			lbls = ['$A$', '$x_0$', '$y_0$', '$\\sigma_x$', '$\\sigma_y$', '$\\theta$']

			theta = popt[5] if abs(popt[3]) > abs(popt[4]) else popt[5] - np.pi/2 # defined via eqs. A2 and A3
			slope = np.tan(theta) # MHz/ms
			theta_err = perr[-1]
			slope_error = (theta_err * (1/np.cos(theta))**2)

			print(f"\t{slope:.4e} MHz/ms +/- {slope_error}")
			print(f"\t{1/slope:.4e} +/- {slope_error/slope**2:.4e} ms/MHz")
			precalc_results += [1/slope, slope_error/slope**2]
			legend_lbl = f"G: $dt/d\\nu =$ {scilabel(1/slope, slope_error/slope**2)} ms/MHz"
		elif model == driftrate.gaussian_dt: # amp, t0, dt, nu0, sigma_t, sigma_nu
			units = ['', 'ms', 'ms/MHz', 'MHz', 'ms', 'MHz']
			lbls = ['$A_{dt}$', '$t_0$', '$d_t$', '$\\nu_0$', '$\\sigma_t$', '$\\sigma_\\nu$']

			print(f"\t{popt[2]:.4e} +/- {perr[2]:.4e} ms/MHz")
			precalc_results += [popt[2], perr[2]]
			legend_lbl = f"$d_t =$ {scilabel(popt[2], perr[2])} ms/MHz"
		elif model == driftrate.gaussian_dnu: # amp, t0, dnu, nu0, w_t, w_nu
			units = ['', 'ms', 'MHz/ms', 'MHz', 'ms', 'MHz']
			lbls = ['$A_{dnu}$', '$t_0$', '$d_\\nu$', '$\\nu_0$', '$w_t$', '$w_\\nu$']

			gnu_dt = popt[2] / (popt[2]**2 + popt[5]**2/popt[4]**2)

			print(f"\t {gnu_dt = :.4e} ms/MHz")
			precalc_results += [gnu_dt, -1] # uncertainty needs to be derived from eq. A4 in Jahns+2023
			legend_lbl = ''#f"$d_{{t,g\\nu}} =$ {scilabel(gnu_dt, -1)} ms/MHz"
			# if perr[0] == np.inf:
			# 	perr = [-1]*len(perr)

		# Output
		model_results += list(popt)+list(perr)
		resultstr = ', '.join([f"{lbl} = {f:.4e} {unit}" for lbl, f, unit in zip(lbls, popt, units)])
		print('\t', resultstr)
		errstr = ', '.join([f"d{lbl} = {f:.4e} {unit}" for lbl, f, unit in zip(lbls, perr, units)])
		print('\t', errstr)

		poptmap = driftrate.makeDataFitmap(
			popt,
			wfall,
			extent,
			model=model
		)

		if popt[0] > 0:
			c = axs[0].contour(
				poptmap,
				[popt[0]/4, popt[0]*0.9],
				colors='w',
				alpha=0.33,
				extent=extent,
			)
			c.collections[0].set_label(legend_lbl)
		else: # Bad fit, plot in red
			c = axs[0].contour(
				poptmap,
				[-popt[0]/4, -popt[0]*0.9],
				colors='r',
				alpha=0.33,
				extent=extent,
			)
			c.collections[0].set_label(legend_lbl)
		bname = filename.split('/')[-1].split('.')[0]
		axs[0].set_title(bname)

	# axs[0].legend(handlelength=0)
	try:
		axs[0].legend(handlelength=0)
	except IndexError as e:
		print("error: weird legend bug")
	axs[1].legend(ncols=1, handlelength=0)

	results_allmethods = list(arrdf.reset_index().iloc[0]) + precalc_results + model_results

	if show: plt.show()
	plt.savefig(f"measurements/collected/{bname}")
	plt.close()
	return results_allmethods

allmethods_columns = [
	## precalc columns
	'dtdnu_ACF (ms/MHz)',
	'dtdnu_ACF_err',
	'dtdnu_g (ms/MHz)',
	'dtdnu_g_err',
	'dt_gdt (ms/MHz)',
	'dt_gdt_err',
	'dt_gdnu (ms/MHz)',
	'dt_gdnu_err',
	## Raw parameter columns with uncertainties
	'A_ACF',
	'x0_ACF',
	'y0_ACF',
	'sigma_x_ACF',
	'sigma_y_ACF',
	'theta_ACF',
	'A_ACF_err',
	'x0_ACF_err',
	'y0_ACF_err',
	'sigma_x_ACF_err',
	'sigma_y_ACF_err',
	'theta_ACF_err',
	##
	'A_g',
	'x0',
	'y0',
	'sigma_x',
	'sigma_y',
	'theta',
	'A_g_err',
	'x0_err',
	'y0_err',
	'sigma_x_err',
	'sigma_y_err',
	'theta_err',
	##
	'A_dt',
	't0_gdt (ms)',
	'd_t (ms/MHz)',
	'nu0_gdt (MHz)',
	'sigma_t (ms)',
	'sigma_nu (MHz)',
	'A_dt_err',
	't0_gdt_err',
	'd_t_err',
	'nu0_gdt_err',
	'sigma_t_err',
	'sigma_nu_err',
	##
	'A_dnu',
	't0_gdnu (ms)',
	'd_nu (MHz/ms)',
	'nu0_gdnu (MHz)',
	'w_t (ms)',
	'w_nu (MHz)',
	'A_dnu_err',
	't0_gdnu_err',
	'd_nu_err',
	'nu0_gdnu_err',
	'w_t_err',
	'w_nu_err',
]


# if __name__ == '__main__':
# 	files = glob.glob('/Users/mchamma/dev/SurveyFRB20121102A/data/snelders2023/figure_1/*.npz')
# 	files = glob.glob('/Users/mchamma/dev/frbdata/FRB20220912A/sheikh2023/npzs/*.npz')
#	files = glob.glob(prefix+'*.npz')
# 	[print(f) for f in sorted(files)]
# 	exit()

