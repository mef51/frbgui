import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import numpy as np
import scipy.optimize
import itertools, glob
import pandas as pd
import os
from tqdm import tqdm

def findCenter(burstwindow):
	freqspectrum = burstwindow.sum(axis=1)[:, None]
	freqi = np.indices(freqspectrum.shape)[0]
	return np.nansum(freqi*freqspectrum) / np.nansum(freqspectrum)

def structureParameter(wfall, dt, tstart, tend):
	"""
	wip. see eq. 1 in gajjar et al. 2018
	dt     - time resolution
	tstart - chan #
	tend   - chan #
	"""
	n = (tend - tstart)
	ts = np.nanmean(wfall, axis=0)
	struct = 0
	for i in enumerate(ts[tstart:tend]):
		struct += abs((ts[i] - ts[i+1]) / dt)

	return struct/n

def subband(wfall, nsub):
	nchan, nsamp = wfall.shape
	sub_factor = nchan // nsub
	return np.nanmean(wfall.reshape(-1, sub_factor, nsamp), axis=1)

def subsample(m, nfreq, ntime):
	""" m : 2x2 array """
	n = np.nanmean(m.reshape(-1, m.shape[0]//nfreq, m.shape[1]), axis=1)
	return np.nanmean(n.reshape(n.shape[0], -1, n.shape[1]//ntime), axis=2)

def subtractbg(wfall, tleft: int=0, tright: int=1):
	return wfall - wfall[:, tleft:tright].mean(axis=1)[:, None]

def moments(data):
	"""Returns (height, x, y, width_x, width_y)
	the gaussian parameters of a 2D distribution by calculating its
	moments """
	total = data.sum()
	X, Y = np.indices(data.shape) # this seems inconsistent with np.meshgrid
	x = (X*data).sum()/total
	y = (Y*data).sum()/total
	col = data[:, int(y)]
	width_x = np.sqrt(abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
	row = data[int(x), :]
	width_y = np.sqrt(abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
	height = data.max()
	# print('using default p0:', height, data.shape[1]/2, data.shape[0]/2, width_x, width_y, 2.0)
	return height, data.shape[1]/2, data.shape[0]/2, width_x, width_y, 2.0
	# return height, x, y, width_x, width_y, 2.0

def twoD_Gaussian(point, amplitude, xo, yo, sigma_x, sigma_y, theta):
	y, x = point
	xo = float(xo)
	yo = float(yo)
	a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
	b = (np.sin(2*theta))/(2*sigma_x**2) - (np.sin(2*theta))/(2*sigma_y**2)
	c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
	g = amplitude*np.exp( - a*((x-xo)**2) - b*(x-xo)*(y-yo) - c*((y-yo)**2))
	return g.ravel()

def fitgaussiannlsq(data, p0=[], sigma=0, bounds=(-np.inf, np.inf)):
	# use curve-fit (non-linear leastsq)
	x = range(0, data.shape[1]); y = range(0, data.shape[0])
	x, y = np.meshgrid(x, y)
	p0 = moments(data) if p0 == [] else p0
	sigma = np.zeros(len(data.ravel())) + sigma
	popt, pcov = scipy.optimize.curve_fit(twoD_Gaussian, (y, x), data.ravel(), p0=p0, sigma=sigma,
										  absolute_sigma=True, bounds=bounds)
	return popt, pcov

def getDataCoords(extents, shape):
	x = np.linspace(extents[0], extents[1], num=shape[1])
	y = np.linspace(extents[2], extents[3], num=shape[0])
	x, y = np.meshgrid(x, y)
	return x, y

def fitdatagaussiannlsq(data, extents, p0=[], sigma=0, bounds=(-np.inf, np.inf)):
	# use curve-fit (non-linear leastsq)
	# x = range(0, data.shape[1]); y = range(0, data.shape[0])
	x, y = getDataCoords(extents, data.shape)
	# p0 = moments(data) if p0 == [] else p0 # moments is in channel space, needs to be rewritten for data space
	p0 = [1,1,1,1,1,1] if p0 == [] else p0
	sigma = np.zeros(len(data.ravel())) + sigma
	popt, pcov = scipy.optimize.curve_fit(twoD_Gaussian, (y, x), data.ravel(), p0=p0, sigma=sigma,
										  absolute_sigma=True, bounds=bounds)
	return popt, pcov

def makeDataFitmap(popt, corr, extents):
	x, y = getDataCoords(extents, corr.shape)
	fitmap = twoD_Gaussian((y, x), *popt).reshape(corr.shape[0], corr.shape[1])
	return fitmap

def _dedisperse(wfall, dm, freq, dt):
	"""Dedisperse a dynamic spectrum.

	Parameters
	----------
	wfall : array_like
		Dynamic spectra of shape (nchan, nsamp).
	dm : float
		Dispersion measure to dedisperse to, in pc cm-3.
	freq : array_like
		Center frequencies of all channels, in MHz. Should have shape nchan.
	dt : float
		Sampling time, in s.

	Returns
	-------
	wfall : array_like
		Dedispersed dynamic spectra of shape (nchan, nsamp).

	"""
	k_dm = 1. / 2.41e-4
	dedisp = np.zeros_like(wfall)

	ref_freq = freq[0]### ORIGINAL
	# ref_freq = freq[-1]
	# print("ref_freq", ref_freq)

	shift = (k_dm * dm * (ref_freq ** -2 - freq ** -2) / dt) ### ORIGINAL (low freq anchor)
	# shift = (k_dm * dm * (freq ** -2 - ref_freq ** -2) / dt)
	shift = shift.round().astype(int)

	for i, ts in enumerate(wfall):
		dedisp[i] = np.roll(ts, shift[i])

	return dedisp

def dedisperse(intensity, DM, nu_low, df_mhz, dt_ms, cshift=0):
	dedispersed = np.copy(intensity)

	shifts = [0 for i in range(0, len(intensity))]
	high_ref_freq = nu_low + len(dedispersed)*df_mhz
	low_ref_freq  = nu_low
	#k_dm = 4.1488064239e6 # kulkarni
	k_dm = 4.14937759336e6 # pulsar community
	for i, row in enumerate(dedispersed): # i == 0 corresponds to bottom of the band
		nu_i = nu_low + i*df_mhz
		# High frequency anchor
		deltat = - k_dm * (nu_i**-2 - high_ref_freq**-2) * DM

		# Low frequency anchor
		#deltat = 4.14937759336e6 * (low_ref_freq**-2 - nu_i**-2) * DM

		channelshift = int(round(deltat/dt_ms))
		dedispersed[i] = np.roll(dedispersed[i], channelshift)

	# optionally center view
	dedispersed = np.roll(dedispersed, cshift, axis=1)

	return dedispersed

def getExtents(wfall, df:float=1.0, dt:float=1.0, lowest_freq:float=1.0):
	extents = (0,
			   dt*wfall.shape[1],
			   lowest_freq,
			   lowest_freq + df*wfall.shape[0])

	corrextents = (-extents[1], extents[1], -(extents[3]-extents[2]), (extents[3]-extents[2]))
	return extents, corrextents

def cropwfall(wfall, twidth=150, pkidx=None):
	twidth = round(twidth)
	if twidth <= 0:
		twidth = wfall.shape[1]
	wfall = wfall.copy()
	ts    = np.nanmean(wfall, axis=0)
	if not pkidx:
		pkidx = np.nanargmax(ts)

	ledge, redge = pkidx-twidth, pkidx+twidth
	if ledge < 0:
		ledge = 0
	if redge == 0: # slicing w/ None takes the whole array
		redge = None

	return wfall[..., ledge:redge]

def autocorr2d(data):

	# Returns a 2D autocorrelation computed via an intermediate FFT

	# Number of data pts
	nx, ny = data.shape[0], data.shape[1]

	padded = np.append(data, np.zeros((nx,ny)), axis=0)
	padded = np.append(padded, np.zeros((2*nx,ny)), axis=1)

	# Perform the FFT
	data_dft = np.fft.fft2(padded)

	# DFT of auto-correlation is simply (conjugate) multiplication
	# Elt-wise multiplication of fft
	data_ac_dft = np.multiply(np.conjugate(data_dft), data_dft)

	# Inverse FFT to return to time
	# Note this array will be half-shifted
	result_shifted = np.fft.ifft2(data_ac_dft)

	# Flip the result array around
	temp_array_a = np.empty(result_shifted.shape)
	temp_array_b = np.empty(result_shifted.shape)

	# Flip in x:
	temp_array_a[0:nx,:] = result_shifted[nx-1:2*nx-1,:]
	temp_array_a[nx:2*nx-1,:] = result_shifted[0:nx-1,:]
	# Flip in y:
	temp_array_b[:,0:ny] = temp_array_a[:,ny-1:2*ny-1]
	temp_array_b[:,ny:2*ny-1] = temp_array_a[:,0:ny-1]

	return temp_array_b[:-1,:-1]#/float(nx*ny)

def processBurst(burstwindow, fres_MHz, tres_ms, lowest_freq, burstkey=1, p0=[], popt_custom=[],
				 bounds=(-np.inf, np.inf), nclip=None, clip=None, plot=False,
				 sigmawindow=(0,50),
				 verbose=True):
	"""
	Given a waterfall of a burst, will use the 2d autocorrelation+gaussian fitting method to
	find the drift and make a plot of the burst and fit.
	returns drift, drift_error, popt, perr, theta,	red_chisq, center_f
	"""

	corr = autocorr2d(burstwindow)
	_, corrextents = getExtents(burstwindow,
									  df=fres_MHz, dt=tres_ms, lowest_freq=lowest_freq)
	# print('wfall info:', f'{np.max(burstwindow) = }, {np.mean(burstwindow) = }, {burstwindow.shape = }, {np.min(burstwindow) = }')
	# print('corr info:', f'{np.max(corr) = }, {np.mean(corr) = }, {corr.shape = }, {np.min(corr) = }')

	if nclip != None or clip != None:
		corr = np.clip(corr, nclip, clip)

	#### Autocorr noise
	autocorr_sigma = np.std( corr[:, sigmawindow[0]:sigmawindow[1]] )

	#### Fit Gaussian to autocorrelation.
	try:
		if popt_custom != []:
			popt, perr = popt_custom, [-1,-1,-1,-1,-1,-1]
		else:
			# popt, pcov = fitgaussiannlsq(corr, p0=p0, sigma=autocorr_sigma, bounds=bounds)
			popt, pcov = fitdatagaussiannlsq(corr, corrextents, p0=p0, sigma=autocorr_sigma, bounds=bounds)
			perr = np.sqrt(np.diag(pcov))
			popt = list(popt); perr = list(perr) # avoid type errors

		if np.isnan(popt).any():
			raise ValueError
		if verbose: print('fit parameters:', popt)
	except (RuntimeError, ValueError):
		if verbose: print('no fit found')
		popt, perr = [-1,-1,-1,-1,-1,-1], [-1,-1,-1,-1,-1,-1]
		if popt_custom != []:
			popt = popt_custom

	x, y = getDataCoords(corrextents, corr.shape)
	fitmap = twoD_Gaussian((y, x), *popt).reshape(corr.shape[0], corr.shape[1])

	# calculate reduced chisquared
	residuals = corr - fitmap
	chisq = np.sum((residuals / autocorr_sigma) ** 2)
	red_chisq = chisq / (corr.shape[0]*corr.shape[1] - len(popt)) # this is chisq/(M-N)

	# Calculate slope
	theta = popt[5] if abs(popt[3]) > abs(popt[4]) else popt[5] - np.pi/2 # defined via eqs. A2 and A3
	slope = np.tan(theta) # MHz/ms
	# print(f'slope calc: {fres_MHz = } {tres_ms = } {theta = } {np.tan(theta) = }')
	theta_err = perr[-1]
	slope_error = (theta_err * (1/np.cos(theta))**2)

	# find center frequency
	center_f = findCenter(burstwindow)*fres_MHz + lowest_freq

	#### Plot
	if plot:
		_plotresult(burstwindow, corr, fitmap, burstkey, center_f, popt, fres_MHz, tres_ms, lowest_freq)

	return (
		slope,
		slope_error,
		popt,
		perr,
		theta,
		red_chisq,
		center_f,
		fitmap
	)

def makeFitmap(popt, corr):
	x, y = np.meshgrid(range(0, corr.shape[1]), range(0, corr.shape[0]))
	fitmap = twoD_Gaussian((y, x), *popt).reshape(corr.shape[0], corr.shape[1])
	return fitmap

# make result headers global
columns = [
	'name',
	'DM',
	'center_f',
	'slope (mhz/ms)',
	'slope error (mhz/ms)',
	'theta',
	'red_chisq',
	'amplitude',
	'xo',
	'yo',
	'sigmax',
	'sigmay',
	'angle',
	'amp_error',
	'xo_error',
	'yo_error',
	'sigmax_error',
	'sigmay_error',
	'angle_error',
	'f_res (mhz)',
	'time_res (s)'
]

def processDMRange(burstname, wfall, burstdm, dmrange, fres_MHz, tres_ms, lowest_freq, p0=[]):
	results = []
	for trialDM in tqdm(dmrange):
		view = np.copy(wfall)
		ddm = trialDM - burstdm
		view = dedisperse(view, ddm, lowest_freq, fres_MHz, tres_ms)

		# bounds = ([-np.inf]*5+ [0], [np.inf]*6) # angle must be positive
		measurement = processBurst(view, fres_MHz, tres_ms, lowest_freq, verbose=False, p0=p0)
		slope, slope_err, popt, perr, theta, red_chisq, center_f, fitmap = measurement
		datarow = [burstname] + [trialDM, center_f, slope, slope_err, theta, red_chisq] + popt + perr + [fres_MHz, tres_ms/1000]
		results.append(datarow)

	df = exportresults(results)
	return results, df

def exportresults(results):
	df = pd.DataFrame(results, columns=columns)
	df = df.set_index('name')
	return df

def plotStampcard(loadfunc, fileglob='*.npy', figsize=(14, 16), nrows=6, ncols=4, twidth=150):
	"""
	Plot bursts and their autocorrelations in a stampcard
	Optionally do the processing to find the slope and plot the resulting fit as well

	loadfunc: a function(filename) that loads the waterfall and returns (subfall, pkidx, wfall)
		where `wfall` is the loaded waterfall,
		`subfall` is the subbanded waterfall that you want to see,
		and `pkidx` is the index in the timeseries of the data with peak intensity
		See `frbprepeaters.loadpsrfits` for an example
	"""
	numfiles = len(glob.glob(fileglob))
	if nrows*ncols != 2*numfiles: # ensure there is enough room for the subplots
		nrows = (2*numfiles // ncols) + 1

	plt.figure(figsize=figsize)
	ploti = itertools.count(start=1, step=1)
	burstnum = 1
	obsdata = None

	for filename in glob.glob(fileglob):
		# print(f'loading {filename}')
		# Handle 2 different types of loadfuncs:
		loadresult = loadfunc(filename)
		if type(loadresult) == tuple and len(loadresult) == 4:
			# eg. loadfunc = frbrepeaters.loadpsrfits
			subfall, pkidx, wfall, obsdata = loadresult
			downf, downt = wfall.shape[0] / subfall.shape[0], wfall.shape[1] / subfall.shape[1]

			df = (obsdata['dfs'][-1] - obsdata['dfs'][0])/len(obsdata['dfs']) * downf  # mhz
			dt = obsdata['dt'][0] / wfall.shape[1]  * downt # s
			dt = dt*1000          # ms
			lowest_freq = obsdata['dfs'][0] # mhz
		else: # eg. loadfunc = np.load
			wfall = loadresult
			subfall = subsample(wfall, 32, wfall.shape[1]//8)
			pkidx = np.nanargmax(np.nanmean(subfall, axis=0))

		twidth = min(twidth, pkidx)-1
		view = subfall[..., pkidx-twidth:pkidx+twidth]
		corr = autocorr2d(view)
		print(f'shape: {view.shape}, \tshape: {corr.shape}')
		if obsdata:
			drift, drift_error, popt, perr, theta, red_chisq, center_f, fitmap = processBurst(view, df, dt, lowest_freq, verbose=False)
			extents, corrextents = getExtents(view, df=df, dt=dt, lowest_freq=lowest_freq)
		else:
			extents, corrextents = None, None

		plt.subplot(nrows, ncols, next(ploti))
		plt.imshow(view, origin='lower', interpolation='none', aspect='auto', extent=extents)
		plt.title(f'Burst #{burstnum}')
		plt.xlabel('time (arb)'), plt.ylabel('freq (arb)')

		plt.subplot(nrows, ncols, next(ploti))
		plt.imshow(corr, origin='lower', interpolation='none', aspect='auto', cmap='gray', extent=corrextents)
		plt.clim(0, np.max(corr)/20)
		plt.title(f'Corr #{burstnum}')
		plt.xlabel('time lag (arb)'), plt.ylabel('freq lag (arb)')
		if obsdata and popt[0] > 0:
			plt.contour(fitmap, [popt[0]/4, popt[0]*0.9], colors='b', alpha=0.75, extent=corrextents, origin='lower')

		burstnum += 1

	plt.tight_layout()
	plt.show()

subburst_suffixes = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k']
def getSubbursts(wfall_cr, df, dt, lowest_freq, regions):
	background = None
	subbursts  = []
	for regname, region in regions.items():
		trange = getExtents(wfall_cr, df=df, dt=dt, lowest_freq=lowest_freq)[0][:2]
		# map time to channel number
		reg_chan = [0, 0]
		for i, edge in enumerate(region):
			reg_chan[i] = round(np.interp(edge, trange, [0, wfall_cr.shape[1]]))

		regiontype =  0 if 'background' in regname else 1
		if regiontype == 0:   # Background
			background = wfall_cr[:, reg_chan[0]:reg_chan[1]]
		elif regiontype == 1: # Burst
			subburst = wfall_cr[:, reg_chan[0]:reg_chan[1]]
			subbursts.append(subburst)

	# pad with background
	# for subburst in subbursts:
	subburstsobj = {}
	suffixes = list(regions.keys())[:-1]
	for subburst, suffix in zip(subbursts, suffixes):
		subburst = np.concatenate((0*background, subburst, 0*background), axis=1)
		subburstsobj[suffix] = subburst
	return subburstsobj

def plotResults(resultsfile, datafiles=[], masks=None, figsize=(14, 16), nrows=6, ncols=4, unitless=False):
	"""
	Similar to plotStampcard but reads all data from the results csv produced by the gui

	Plots fit results by burst at each DM
	"""
	resultsdf = pd.read_csv(resultsfile).set_index('name')
	plt.figure(figsize=figsize)
	ploti = itertools.count(start=1, step=1)
	outputfile = f"{resultsfile.split('.')[0].split('/')[-1]}.pdf"
	print(outputfile)

	try:
		pdf = matplotlib.backends.backend_pdf.PdfPages(outputfile)
	except PermissionError as e:
		return False

	if type(masks) == str: # filename
		if os.path.exists(masks):
			masks = np.load(masks, allow_pickle=True)[0]

	pname = ''
	for name, row in resultsdf.iterrows():
		if pname != name: print('plotting', name)  # print once
		pname = name
		if '_' in name:
			subname = name
			name, suffix = name.split('_')
			regcols = [col for col in row.index if 'reg' in col]
			regcols.append('background') if 'background' in row.index else None
			regions = {suffix: [row[f'regstart_{suffix}'], row[f'regend_{suffix}']]}
			regions['background'] = [0, row['background']]
		else:
			regions = None

		file = [f for f in datafiles if name in f][0]
		data = np.load(file, allow_pickle=True)

		wfall = data['wfall']
		storedshape = wfall.shape

		df, dt_ms = row['f_res (mhz)'], row['time_res (s)']*1000
		lowest_freq = min(data['dfs']) # off by half a channel probably

		# apply masks
		if masks is not None and masks != []:
			mask = [masks[k] for k in masks.keys() if name in k]
			if mask == []: mask = [[]]
			mask = mask[0]
			for m in mask:
				if m < len(wfall):
					wfall[m] = 0

		# Check if the waterfall was subsampled before measuring and subsample if so
		if 'downf' in row and 'downt' in row:
			wfall = subsample(wfall, int(wfall.shape[0]/row['downf']), int(wfall.shape[1]/row['downt']))
		if 'tsamp_width' in row:
			wfall = cropwfall(wfall, twidth=row['tsamp_width']) # crop after subsampling.

		if 'subbg_start (ms)' in row and not np.isnan(row['subbg_start (ms)']):
			tleft, tright = row['subbg_start (ms)'], row['subbg_end (ms)']
			timerange = [0, round(dt_ms*wfall.shape[1])]
			tleft  = round(np.interp(tleft, timerange, [0, wfall.shape[1]]))
			tright = round(np.interp(tright, timerange, [0, wfall.shape[1]]))
			wfall = subtractbg(wfall, tleft, tright)

		extents, corrextents = getExtents(wfall, df=df, dt=dt_ms, lowest_freq=lowest_freq)
		if regions:
			wfall = getSubbursts(wfall, df, dt_ms, lowest_freq, regions)[suffix]
			extents, corrextents = getExtents(wfall, df=df, dt=dt_ms, lowest_freq=lowest_freq)

		# dedisperse
		ddm = row['DM'] - data['DM']
		wfall = dedisperse(wfall, ddm, lowest_freq, df, dt_ms)

		corr = autocorr2d(wfall)
		popt = [row['amplitude'], row['xo'], row['yo'], row['sigmax'], row['sigmay'], row['angle']]

		fitmap = makeFitmap(popt, corr)

		currentplot = next(ploti)
		bname = name if not regions else subname
		aspect = 'auto'
		if unitless:
			corrextents, extent = None, None

		plt.subplot(nrows, ncols, currentplot)
		# print(f"{df/dt_ms = } {wfall.shape[0]/wfall.shape[1] = } {corr.shape[0]/corr.shape[1] = }")
		plt.imshow(wfall, origin='lower', interpolation='none', aspect=aspect, extent=extents)
		plt.axhline(y=row['center_f'], c='k', ls='--', lw=1)
		if popt[0] > 0:
			# add duration and slope markers to waterfall
			pkidx = np.nanargmax(np.nanmean(wfall, axis=0))
			plt.errorbar(pkidx*dt_ms, row['center_f'], xerr=row['tau_w_ms']/2, color='#4b0082',
					 	 linewidth=2, capthick=2, capsize=5, fmt='none')
			xo = pkidx*dt_ms
			x = np.array([xo - row['tau_w_ms'], xo + row['tau_w_ms']])
			plt.plot(x, row['slope (mhz/ms)']*x + row['center_f'] - row['slope (mhz/ms)']*xo, 'b--', lw=1.5)
			plt.xlim(extents[:2]); plt.ylim(extents[2:])
		plt.title(f'{bname}: DM = {round(row["DM"], 2)}')
		plt.xlabel('Time (ms)'), plt.ylabel('Freq (MHz)')

		currentplot = next(ploti)
		plt.subplot(nrows, ncols, currentplot)
		plt.imshow(corr, origin='lower', interpolation='none', aspect=aspect, cmap='gray', extent=corrextents)
		plt.clim(0, np.max(corr)/20)
		plt.title(f'Corr {bname}: DM = {round(row["DM"], 2)}')
		plt.xlabel('time lag (ms)'), plt.ylabel('freq lag (MHz)')
		if popt[0] > 0:
			corrextents = np.array(plt.gca().get_ylim() + plt.gca().get_xlim()) + 0.5 if unitless else corrextents
			plt.contour(fitmap, [popt[0]/4, popt[0]*0.9], colors='b', alpha=0.75, extent=corrextents)
			# plot line of corresponding slope
			xo, yo = (row['xo'] - fitmap.shape[1]/2)*dt_ms, (row['yo'] - fitmap.shape[0]/2)*dt_ms
			lw = 0.8
			if unitless:
				xo, yo = row['xo'], row['yo']
				lw = 1

			def getx(slope):
				if unitless:
					return np.array([corrextents[2] / slope + xo, corrextents[3] / slope + xo])
				else:
					return np.array([corrextents[2] / slope + xo*(1+lw), corrextents[3] / slope + xo*(1+lw)])

			x = getx(row['slope (mhz/ms)'])
			units = 1 if unitless else df/dt_ms
			slope1 = np.tan(popt[5])*units         # row['slope1']
			slope2 = np.tan(popt[5]-np.pi/2)*units # row['slope2']

			# print(f"{round(row['DM'], 2) = } {slope1 = } {slope2 = }")
			x1 = getx(slope1)
			x2 = getx(slope2)
			xp = np.array(corrextents[:2])

			y  = row['slope (mhz/ms)']*(x-xo*(1+lw)) # (1+lw) to account for line thickness when centering
			y1 = slope1*(x1-xo) if unitless else slope1*(x1-xo*(1+lw))
			y2 = slope2*(x2-xo) if unitless else slope2*(x2-xo*(1+lw))

			plt.plot(x, y, 'c--', lw=lw)
			plt.plot(x1, y1, 'm--', lw=lw)
			plt.plot(x2, y2, 'g--', lw=lw)
			plt.xlim(corrextents[:2]); plt.ylim(corrextents[2:])
			# plt.plot(xp, -(1/row['slope (mhz/ms)'])*(xp-xo*(1+lw)), 'w--', lw=lw) # perpendicular axis

		if currentplot == nrows*ncols:
			# save current page of pdf, start new page, reset counter
			plt.tight_layout()
			# plt.show()
			pdf.savefig(plt.gcf())
			plt.close()
			plt.figure(figsize=figsize)
			ploti = itertools.count(start=1, step=1)

	plt.tight_layout()
	pdf.savefig(plt.gcf())
	pdf.close()
	plt.close()
	return True

def _plotresult(burstwindow, corr, fitmap, burstkey, center_f, popt, freq_res, time_res,
				lowest_freq, ploti=None):
	fontsize = 22
	cmap = plt.get_cmap('gray')
	cmap.set_bad(color = 'w', alpha = 1.)

	extents = (0,
			   time_res*burstwindow.shape[1],
			   lowest_freq - freq_res/2.,
			   lowest_freq + freq_res*burstwindow.shape[0])

	corrextents = (-extents[1], extents[1], -(extents[3]-extents[2])*2, (extents[3]-extents[2])*2)

	# extents, corrextents = None, None

	nrows = 7
	aspect = 'auto'
	if ploti == None:
		plt.figure(figsize=(15, 12))
		plt.subplot(121)
	else:
		plt.subplot(nrows, 2, next(ploti))
	plt.title("Burst #{}".format(burstkey), fontsize=fontsize)
	plt.imshow(burstwindow, interpolation='none', aspect=aspect, cmap=cmap, extent=extents, origin='lower')
	# plt.axhline(y=center_f, c='k', ls='--', lw=3)
	plt.xlabel("Time (ms)")
	plt.ylabel("Frequency (MHz)")

	if ploti == None:
		plt.subplot(122)
	else:
		plt.subplot(nrows, 2, next(ploti))
	plt.title("Correlation #{}".format(burstkey), fontsize=fontsize)
	plt.imshow(corr, interpolation='none', aspect=aspect, cmap='gray', extent=corrextents, origin='lower')
	plt.xlabel("Time Shift (ms)")
	plt.ylabel("Frequency Shift (MHz)")
	plt.clim(0, np.max(corr)/20)

	if popt[0] > 0:
		plt.contour(fitmap, [popt[0]/4, popt[0]*0.9], colors='b', alpha=0.75, extent=corrextents,
					origin='lower')

