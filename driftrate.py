import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import numpy as np
import scipy.optimize
import itertools, glob
import pandas as pd
import os
from tqdm import tqdm
import string

def findCenter(burstwindow):
	"""
	Find the center *index* with uncertainties of the peak of the waterfall.
	This is done by time averaging the waterfall to obtain a spectrum and then finding the average channel weighted by the square of the intensity.

	This is calculated as

	.. math::
		\\bar{\\nu_i} = \\frac{\sum_{\\nu_i} \\nu_i I^2(\\nu_i)}{\sum_{\\nu_i} I^2(\\nu_i)}

	The full uncertainty (after taking into account units and standard error propagation) is given by

	.. math::
		\delta \\bar{\\nu} = \sqrt{\\nu_{\\text{res}}^2/12 + 4\sigma_I^2\Big(\\nu_{\\text{res}}\\frac{\sum_{\\nu} (\\nu - \\bar{\\nu})I(\\nu) }{\sum_{\\nu} I^2(\\nu)}\Big)^2}

	Args:
		burstwindow (np.ndarray): 2d array of the waterfall

	Returns:
		tuple: (meanfreqi, errorsum) - tuple of the mean peak frequency index and summation term of the uncertainty.

		You may then calculate the mean frequency in MHz with

		.. code-block:: python

			center_f = meanfreqi*fres_MHz + lowest_freq

		The uncertainty on ``center_f`` is then

		.. code-block:: python

			center_f_err = np.sqrt(fres_MHz**2/12 + 4*wfallsigma**2*(errorsum*fres_MHz)**2)

		``wfallsigma`` (:math:`\sigma_I`) is the standard deviation of the intensity noise. This can be sampled from the noise in the waterfall that doesn't include the burst. For example:

		.. code-block:: python

			wfallsigma = np.std( burstwindow[:, 0:burstwindow.shape[1]//20] )
	"""
	freqspectrum = burstwindow.sum(axis=1)[:, None]
	freqi = np.indices(freqspectrum.shape)[0]
	meanfreqi = np.nansum(freqi*(freqspectrum**2)) / np.nansum(freqspectrum**2)
	errorsum = (np.nansum((freqi - meanfreqi)*(freqspectrum)) / np.nansum(freqspectrum**2))
	return meanfreqi, errorsum

def structureParameter(wfall, dt, tstart, tend):
	"""
	wip. Compute the structure parameter of a burst.

	See eq. 1 in gajjar et al. 2018

	Args:
		wfall (np.ndarray): burst waterfall
		dt (float): time resolution
		tstart (int): time channel to start at
		tend (int): time channel to end at

	Return:
		float: the structure parameter
	"""
	n = (tend - tstart)
	ts = np.nanmean(wfall, axis=0)
	struct = 0
	for i in enumerate(ts[tstart:tend]):
		struct += abs((ts[i] - ts[i+1]) / dt)

	return struct/n

def subband(wfall, nsub):
	"""Downsample frequency channels of a waterfall.

	See :py:meth:`subsample` for more general method.

	Args:
		wfall (np.ndarray): 2d array
		nsub (int): number of channels. nsub should evenly divide ``wfall.shape[0]``

	Returns:
		np.ndarray: 2d subbanded array.
	"""
	nchan, nsamp = wfall.shape
	sub_factor = nchan // nsub
	return np.nanmean(wfall.reshape(-1, sub_factor, nsamp), axis=1)

def subsample(m, nfreq, ntime):
	"""Subsample a waterfall in time and frequency

	Args:
		m (np.ndarray): 2d array to subsample.
		nfreq (int): the number of frequency channels desired. Should evenly divide ``m.shape[0]``
		ntime (int): the number of time channels desired. Should evenly divide ``m.shape[1]``

	Returns:
		np.ndarray: the subsampled array
	"""
	n = np.nanmean(m.reshape(-1, m.shape[0]//nfreq, m.shape[1]), axis=1)
	return np.nanmean(n.reshape(n.shape[0], -1, n.shape[1]//ntime), axis=2)

def subtractbg(wfall, tleft: int=0, tright: int=1):
	"""Subtract a background sample from a waterfall

	This will compute the mean along the time axis to produce an array of noise as a
	function of frequency and subtract that array from the entire waterfall.

	Avoid sampling the burst if possible to improve accuracy of measurements.

	Args:
		wfall (np.ndarray): 2d array of burst waterfall
		tleft (int): time channel to start sample from
		tright (int): time channel to end sample from

	Returns:
		np.ndarray: the waterfall subtracted by the background sample
	"""

	return wfall - wfall[:, tleft:tright].mean(axis=1)[:, None]

def moments(data):
	"""Deprecated function. Returns (height, x, y, width_x, width_y)
	initial gaussian parameters of a 2D distribution by calculating its
	moments using integer channel coordinates.

	Args:
		data (np.ndarray): the data to compute moments of

	Returns:
		list: (height, x, y, sigma_x, sigma_y)
	"""
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

def gaussian_dt(point, amplitude, t0, dt, nu0, sigma_t, sigma_nu):
	"""Gaussian model in terms of the temporal drift $d_t$ in ms/MHz

	Equation 2 of Jahns et al. (2023)
	$sigma_t$ is the burst width measured at $t_0$
	$sigma_\\nu$ is the overall bandwidth of the summed 1D spectrum
	"""
	nu, t = point
	g = amplitude * np.exp(-(t - t0 - dt*(nu-nu0))**2/(2*sigma_t**2) - (nu-nu0)**2/(2*sigma_nu**2))
	return g.ravel()

def gaussian_dnu(point, amplitude, t0, dnu, nu0, w_t, w_nu):
	"""Gaussian model in terms of the frequency drift $d_\\nu$ in MHZ/ms

	Equation 10 of Jahns et al. (2023)
	$w_t$ is equivalent to the duration measured from the summed 1D time series
	$w_\\nu$ is the bandwidth at $t_0$
	"""
	nu, t = point
	g = amplitude * np.exp(-(t-t0)**2 / (2*w_t**2) - (nu - nu0 - dnu*(t-t0))**2/(2*w_nu**2))
	return g.ravel()

def twoD_Gaussian(point, amplitude, xo, yo, sigma_x, sigma_y, theta):
	"""2D rotatable Gaussian Model function. Used as the model function in scipy.optimize.curve_fit

	Args:
		point (tuple): (y, x) point to evaluate the model at
		amplitude (float): amplitude of the 2D gaussian
		xo (float): central x position of the gaussian
		yo (float): central y position of the gaussian
		sigma_x (float): standard deviation in x-direction
		sigma_y (float): standard deviation in y-direction
		theta (float): angle of gaussian

	Returns:
		list: 1d array of raveled 2d gaussian data
	"""
	y, x = point
	xo = float(xo)
	yo = float(yo)
	a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
	b = (np.sin(2*theta))/(2*sigma_x**2) - (np.sin(2*theta))/(2*sigma_y**2)
	c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
	g = amplitude*np.exp( - a*((x-xo)**2) - b*(x-xo)*(y-yo) - c*((y-yo)**2))
	return g.ravel()

def twoD_Gaussian_floor(point, amplitude, xo, yo, sigma_x, sigma_y, theta, floor):
	return twoD_Gaussian(point, amplitude, xo, yo, sigma_x, sigma_y, theta) + floor

def getDataCoords(extents, shape):
	x = np.linspace(extents[0], extents[1], num=shape[1])
	y = np.linspace(extents[2], extents[3], num=shape[0])
	x, y = np.meshgrid(x, y)
	return x, y

def fitdatagaussiannlsq(
	data,
	extents,
	p0=[],
	sigma=0,
	bounds=(-np.inf, np.inf),
	model=twoD_Gaussian
):
	"""Fit 2d gaussian to data with the non-linear least squares algorithm implemented by `scipy.optimize.curve_fit <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html>`_.

	Uses the physical units as the data's coordinates.

	Args:
		data (np.ndarray): 2D array of the data to fit on
		extents (tuple): Tuple of the initial time, final time, initial frequency, and final frequency (ti, tf, nu_i, nu_f). See :py:meth:`getExtents()`.
		p0 (list, optional): initial guess to used, input as a list of floats corresponding to [amplitude, x0, y0, sigmax, sigmay, theta]. See :py:meth:`twoD_Gaussian`.
		sigma (float, optional): uncertainty to use during fitting. See `Scipy Documentation <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html>`_ for more information
		bounds (tuple, optional): lower and upper bounds on parameters. See `Scipy Documentation <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html>`_ for more information
		model (function, optional): function to use. ``twoD_Gaussian`` by default.

	Returns:
		tuple: tuple of (popt, pcov) where popt is the list of gaussian parameters found and pcov is the covariance matrix.

	"""
	# use curve-fit (non-linear leastsq)
	# x = range(0, data.shape[1]); y = range(0, data.shape[0])
	data = data / np.max(data) # normalize
	x, y = getDataCoords(extents, data.shape)
	p0 = [1,0,0,1,100,0] if p0 == [] else p0
	sigma = np.zeros(len(data.ravel())) + sigma
	print(">> fitdatagaussianlsq", model.__name__)
	popt, pcov = scipy.optimize.curve_fit(
		model,
		(y, x),
		data.ravel(),
		p0=p0,
		sigma=sigma,
		absolute_sigma=True,
		bounds=bounds
	)
	return popt, pcov

def makeDataFitmap(popt, corr, extents, model=twoD_Gaussian):
	"""Make a fitmap using physical coordinates

	Args:
		popt (list): List of 2D gaussian model parameters. e.g. [amplitude, xo, yo, sigmax, sigmay, theta)]
		corr (np.ndarray): the 2d autocorrelation. Can simply pass a 2d array of the same size
		extents (list): Return value of :py:meth:`getExtents`
		model (function, optional): the model function corresponding to ``popt``

	Returns:
		np.ndarray: a 2d array of the 2d gaussian specified by popt

	"""
	x, y = getDataCoords(extents, corr.shape)
	fitmap = model((y, x), *popt).reshape(corr.shape[0], corr.shape[1])
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
	a_dm = 1. / 2.41e-4
	dedisp = np.zeros_like(wfall)

	ref_freq = freq[0]### ORIGINAL
	# ref_freq = freq[-1]
	# print("ref_freq", ref_freq)

	shift = (a_dm * dm * (ref_freq ** -2 - freq ** -2) / dt) ### ORIGINAL (low freq anchor)
	# shift = (a_dm * dm * (freq ** -2 - ref_freq ** -2) / dt)
	shift = shift.round().astype(int)

	for i, ts in enumerate(wfall):
		dedisp[i] = np.roll(ts, shift[i])

	return dedisp

a_dm = 4.14937759336e6 # MHz^2 ms cm^3/pc
def dedisperse(intensity, DM, nu_low, df_mhz, dt_ms, cshift=0, a_dm=a_dm):
	"""Incoherently dedisperse a dynamic spectrum to the specified DM.

	Computes the time delay at each frequency and correspondingly rolls the data

	.. math::
		\Delta t = - a\Big(\\frac{1}{\\nu^2_i} - \\frac{1}{\\nu^2_\\text{high}}\Big)\\text{DM}

	By default uses the pulsar community dispersion constant of

	:math:`a = 4.14937759336e6 \quad\\text{MHz}^2 \\text{cm}^3 \\text{pc}^-1 \\text{ms}`.

	Args:
		intensity (np.ndarray): the dynamic spectrum (or waterfall) to dedisperse. ``intensity[0]`` should correspond to the lowest frequency.
		DM (float): the dispersion measure in pc/cm^3
		nu_low (float): the frequency at the bottom of the band
		df_mhz (float): the frequency resolution of the channels in MHZ
		dt_ms (float): the time resolution of channels in ms
		cshift (int, optional): additional horizontal shift to add after dedispersion in number of channels
		a_dm (float): the dispersion constant to use when dedispersing. See units above.

	Returns:
		np.ndarray: the dedispersed 2d intensity array
	"""
	dedispersed = np.copy(intensity)

	high_ref_freq = nu_low + len(dedispersed)*df_mhz # half channel?
	low_ref_freq  = nu_low
	# a_dm = 4.1488064239e6 # kulkarni
	for i, row in enumerate(dedispersed): # i == 0 corresponds to bottom of the band
		nu_i = nu_low + i*df_mhz
		# High frequency anchor
		deltat = - a_dm * (nu_i**-2 - high_ref_freq**-2) * DM

		# Low frequency anchor
		#deltat = 4.14937759336e6 * (low_ref_freq**-2 - nu_i**-2) * DM

		channelshift = int(round(deltat/dt_ms))
		dedispersed[i] = np.roll(dedispersed[i], channelshift)

	# optionally center view
	dedispersed = np.roll(dedispersed, cshift, axis=1)

	return dedispersed

def getExtents(wfall, df:float=1.0, dt:float=1.0, lowest_freq:float=1.0, lowest_time:float=0.0):
	"""Given a waterfall's time and frequency resolutions, return the extents of the axes as well as the axes of its autocorrelation.

	Convenience function for plt.imshow's extent keyword

	Args:
		wfall (np.ndarray): 2D array of the waterfall
		df (float): frequency resolution
		dt (float): time resolution
		lowest_freq (float): lowest frequency in the band

	Returns:
		tuple: (extents, corrextents) where extents is a list of the waterfall extents and corrextents is a list of the corresponding autocorrelation extents
	"""
	extents = (lowest_time,
			   lowest_time + dt*wfall.shape[1],
			   lowest_freq - df/2,
			   lowest_freq + df*wfall.shape[0] + df/2)

	# corrextents = (-extents[1], extents[1], -(extents[3]-extents[2]), (extents[3]-extents[2]))
	corrextents = (-(extents[1]-extents[0]), (extents[1]-extents[0]), -(extents[3]-extents[2]), (extents[3]-extents[2]))
	return extents, corrextents

def cropwfall(wfall, twidth=150, pkidx=None):
	"""
	Crops a waterfall and attempts to center the burst in the array.

	Args:
		wfall (np.ndarray): the array to be cropped of size M x N
		twidth (int): the number of channels on either side of the burst. Total width is therefore 2*twidth
		pkidx (int, optional): the index of the waterfall to center. If None will use np.argmax of the timeseries as the index to center.

	Returns:
		np.ndarray: cropped waterfall of size M x (2*twidth)
	"""
	twidth = round(twidth)
	if twidth <= 0 or twidth > wfall.shape[1]:
		twidth = wfall.shape[1]
	wfall = wfall.copy()
	ts    = np.nanmean(wfall, axis=0)
	if not pkidx:
		pkidx = np.nanargmax(ts)

	# try making a window with pkidx in the middle, and shift if window is too large
	# window will always be of width 2*twidth
	ledge, redge = pkidx-twidth, pkidx+twidth
	if ledge < 0:
		if redge + abs(ledge) <= wfall.shape[1]:
			redge += abs(ledge)
		ledge = 0
	if redge > wfall.shape[1]:
		if ledge - (redge - wfall.shape[1]) >= 0:
			ledge -= redge
		redge = wfall.shape[1]

	return wfall[..., ledge:redge]

def updatenpz(npz, field, val):
	"""Update a field in a npz file and resave it.

	Args:
		npz (str): filename of .npz file
		field (str): the field you would like to update
		val (Any): the value to update the field to

	Returns:
		None

	"""
	data = np.load(npz)
	newdata = {}
	for key in data.files:
		if key == field:
			newdata[key] = val
		else:
			newdata[key] = data[key]
	np.savez(npz, **newdata)

def autocorr2d(data, maskcorrpeak=True):
	"""Returns a 2D autocorrelation computed via an intermediate FFT

	Args:
		data (np.ndarray): 2D array of size MxN
		maskcorrpeak (bool, optional): Set to False to keep the zero lag peak

	Returns:
		np.ndarray: 2D array of size 2M-1 x 2N-1. The autocorrelation of ``data``.
	"""

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

	corr = temp_array_b[:-1,:-1]#/float(nx*ny)
	if maskcorrpeak:
		ind = np.unravel_index(np.argmax(corr), corr.shape)
		if corr[ind] == np.max(corr): # double check
			corr[ind] = 0
			if np.max(corr) > 5*np.median(corr): # Repeat once
				ind = np.unravel_index(np.argmax(corr), corr.shape)
				corr[ind] = 0

	return corr

def _fakeburst(shape, popt, noisefrac=0.2):
	x, y = getDataCoords([0, shape[0], 0, shape[1]], shape)
	wfall = twoD_Gaussian((y, x), *popt).reshape(shape[0], shape[1])
	wfall += noisefrac*np.random.random(wfall.shape)
	wfall = subtractbg(wfall, 0, 0.25*shape[1])
	corr = autocorr2d(wfall)
	return wfall, corr

def processBurst(
	burstwindow,
	fres_MHz,
	tres_ms,
	lowest_freq,
	burstkey=1,
	p0=[],
	popt_custom=[],
	bounds=(-np.inf, np.inf),
	nclip=None,
	clip=None,
	plot=False,
	corrsigma=None,
	wfallsigma=None,
	maskcorrpeak=True,
	verbose=True,
	lowest_time=0,
	usefloor=False
):
	"""
	Given a waterfall of a burst, will use the 2d autocorrelation+gaussian fitting method
	to perform spectro-temporal measurements of the burst

	Can optionally plot the measurement.

	Args:
		burstwindow (np.ndarray): 2d array of the burst waterfall
		fres_MHz (float): frequency resolution in MHz
		tres_ms (float): time resolution in ms
		lowest_freq (float): lowest frequency in MHz
		burstkey (int, optional): Burst number. Used in plot title
		p0 (list, optional): Initial 2d gaussian guess to use. [amplitude, x0, y0, sigmax, sigmay, theta]
		popt_custom (list, optional): Override the fit and plot your own 2D gaussian by placing your popt list here
		bounds (tuple, optional): parameter bounds. See :py:meth:`fitdatagaussianlsq`
		nclip (float, optional): minimum clip value of autocorrelation. Applied before fitting.
		clip (float, optional): maximum clip value of autocorrelation. Applied before fitting.
		plot (bool, optional): If true will display a diagnostic plot of the fit
		corrsigma (float, optional): Standard deviation of noise of autocorrelation. Used when fitting.
		wfallsigma (float, optional): Standard deviation of noise of waterfall. Used when fitting.
		verbose (bool, optional): Set to False to limit console output
		maskcorrpeak (bool, optional): Set to False to keep the zero lag correlation peak
		lowest_time (float, optional): starting time (ms) of waterfall. Default 0.
		usefloor (bool, optional): Whether or not to add a constant floor as a parameter in the gaussian ACF fit.

	Returns:
		tuple: (slope, slope_error, popt, perr, theta, red_chisq, center_f, center_f_err, fitmap)

	"""

	corr = autocorr2d(burstwindow, maskcorrpeak=maskcorrpeak)
	_, corrextents = getExtents(
		burstwindow,
		df=fres_MHz,
		dt=tres_ms,
		lowest_freq=lowest_freq,
		lowest_time=lowest_time
	)
	# print('wfall info:', f'{np.max(burstwindow) = }, {np.mean(burstwindow) = }, {burstwindow.shape = }, {np.min(burstwindow) = }')
	# print('corr info:', f'{np.max(corr) = }, {np.mean(corr) = }, {corr.shape = }, {np.min(corr) = }')

	if nclip != None or clip != None:
		corr = np.clip(corr, nclip, clip)

	#### Autocorr noise
	if not corrsigma:
		corrsigma = (0, 50)
		autocorr_sigma = np.std( corr[:, corrsigma[0]:corrsigma[1]] )
	elif type(corrsigma) == tuple or type(corrsigma) == list:
		autocorr_sigma = np.std( corr[:, corrsigma[0]:corrsigma[1]] )
	else:
		autocorr_sigma = corrsigma

	#### Fit Gaussian to autocorrelation.
	try:
		if popt_custom != []:
			popt, perr = popt_custom, [-1,-1,-1,-1,-1,-1]
		elif not usefloor:
			popt, pcov = fitdatagaussiannlsq(
				corr,
				corrextents,
				p0=p0,
				sigma=autocorr_sigma,
				bounds=bounds
			)
			perr = np.sqrt(np.diag(pcov))
			popt = list(popt); perr = list(perr) # avoid type errors
		elif usefloor:
			popt, pcov = fitdatagaussiannlsq(
				corr,
				corrextents,
				p0=[1,0,0,1,100,0,0] if p0 == [] else p0+[0],
				sigma=autocorr_sigma,
				bounds=bounds,
				model=twoD_Gaussian_floor
			)
			perr = np.sqrt(np.diag(pcov))
			popt = list(popt); perr = list(perr) # avoid type errors


		if np.isnan(popt).any():
			raise ValueError
		if verbose: print('fit parameters:', popt)
	except (RuntimeError, ValueError):
		if verbose: print('no fit found')
		popt, perr = [-1,-1,-1,-1,-1,-1], [-1,-1,-1,-1,-1,-1]
		if usefloor:
			popt.append(-1)
			perr.append(-1)
		if popt_custom != []:
			popt = popt_custom

	x, y = getDataCoords(corrextents, corr.shape)
	if not usefloor:
		fitmap = twoD_Gaussian((y, x), *popt).reshape(corr.shape[0], corr.shape[1])
	elif usefloor:
		print(popt)
		fitmap = twoD_Gaussian_floor((y, x), *popt).reshape(corr.shape[0], corr.shape[1])

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
	center_i, errorsumi = findCenter(burstwindow)
	center_f = center_i*fres_MHz + lowest_freq
	errorsum = errorsumi*fres_MHz
	if not wfallsigma:
		wfallsigma = np.std( burstwindow[:, 0:burstwindow.shape[1]//20] ) # use the first 5% of channels
	center_f_err = np.sqrt(fres_MHz**2/12 + 4*wfallsigma**2*errorsum**2)

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
		center_f_err,
		fitmap
	)

def makeFitmap(popt, corr):
	"""Deprecated function. See :py:meth:`makeDataFitmap`"""
	x, y = np.meshgrid(range(0, corr.shape[1]), range(0, corr.shape[0]))
	fitmap = twoD_Gaussian((y, x), *popt).reshape(corr.shape[0], corr.shape[1])
	return fitmap

# make result headers global
columns = [
	'name',
	'DM',
	'center_f',
	'center_f_err',
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
	# 'pkidx'
]

def processDMRange(burstname, wfall, burstdm, dmrange, fres_MHz, tres_ms, lowest_freq, p0=[],
				   corrsigma=None, wfallsigma=None, tqdmout=None, progress_cb=None, lowest_time=0):
	"""Process burst measurements over a range of DMs.

	Dedisperses to each DM in ``dmrange`` and calls :py:meth:`processBurst` on the resulting waterfall. Returns a list of all measurements as well as a dataframe. Columns of the dataframe are given by ``driftrate.columns``

	Args:
		burstname (str): name of burst to use in results dataframe
		wfall (np.ndarray): 2d array of burst waterfall
		burstdm (float): the DM of the burst in ``wfall``
		dmrange (list[float]): a list of the DMs you would measurements at
		fres_MHz (float): frequency resolution in MHz
		tres_ms (float): time resolution in milliseconds
		lowest_freq (float): lowest frequency of the waterfall
		p0 (list, optional): optionally provide an initial guess for the 2d gaussian fit with list of [amplitude, x0, y0, sigmax, sigmay, theta]
		corrsigma (float, optional): standard deviation of the burst's autocorrelation
		wfallsigma (float, optional): standard deviation of the burst waterfall
		tqdmout (object, optional): output for TQDM's progress bar
		progress_cb (function, optional): callback to run after each DM is processed. Called with ( #of DMs processed) / len(dmrange) and a string of "{#}/{len(dmrange)}"
		lowest_time (float, optional): the starting time of the waterfall, passed to :py:meth:`processBurst`.

	Returns:
		tuple: (list of measurement results, dataframe of measurement results)

	"""
	results = []
	prog = 0
	for trialDM in tqdm(dmrange, file=tqdmout):
		prog += 1
		view = np.copy(wfall)
		ddm = trialDM - burstdm
		view = dedisperse(view, ddm, lowest_freq, fres_MHz, tres_ms)

		# bounds = ([-np.inf]*5+ [0], [np.inf]*6) # angle must be positive
		measurement = processBurst(view, fres_MHz, tres_ms, lowest_freq, verbose=False, p0=p0,
									corrsigma=corrsigma, wfallsigma=wfallsigma,
									lowest_time=lowest_time)
		slope, slope_err, popt, perr, theta, red_chisq, center_f, center_f_err, fitmap = measurement
		datarow = [burstname] + [trialDM, center_f, center_f_err, slope, slope_err, theta, red_chisq] + popt + perr + [fres_MHz, tres_ms/1000]
		results.append(datarow)
		p0=popt
		if progress_cb:
			progress_cb(prog/len(dmrange), f"{prog}/{len(dmrange)}")

	df = exportresults(results)
	return results, df

def exportresults(results):
	"""Creates a dataframe of measurement results.

	See ``driftrate.columns``.

	Args:
		results (list): the results list
	"""
	df = pd.DataFrame(results, columns=columns)
	df = df.set_index('name')
	return df

def plotStampcard(loadfunc, fileglob='*.npy', figsize=(14, 16), nrows=6, ncols=4, twidth=150):
	"""	Plot bursts and their autocorrelations in a stampcard

	Optionally do custom processing to find the slope and plot the resulting fit as well.

	Args:
		loadfunc (function): a function that accepts "filename" to a waterfall as an argument and loads the waterfall and returns (subfall, pkidx, wfall) where ``wfall`` is the loaded waterfall, ``subfall`` is the subbanded waterfall that you want to see, and ``pkidx`` is the index in the timeseries of the data with peak intensity. Implement this yourself to load your data files
		fileglob (str): a glob that matches on your data files (e.g. "\*.npy")
		figsize (tuple[float]): Tuple of (width, height) for the resulting figure
		nrows (int): number of rows in the stampcard
		ncols (int): number of columns in the stampcard
		twidth (int): number of channels on either side of the burst to plot

	Returns:
		None: shows a figure. Useful in jupyterlab for example.
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
			drift, drift_error, popt, perr, theta, red_chisq, center_f, center_f_err, fitmap = processBurst(view, df, dt, lowest_freq, verbose=False)
			extents, corrextents = getExtents(view, df=df, dt=dt, lowest_freq=lowest_freq)
		else:
			extents, corrextents = None, None

		plt.subplot(nrows, ncols, next(ploti))
		plt.imshow(view, origin='lower', interpolation='none', aspect='auto', extent=extents)
		plt.title(f'Burst #{burstnum}')
		plt.xlabel('time (arb)'), plt.ylabel('freq (arb)')

		plt.subplot(nrows, ncols, next(ploti))
		plt.imshow(corr, origin='lower', interpolation='none', aspect='auto', cmap='gray', extent=corrextents)
		plt.clim(0, np.max(corr))
		plt.title(f'Corr #{burstnum}')
		plt.xlabel('time lag (arb)'), plt.ylabel('freq lag (arb)')
		if obsdata and popt[0] > 0:
			plt.contour(fitmap, [popt[0]/4, popt[0]*0.9], colors='b', alpha=0.75, extent=corrextents, origin='lower')

		burstnum += 1

	plt.tight_layout()
	plt.show()

subburst_suffixes = list(string.ascii_lowercase)
subburst_suffixes += [''.join(p) for p in list(zip(list(string.ascii_lowercase), list(string.ascii_lowercase)))]
def getSubbursts(wfall_cr, df, dt, lowest_freq, regions):
	"""Split a waterfall into regions.

	Specify regions with a dictionary. For example

	.. code-block:: python

		regions = {"a": [0, 3], "b": [3, 5.5]}

	where the arrays are the timestamps in milliseconds.

	Regions should be named as "a", "b", "c", etc. See ``driftrate.subburst_suffixes``.

	Meant for use by frbgui.

	Args:
		wfall_cr (np.ndarray): waterfall to split
		df (float): frequency resolution
		dt(float): time resolution
		lowest_freq (float): lowest frequency
		regions (dict): Dictionary of {"a": [start_ms, end_ms], etc..}

	Returns:
		dict: A dictionary of cropped waterfalls. Follows {"a": np.ndarray}
	"""
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

def readRegions(resultsdf):
	"""Parse regions out of a results csv produced by FRBGui

	Args:
		resultsdf (pd.DataFrame): dataframe of results. Can load with ``pd.read_csv(filename).set_index('name')``
	"""
	if 'background' not in resultsdf.columns:
		return {}
	regionsdf = resultsdf.loc[~np.isnan(resultsdf['background'])]
	bursts_with_regions = regionsdf.index.unique()

	regionsobj = {}
	for burst in bursts_with_regions:
		name, suffix = '_'.join(burst.split('_')[:-1]), burst.split('_')[-1]
		if suffix not in subburst_suffixes:
			regionsobj[burst] = {}
			regionsobj[burst]['background'] = [0, float(resultsdf.loc[burst][['background']].iloc[0])]
		else:
			regionsobj[name][suffix] = list(resultsdf.loc[name][[f'regstart_{suffix}', f'regend_{suffix}']].iloc[0])
	return regionsobj

def scilabel(num, err):
	"""Utility for pretty scientific notation labels

	e.g. (6.0 +/- 0.3) $\\times 10^2$ MHz/ms

	Args:
		num (float): the number to label
		err (float): the uncertainty of ``num``

	Returns:
		str: the formatted label string
	"""
	sign = '' if num > 0 else '-'
	num = abs(num)
	return (
		f'({sign}{num/(10**np.floor(np.log10(num))):.1f}'
		f'$\\pm$ {err/(10**np.floor(np.log10(num))):.1f})'
		f'$\\times 10^{{{np.floor(np.log10(num)):.0f}}}$'
	)

def plotResults(resultsfile, datafiles=[], masks=None, figsize=(14, 16), nrows=6, ncols=4, clip=1,
				snrLines=False, show=False):
	"""Given a results CSV produced by FRBGui will plot fit results by burst at each DM
	in a stampcard and save a PDF with the same name.

	Args:
		resultsfile (str): Filename of CSV file
		datafiles (list[str]): list of burst filenames in FRBGui's .npz :ref:`burstformat`
		masks (str): filename to FRBGui maskfile
		figsize (tuple(float)): (width, height) of figure
		nrows (int): number of rows in stampcard
		ncols (int): number of cols in stampcard
		clip (int): factor to clip the color scale of the autocorrelation figure for improving display SNR
		snrLines (bool): If true plots 1 sigma lines around the autocorrelation
		show (bool): Show the figure in a window if true. Displays a pdf otherwise.

	Returns:
		bool: True if completed. Saves a PDF file when completed.
	"""
	resultsdf = pd.read_csv(resultsfile).set_index('name')
	plt.figure(figsize=figsize)
	ploti = itertools.count(start=1, step=1)
	outputfile = f"{resultsfile.split('.')[0]}.pdf"
	print(outputfile)

	try:
		pdf = matplotlib.backends.backend_pdf.PdfPages(outputfile)
	except PermissionError as e:
		return "Permission Denied"

	if type(masks) == str: # filename
		if os.path.exists(masks):
			masks = np.load(masks, allow_pickle=True)[0]

	pname = ''
	for name, row in resultsdf.iterrows():
		if pname != name: print('plotting', name)  # print once
		pname = name
		ismulti = any([suffix in str(name)[-2:] for suffix in subburst_suffixes])
		if 'background' in row.index and not np.isnan(row['background']) and ismulti:
			subname = name
			name, suffix = '_'.join(name.split('_')[:-1]), name.split('_')[-1]
			regcols = [col for col in row.index if 'reg' in col]
			regcols.append('background') if 'background' in row.index else None
			regions = {suffix: [row[f'regstart_{suffix}'], row[f'regend_{suffix}']]}
			regions['background'] = [0, row['background']]
		else:
			regions = None

		file = [f for f in datafiles if name in os.path.basename(f)][0]
		data = np.load(file, allow_pickle=True)

		wfall = data['wfall']
		storedshape = wfall.shape

		df, dt_ms = row['f_res (mhz)'], row['time_res (s)']*1000
		lowest_freq = min(data['dfs']) # off by half a channel probably

		# apply masks
		if masks is not None and masks != []:
			mask = [masks[k] for k in masks.keys() if name in k]
			if mask != []:
				mask = mask[0]['chans']
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

		fitmap = makeDataFitmap(popt, corr, corrextents)

		currentplot = next(ploti)
		bname = name if not regions else subname
		if len(bname) > 15: # limit num characters
			bname = '[...]'+bname[-15:]
		aspect = 'auto'

		plt.subplot(nrows, ncols, currentplot)
		# print(f"{df/dt_ms = } {wfall.shape[0]/wfall.shape[1] = } {corr.shape[0]/corr.shape[1] = }")
		plt.imshow(wfall, origin='lower', interpolation='none', aspect=aspect, extent=extents)
		plt.axhline(y=row['center_f'], c='k', ls='--', lw=1)
		if popt[0] > 0:
			# add duration and slope markers to waterfall
			pkidx = np.nanargmax(np.nanmean(wfall, axis=0))
			plt.errorbar(pkidx*dt_ms, row['center_f'],
				xerr=row['tau_w_ms']/2, yerr=row['bandwidth (mhz)']/2,
				color='#4b0082', linewidth=2, capthick=2, capsize=5, fmt='none')
			xo = pkidx*dt_ms
			x = np.array([xo - row['tau_w_ms'], xo + row['tau_w_ms']])
			plt.plot(x, row['slope (mhz/ms)']*x + row['center_f'] - row['slope (mhz/ms)']*xo, 'b--', lw=1.5)
			plt.xlim(extents[:2]); plt.ylim(extents[2:])
		plt.title(f'{bname}: DM = {round(row["DM"], 2)}')
		plt.xlabel('Time (ms)'), plt.ylabel('Freq (MHz)')

		currentplot = next(ploti)
		plt.subplot(nrows, ncols, currentplot)
		plt.imshow(corr, origin='lower', interpolation='none', aspect=aspect, cmap='gray', extent=corrextents)
		plt.clim(0, np.max(corr)/clip)
		plt.title(f'Corr {bname}: DM = {round(row["DM"], 2)}')
		plt.xlabel('time lag (ms)'), plt.ylabel('freq lag (MHz)')
		if popt[0] > 0:
			plt.contour(fitmap, [popt[0]/4, popt[0]*0.9], colors='b', alpha=0.75, extent=corrextents)
			# find and display snr
			a = round(np.interp(popt[1]-popt[3], [corrextents[0], corrextents[1]], [0, corr.shape[1]]))
			b = round(np.interp(popt[1]+popt[3], [corrextents[0], corrextents[1]], [0, corr.shape[1]]))
			c = round(np.interp(popt[2]-popt[4], [corrextents[2], corrextents[3]], [0, corr.shape[0]]))
			d = round(np.interp(popt[2]+popt[4], [corrextents[2], corrextents[3]], [0, corr.shape[0]]))
			if snrLines:
				plt.axvline(x=corrextents[0] + 2*popt[3], c='w', ls='--')
				plt.axhline(y=corrextents[2] + 2*popt[4], c='w', ls='--')

			# snr = abs(corr[c:d, a:b].sum() / corr[0:(d-c), 0:(b-a)].sum())
			# snr = 2*popt[0]*(1-np.exp(-1/2)) / (corr[0:(d-c), 0:(b-a)].sum() / ((d-c)*(b-a)))
			snr = abs(2*popt[0]*(1-np.exp(-1/2)) / (corr[0:(d-c), 0:(b-a)].sum() / (4*popt[4]*popt[3])))
			if snr < 100000:
				plt.text(corrextents[0]*0.95, corrextents[2]*0.95, f"snr: {snr:.1f}", color='w')
			else:
				plt.text(corrextents[0]*0.95, corrextents[2]*0.95, f"snr: > 100000", color='w')

			# plot line of corresponding slope
			xo, yo = row['xo'], row['yo']
			slope = row['slope (mhz/ms)']
			x = np.array([corrextents[2] / slope + xo, corrextents[3] / slope + xo])
			plt.plot(x, slope*(x-xo), 'g--', lw=1.5)
			plt.xlim(corrextents[:2]); plt.ylim(corrextents[2:])

		if currentplot == nrows*ncols:
			# save current page of pdf, start new page, reset counter
			plt.tight_layout()
			if show: plt.show()
			pdf.savefig(plt.gcf())
			plt.close()
			plt.figure(figsize=figsize)
			ploti = itertools.count(start=1, step=1)

	plt.tight_layout()
	pdf.savefig(plt.gcf())
	if show: plt.show()
	pdf.close()
	plt.close()
	return True

def _plotresult(burstwindow, corr, fitmap, burstkey, center_f, popt, freq_res, time_res,
				lowest_freq, ploti=None):
	fontsize = 'large'
	extents = (
		0,
	   time_res*burstwindow.shape[1],
	   lowest_freq - freq_res/2.,
	   lowest_freq + freq_res*burstwindow.shape[0]
   )

	corrextents = (-extents[1], extents[1], -(extents[3]-extents[2])*2, (extents[3]-extents[2])*2)

	# extents, corrextents = None, None
	nrows = 7
	aspect = 'auto'
	if ploti == None:
		plt.figure()
		plt.subplot(121)
	else:
		plt.subplot(nrows, 2, next(ploti))
	plt.title("Burst #{}".format(burstkey), fontsize=fontsize)
	plt.imshow(
		burstwindow,
		interpolation='none',
		aspect=aspect,
		extent=extents,
		origin='lower'
	)
	# plt.axhline(y=center_f, c='k', ls='--', lw=3)
	plt.xlabel("Time (ms)")
	plt.ylabel("Frequency (MHz)")

	if ploti == None:
		plt.subplot(122)
	else:
		plt.subplot(nrows, 2, next(ploti))
	plt.title("Correlation #{}".format(burstkey), fontsize=fontsize)
	plt.imshow(
		corr,
		interpolation='none',
		aspect=aspect,
		cmap='gray',
		extent=corrextents,
		origin='lower'
	)
	plt.xlabel("Time Shift (ms)")
	plt.ylabel("Frequency Shift (MHz)")
	plt.clim(0, np.max(corr)/1)

	if popt[0] > 0:
		plt.contour(
			fitmap,
			[popt[0]/4, popt[0]*0.9],
			colors='b',
			alpha=0.75,
			extent=corrextents,
			origin='lower'
		)
	plt.tight_layout()

