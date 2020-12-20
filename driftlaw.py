import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.odr

def computeModelDetails(frame):
	tauwerror_expr = lambda r: 1e3*r['time_res']*np.sqrt(r['max_sigma']**6*r['min_sigma_error']**2*np.cos(r['angle']-np.pi/2)**4 + r['angle_error']**2*r['max_sigma']**2*r['min_sigma']**2*(-r['max_sigma']**2 + r['min_sigma']**2)**2*np.cos(r['angle']-np.pi/2)**2*np.sin(r['angle']-np.pi/2)**2 + r['max_sigma_error']**2*r['min_sigma']**6*np.sin(r['angle']-np.pi/2)**4)/(r['max_sigma']**2*np.cos(r['angle']-np.pi/2)**2 + r['min_sigma']**2*np.sin(r['angle']-np.pi/2)**2)**1.5

	frame['drift_abs'] = -1*(frame['drift (mhz/ms)'])
	frame['drift_over_nuobs'] = frame[['drift_abs','center_f']].apply(lambda row: row['drift_abs'] / row['center_f'], axis=1)
	frame['min_sigma'] = frame[['sigmax','sigmay']].apply(lambda row: min(abs(row['sigmax']), abs(row['sigmay'])), axis=1)
	frame['max_sigma'] = frame[['sigmax','sigmay']].apply(lambda row: max(abs(row['sigmax']), abs(row['sigmay'])), axis=1)
	# the following two lines assume that if sigmax > sigmay, then sigmax_error > sigmay_error, which is true (so far) for this dataset
	frame['min_sigma_error'] = frame[['sigmax_error','sigmay_error']].apply(lambda row: min(row['sigmax_error'], row['sigmay_error']), axis=1)
	frame['max_sigma_error'] = frame[['sigmax_error','sigmay_error']].apply(lambda row: max(row['sigmax_error'], row['sigmay_error']), axis=1)

	frame['sigma_t']   = frame[['min_sigma','time_res']].apply(lambda row: row['min_sigma']*row['time_res'], axis=1)

	frame['tau_w'] = frame[['time_res', 'min_sigma', 'max_sigma', 'angle']].apply(
		lambda r: r['time_res']*r['min_sigma']*r['max_sigma'] / np.sqrt( np.abs((np.sin(r['angle']-np.pi/2)*r['min_sigma'])**2 + (np.cos(r['angle']-np.pi/2)*r['max_sigma'])**2 )),
		axis=1
	)

	# this error is in ms
	frame['tau_w_error'] = frame[['tau_w', 'time_res', 'min_sigma', 'max_sigma', 'min_sigma_error', 'max_sigma_error', 'angle', 'angle_error']].apply(
		tauwerror_expr,
		axis=1
	)

	frame['sigma_t_ms'] = frame['sigma_t']*1e3
	frame['tau_w_ms'] = frame['tau_w']*1e3

	## Redshift corrections
	frame['drift_z'] = frame[['drift_over_nuobs', 'z']].apply(lambda row: row['drift_over_nuobs']*(1+row['z']), axis=1)
	frame['tau_w_ms_z'] = frame[['tau_w_ms', 'z']].apply(lambda row: row['tau_w_ms']/(1+row['z']), axis=1)

	return frame

def cleanAngle(row):
	angle = row['angle']
	if angle < 0 or angle > np.pi:
		if angle > np.pi:
			return angle % (np.pi)
		elif angle < 0 and angle > -np.pi:
			return angle + np.pi
		elif angle < 0 and angle < -np.pi:
			angle = angle % (2*np.pi)
			if angle > np.pi:
				return angle - np.pi
			else:
				return angle
	else:
		return angle

def _plotAnglevsDM(frames, annotate=False, save=False, drops=[]):
	thetamodel = scipy.odr.Model(atanmodel)
	offsetmodel = scipy.odr.Model(offset_atanmodel)

	for frame in frames:
		frame = computeModelDetails(frame)
		frame['angle_clean'] = frame[['angle']].apply(cleanAngle, axis=1) - (np.pi/2)

	def errorexpr(frame):
		ex = frame['tau_w_error']
		ey = frame['angle_error']
		return ex, ey

	markersize = 125 #100
	fontsize = 25 #18
	annotsize = 14
	logscale = False
	figsize = (15, 8)
	ax = frames[0].drop(drops).plot.scatter(x='tau_w_ms', y='angle_clean',
							   xerr=errorexpr(frame[0])[0],
							   yerr=errorexpr(frame[0])[0],
							   figsize=figsize, s=markersize, c='b', colorbar=False, fontsize=fontsize, logy=logscale, logx=logscale, marker='X', edgecolors='k',
							   label='$\\Delta$DM = 1/2 pc/cm$^3$')

	markers = ['o', 'p', 's']
	for frame, c, label, mark in zip(frames[:3], ['r', 'c', 'g'], ['$\\Delta$DM = 0 pc/cm$^3$', '$\\Delta$DM = -1 pc/cm$^3$', '$\\Delta$DM = -2 pc/cm$^3$'], markers):
		frame.drop(drops).plot.scatter(ax=ax, x='tau_w_ms', y='angle_clean',
									xerr=errorexpr(frame)[0],
									yerr=errorexpr(frame)[1],
									figsize=figsize, s=markersize, c=c, colorbar=False, fontsize=fontsize, logy=logscale, logx=logscale, marker=mark, edgecolors='k',
									label=label)

	ax.set_xlabel('Sub-burst Duration $t_\\mathrm{w}$ (ms)', size=fontsize)
	#ax.set_ylabel('-$\pi/2 + $ Gaussian2d angle (rad)', size=fontsize)
	ax.set_ylabel('Sub-burst Drift Angle $\\theta$ (rad)', size=fontsize)

	## Find Fits
	lstyles = ['-', '--', '-.', ':']
	for frame, drops, pcol, beta, lstyle in zip(frames, [[15], [15], [15], [15]], ['r', 'c', 'g', 'b'], [-6, -4, -3, -9], lstyles):
		if frame.equals(frames[0]):
			model = thetamodel
		else:
			model = offsetmodel
		#model = thetamodel

		datafitter = scipy.odr.RealData(frame.drop(drop)['tau_w_ms'],
								frame.drop(drop)['angle_clean'],
								sx=errorexpr(frame)[0],
								sy=errorexpr(frame)[1])
		anglefitter = scipy.odr.ODR(datafitter, model, beta0=[1])
		anglefitter.set_job(fit_type=0)
		anglefit = anglefitter.run()

		tws = np.linspace(0, 8.5, num=80)
		print(anglefit.beta)
		#print(anglefit.beta[0])
		if model == thetamodel:
			plt.plot(tws, np.arctan(tws/anglefit.beta[0]), c=pcol, label="$\\tan^{{-1}}(t_\\mathrm{{w}}/{:.2f})$".format(anglefit.beta[0]), linestyle=lstyle)
		elif model == offsetmodel:
			plt.plot(tws, np.arctan(tws/zero_ddm_fit) + anglefit.beta[0], c=pcol, label="$\\tan^{{-1}}(t_\\mathrm{{w}}/{:.2f}) {:+.2f}$ rad".format(zero_ddm_fit, anglefit.beta[0]), linestyle=lstyle)

	## Point Annotations
	if annotate:
		for k, v in frames[0].iterrows():
			ax.annotate(int(k) if k != 15.5 else k, (v['tau_w_ms'], v['angle_clean']), xytext=(-3,5), textcoords='offset points', weight='bold', size=annotsize)

	ax.set_xlim(0, 8.5)
	plt.title("Fit Angles for FRB180916 at different DMs", size=25)
	plt.legend(fontsize='xx-large')
	if save:
		for fformat in ['png', 'pdf', 'eps']: plt.savefig('angleatdifferentDMs.{}'.format(fformat))

def atanmodel(B, x):
	return np.arctan(x/B[0])

def offset_atanmodel(B, x, zero_ddm_fit=6.554):
	return np.arctan(x/zero_ddm_fit) + B[0]

def reciprocal(x, a):
	return a/x

def reciprocal_log(x, b):
	return -x+b

def log_log(x, k, b):
	return k*x+b

def reciprocal_odr(B, x):
	return B[0]/x

def reciprocal_odr_log(B, x):
	return -x+B[0]

def fitreciprocal(x, data, sigma=1):
	guess = [522]
	abs_sigma = True
	if (type(sigma) == int) and (sigma == 1):
		abs_sigma = False
	sigma = np.zeros(len(data.ravel())) + sigma

	sigma = np.zeros(len(data.ravel())) + sigma
	popt, pcov = scipy.optimize.curve_fit(reciprocal, x, data, p0=guess, sigma=sigma, absolute_sigma=abs_sigma)
	return popt, pcov

def fitreciprocal_log(x, data, sigma=1, loglog=False):
	guess = [522]
	abs_sigma = True
	if (type(sigma) == int) and (sigma == 1):
		abs_sigma = False
	sigma = np.zeros(len(data.ravel())) + sigma

	if loglog:
		guess = [1,1]
		popt, pcov = scipy.optimize.curve_fit(log_log, x, data, p0=guess, sigma=sigma, absolute_sigma=abs_sigma)
	else:
		popt, pcov = scipy.optimize.curve_fit(reciprocal_log, x, data, p0=guess, sigma=sigma, absolute_sigma=abs_sigma)
	return popt, pcov
