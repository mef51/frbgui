import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.odr

def computeModelDetails(frame):
	""" Takes a dataframe and computes columns related to the dynamical frb model """

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
	if 'z' in frame.index:
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

def _plotDriftVsDuration(frames=[], labels=[]):
	""" wip """
	plt.rcParams["errorbar.capsize"] = 4
	plt.rcParams["font.family"] = "serif"

	markersize = 125#100
	fontsize = 25 #18
	annotsize = 14
	logscale = True
	filename = 'log_drift_over_nu_obsvsduration' if logscale else 'drift_over_nu_obsvsduration'
	figsize = (17, 8)
	figsize = (17, 9)
	# figsize = (14, 10)

	def errorexpr(frame):
		ex = np.sqrt(frame['red_chisq'])*frame['tau_w_error']
		ey = np.sqrt(frame['red_chisq'])*frame['drift error (mhz/ms)']/frame['center_f']
		return ex, ey

	ax = frames[0].plot.scatter(x='tau_w_ms', y='drift_over_nuobs',
			xerr=errorexpr(frames[0])[0],
			yerr=errorexpr(frames[0])[1],
			figsize=figsize, s=markersize, c='color', colorbar=False, fontsize=fontsize,
			logy=logscale, logx=logscale, marker='o', edgecolors='k',
			label=labels[0])

	for frame, lbl in zip(frames[1:], labels[1:]):
		frame.plot.scatter(ax=ax, x='tau_w_ms', y='drift_over_nuobs',
			xerr=errorexpr(frame)[0],
			yerr=errorexpr(frame)[1],
			figsize=figsize, s=markersize, c='color', colorbar=False, fontsize=fontsize,
			logy=logscale, logx=logscale, marker='p', edgecolors='k',
			label=lbl)

	return
	######

	ax = michillibursts.plot.scatter(x='tau_w_ms', y='drift_over_nuobs',
								   xerr=np.sqrt(michillibursts['red_chisq'])*michillibursts['tau_w_error'],
								   yerr=np.sqrt(michillibursts['red_chisq'])*michillibursts['drift error (mhz/ms)']/michillibursts['center_f'],
								   figsize=figsize, s=markersize, c='color', colorbar=False, fontsize=fontsize, logy=logscale, logx=logscale, marker='o', edgecolors='k',
								   label='FRB121102 Michilli et al. (2018)')

	selectbursts180916.plot.scatter(ax=ax, x='tau_w_ms', y='drift_over_nuobs',
								   xerr=np.sqrt(selectbursts180916['red_chisq'])*selectbursts180916['tau_w_error'],
								   yerr=np.sqrt(selectbursts180916['red_chisq'])*selectbursts180916['drift error (mhz/ms)']/selectbursts180916['center_f'],
								   figsize=figsize, s=markersize, c='color', colorbar=False, fontsize=fontsize, logy=logscale, logx=logscale, marker='p', edgecolors='k',
								   label='FRB180916.J0158+65 bursts')
	selectbursts180814.plot.scatter(ax=ax, x='tau_w_ms', y='drift_over_nuobs',
								   xerr=np.sqrt(selectbursts180814['red_chisq'])*selectbursts180814['tau_w_error'],
								   yerr=np.sqrt(selectbursts180814['red_chisq'])*selectbursts180814['drift error (mhz/ms)']/selectbursts180814['center_f'],
								   figsize=figsize, s=markersize+50, c='color', colorbar=False, fontsize=fontsize, logy=logscale, logx=logscale, marker='X', edgecolors='k',
								   label='FRB180814.J0422+73 bursts')
								   #label='FRB180814.J0422+73 bursts @ DM={} pc/cm$^3$'.format(dms180814[dm_idx]))
	burstsSGR1935.plot.scatter(ax=ax, x='tau_w_ms', y='drift_over_nuobs',
								   xerr=np.sqrt(burstsSGR1935['red_chisq'])*burstsSGR1935['tau_w_error'],
								   yerr=np.sqrt(burstsSGR1935['red_chisq'])*burstsSGR1935['drift error (mhz/ms)']/burstsSGR1935['center_f'],
								   figsize=figsize, s=markersize, c='color', colorbar=False, fontsize=fontsize, logy=logscale, logx=logscale, marker='p', edgecolors='k',
								   label='SGR1935+2154 bursts')

	otherbursts.head(5).plot.scatter(ax=ax, x='tau_w_ms', y='drift_over_nuobs',
							 xerr=np.sqrt(otherbursts['red_chisq'])*otherbursts['tau_w_error'],
							 yerr=np.sqrt(otherbursts['red_chisq'])*otherbursts['drift error (mhz/ms)']/otherbursts['center_f'], edgecolors='k',
							 figsize=figsize, s=markersize, c='color', colorbar=False, marker='d', label='FRB121102 Gajjar et al. (2018)')
	otherbursts.tail(1).plot.scatter(ax=ax, x='tau_w_ms', y='drift_over_nuobs',
							 xerr=np.sqrt(otherbursts['red_chisq'])*otherbursts['tau_w_error'],
							 yerr=np.sqrt(otherbursts['red_chisq'])*otherbursts['drift error (mhz/ms)']/otherbursts['center_f'], edgecolors='k',
							 figsize=figsize, s=markersize, c='color', colorbar=False, marker='s', label='FRB121102 Josephy et al. (2019)')

	# for k, v in otherbursts.iterrows():
	#     ax.annotate(k+'_c', (v['tau_w_ms'], v['drift_over_nuobs']), xytext=(-3,5), textcoords='offset points', weight='bold', size=annotsize)
	# for k, v in michillibursts.iterrows():
	#     ax.annotate(k, (v['tau_w_ms'], v['drift_over_nuobs']), xytext=(-3,5), textcoords='offset points', weight='bold', size=annotsize)
	# for k, v in selectbursts180916.iterrows():
	#     ax.annotate(int(k) if k != 15.5 else k, (v['tau_w_ms'], v['drift_over_nuobs']), xytext=(-3,5), textcoords='offset points', weight='bold', size=annotsize)
	for k, v in burstsSGR1935.iterrows():
		if v['drift_over_nuobs'] > 0 or not logscale:
			ax.annotate(k, (v['tau_w_ms'], v['drift_over_nuobs']), xytext=(-3,5), textcoords='offset points', weight='bold', size=annotsize)
	# for k, v in selectbursts180814.iterrows():
	#     if v['drift_over_nuobs'] > 0:
	#         ax.annotate(k, (v['tau_w_ms'], v['drift_over_nuobs']), xytext=(-3,5), textcoords='offset points', weight='bold', size=annotsize)

	if not logscale:
		pass
		# ax.set_xlim(-0.1, 20)
		# ax.set_ylim(-0.2, 4)
	elif logscale:
		ax.set_xlim(0.04, 20)
		ax.set_ylim(10**-3, 10**1)

	# ax.set_title('Sub-burst Drift Rate vs. Burst Duration (fit to Michilli bursts)', size=fontsize)
	ax.set_xlabel('Sub-burst Duration $t_\\mathrm{w}$ (ms)', size=fontsize)
	ax.set_ylabel('Sub-burst Drift Rate $\,(1/\\nu_{\\mathrm{obs}}) \left|\\frac{d\\nu_\\mathrm{obs}}{dt_\\mathrm{D}}\\right|$ (ms$^{-1}$)', size=fontsize)

	def driftnu_error(frame):
		sx = np.log((frame['tau_w_ms'] + np.sqrt(frame['red_chisq'])*frame['tau_w_error']) / frame['tau_w_ms'])
		sy = np.log((frame['drift_over_nuobs'] + np.sqrt(frame['red_chisq'])*(frame['drift error (mhz/ms)'])) / frame['drift_over_nuobs'])
		return sx, sy

	# ODR fit log
	num_to_fit = 24 #23 to exlude chime
	fitdata_log = scipy.odr.RealData(np.log(selectbursts121102.head(num_to_fit)['tau_w_ms']),
								 np.log(selectbursts121102.head(num_to_fit)['drift_over_nuobs']),
								 sx=driftnu_error(selectbursts121102.head(num_to_fit))[0],
								 sy=driftnu_error(selectbursts121102.head(num_to_fit))[1])
								 #sx=np.log(np.sqrt(selectbursts121102.head(num_to_fit)['red_chisq'])*selectbursts121102.head(num_to_fit)['tau_w_error']),
								 #sy=np.log(np.sqrt(selectbursts121102.head(num_to_fit)['red_chisq'])*selectbursts121102.head(num_to_fit)['drift error (mhz/ms)']/selectbursts121102.head(num_to_fit)['center_f']))
	odrfitter_log = scipy.odr.ODR(fitdata_log, fit_model_log, beta0=[500])
	odrfitter_log.set_job(fit_type=0)
	odrfit_log = odrfitter_log.run()

	# ODR fit log 180916
	fitdata_log_180916 = scipy.odr.RealData(np.log(selectbursts180916['tau_w_ms']),
								 np.log(selectbursts180916['drift_over_nuobs']),
								 sx=driftnu_error(selectbursts180916)[0],
								 sy=driftnu_error(selectbursts180916)[1])
								 #sx=np.log(np.sqrt(selectbursts180916['red_chisq'])*selectbursts180916['tau_w_error']),
								 #sy=np.log(np.sqrt(selectbursts180916['red_chisq'])*selectbursts180916['drift error (mhz/ms)']/selectbursts180916['center_f'] ))
	odrfitter_log180916 = scipy.odr.ODR(fitdata_log_180916, fit_model_log, beta0=[1000])
	odrfitter_log180916.set_job(fit_type=0)
	odrfit_log180916 = odrfitter_log180916.run()

	# ODR fit log 180814
	fitdata_log_180814 = scipy.odr.RealData(np.log(selectbursts180814['tau_w_ms']),
								 np.log(selectbursts180814['drift_over_nuobs']),
								 sx=driftnu_error(selectbursts180814)[0],
								 sy=driftnu_error(selectbursts180814)[1])
								 #sx=np.log(np.sqrt(selectbursts180814['red_chisq'])*selectbursts180814['tau_w_error']),
								 #sy=np.log(np.sqrt(selectbursts180814['red_chisq'])*selectbursts180814['drift error (mhz/ms)']/selectbursts180814['center_f'] ))
	odrfitter_log180814 = scipy.odr.ODR(fitdata_log_180814, fit_model_log, beta0=[1000])
	odrfitter_log180814.set_job(fit_type=0)
	odrfit_log180814 = odrfitter_log180814.run()

	### Plot fits
	x = np.linspace(0, 20, num=1200)
	opts  = [np.exp(odrfit_log.beta), np.exp(odrfit_log180916.beta), np.exp(odrfit_log180814.beta)]
	errs  = [opts[0]*(np.exp(odrfit_log.sd_beta)-1), opts[1]*(np.exp(odrfit_log180916.sd_beta)-1), opts[2]*(np.exp(odrfit_log180814.sd_beta)-1)]
	#errs  = [np.exp(odrfit_log.sd_beta), np.exp(odrfit_log180916.sd_beta), np.exp(odrfit_log180814.sd_beta)]
	#names = ['Fit from FRB121102 bursts', 'Fit from FRB180916.J0158+65 bursts', 'Fit from FRB180814.J0422+73 bursts']

	names = ['FRB121102 fit', 'FRB180916.J0158+65 fit', 'FRB180814.J0422+73 fit']
	ls    = ['r-', 'b--', 'g-.']
	for opt, err, name, l in zip(opts, errs, names, ls):
		print(opt, err)
		lstr = '{} ({:.3f} $\pm$ {:.3f}) / $t_w$'.format(name, opt[0], err[0])
		plt.plot(x, opt[0]/x, l, label=lstr)

	handles, labels = ax.get_legend_handles_labels()
	#handles = [handles[0], handles[2], handles[4], handles[5], handles[1], handles[3]]
	#labels = [labels[0], labels[2], labels[4], labels[5], labels[1], labels[3]]
	# handles = [handles[2], handles[5], handles[6], handles[4], handles[3], handles[0], handles[1]]
	# labels = [labels[2], labels[5], labels[6], labels[4], labels[3], labels[0], labels[1]]
	handles = [handles[3], handles[6], handles[7], handles[5], handles[4], handles[0], handles[1], handles[2]]
	labels = [labels[3], labels[6], labels[7], labels[5], labels[4], labels[0], labels[1], labels[2]]
	plt.legend(handles, labels, fontsize='xx-large')
	# 2 digits, ytitle

	# plt.title("Non Redshift Corrected", size=20)
	#plt.title("FRB121102, FRB180916.J0158+65, and FRB180814.J0422+73", size=20)
	plt.tight_layout()
	# plt.savefig('180814_dm{}.png'.format(dms180814[dm_idx]))
	# print('180814_dm{}.png'.format(dms180814[dm_idx]))
	# for f in ['png', 'pdf', 'eps']: plt.savefig('figures/{}.{}'.format(filename, f))
	# for f in ['png']: plt.savefig('figures/{}.{}'.format(filename, f))
