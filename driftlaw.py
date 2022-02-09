import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.odr
import itertools, re

def computeModelDetails(frame, channelSpaceDuration=False):
	""" Takes a dataframe and computes columns related to the dynamical frb model """
	frame = frame.copy()
	tauwerror_expr = lambda r: 1e3*r['time_res (s)']*np.sqrt(r['max_sigma']**6*r['min_sigma_error']**2*np.cos(r['theta']-np.pi/2)**4 + r['angle_error']**2*r['max_sigma']**2*r['min_sigma']**2*(-r['max_sigma']**2 + r['min_sigma']**2)**2*np.cos(r['theta']-np.pi/2)**2*np.sin(r['theta']-np.pi/2)**2 + r['max_sigma_error']**2*r['min_sigma']**6*np.sin(r['theta']-np.pi/2)**4)/(r['max_sigma']**2*np.cos(r['theta']-np.pi/2)**2 + r['min_sigma']**2*np.sin(r['theta']-np.pi/2)**2)**1.5

	frame['slope_abs'] = -1*(frame['slope (mhz/ms)']) # multiply be negative 1 because of later measurement exclusions
	frame['slope_over_nuobs'] = frame[['slope_abs','center_f']].apply(lambda row: row['slope_abs'] / row['center_f'], axis=1)
	frame['slope_over_nuobs_err'] = np.sqrt(frame['red_chisq'])*frame['slope error (mhz/ms)']/frame['center_f']
	frame['recip_slope_over_nuobs'] = 1/frame['slope_over_nuobs']
	frame['slope_abs_nuobssq'] = frame['slope_abs']/frame['center_f']**2/1000 # unitless
	frame['min_sigma'] = frame[['sigmax','sigmay']].apply(lambda row: min(abs(row['sigmax']), abs(row['sigmay'])), axis=1)
	frame['max_sigma'] = frame[['sigmax','sigmay']].apply(lambda row: max(abs(row['sigmax']), abs(row['sigmay'])), axis=1)
	# the following two lines assume that if sigmax > sigmay, then sigmax_error > sigmay_error, which is true (so far) for this dataset
	frame['min_sigma_error'] = frame[['sigmax_error','sigmay_error']].apply(lambda row: min(row['sigmax_error'], row['sigmay_error']), axis=1)
	frame['max_sigma_error'] = frame[['sigmax_error','sigmay_error']].apply(lambda row: max(row['sigmax_error'], row['sigmay_error']), axis=1)

	frame['sigma_t']   = frame[['min_sigma','time_res (s)']].apply(lambda row: row['min_sigma']*row['time_res (s)'], axis=1)

	frame['tau_w_old'] = frame[['time_res (s)', 'min_sigma', 'max_sigma', 'theta']].apply(
		lambda r: r['time_res (s)']*r['min_sigma']*r['max_sigma'] / np.sqrt( np.abs((np.sin(r['theta']-np.pi/2)*r['min_sigma'])**2 + (np.cos(r['theta']-np.pi/2)*r['max_sigma'])**2 )),
		axis=1
	)

	frame['tau_w'] = frame[['min_sigma', 'max_sigma', 'theta']].apply(
		lambda r: r['min_sigma']*r['max_sigma'] / np.sqrt( np.abs((np.sin(r['theta']-np.pi/2)*r['min_sigma'])**2 + (np.cos(r['theta']-np.pi/2)*r['max_sigma'])**2 )),
		axis=1
	)
	frame['t_from_sigma'] = frame['min_sigma']*np.sin(frame['theta'])

	# this error is in ms
	frame['tau_w_error'] = frame[['tau_w', 'time_res (s)', 'min_sigma', 'max_sigma', 'min_sigma_error', 'max_sigma_error', 'theta', 'angle_error']].apply(
		tauwerror_expr,
		axis=1
	)

	frame['sigma_t_ms'] = frame['sigma_t']*1e3
	frame['tau_w_ms'] = frame['tau_w'] # alias, units are implicit when solving in data space
	frame['tau_w_ms_old'] = frame['tau_w_old']*1e3
	if channelSpaceDuration: # hack to force channel space equation for old measurements
		frame['tau_w_ms'] = frame['tau_w_ms_old']

	## Redshift corrections
	if 'z' in frame.index:
		frame['slope_z'] = frame[['slope_over_nuobs', 'z']].apply(lambda row: row['slope_over_nuobs']*(1+row['z']), axis=1)
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

def modelerror(frame):
	ex = np.sqrt(frame['red_chisq'])*frame['tau_w_error']
	ey = np.sqrt(frame['red_chisq'])*frame['slope error (mhz/ms)']/frame['center_f']
	return ex, ey

def rangeerror(frame):
	"""
	These ranges are not errors in the statistical sense. they are the min/max values, which should
	be larger than the real errors. So this is extremely conservative while also being easier
	to compute.

	The strange shape of the returned value is due to a quirk in the way pandas handles asymmetric
	errors.
	"""
	ex = [np.array([frame['tau_w_ms'] - frame['tw_min'], frame['tw_max'] - frame['tau_w_ms']])]
	ey = [np.array([frame['slope_over_nuobs'] - frame['slope_nu_min'], frame['slope_nu_max'] - frame['slope_over_nuobs']])]
	return ex, ey

def log_error(frame):
	""" see modelerror() """
	sx = np.log((frame['tau_w_ms'] + np.sqrt(frame['red_chisq'])*frame['tau_w_error']) / frame['tau_w_ms'])
	sy = np.log((frame['slope_over_nuobs'] + np.sqrt(frame['red_chisq'])*(frame['slope error (mhz/ms)'])) / frame['slope_over_nuobs'])
	return sx, sy

def rangelog_error(frame):
	""" The range errors are asymmetric. Average the error """
	ex, ey = rangeerror(frame)
	ex = np.log((frame['tau_w_ms'] + (ex[0][0]+ex[0][1])/2 ) / frame['tau_w_ms'])
	ey = np.log((frame['slope_over_nuobs'] + (ey[0][0]+ey[0][1])/2) / frame['slope_over_nuobs'])
	return ey, ey
	# return np.log(np.maximum(ex[0][0], ex[0][1])), np.log(np.maximum(ey[0][0], ey[0][1]))

def rangeerror_odr(frame):
	""" The range errors are asymmetric. Take the largest error """
	ex, ey = rangeerror(frame)
	return np.maximum(ex[0][0], ex[0][1]), np.maximum(ey[0][0], ey[0][1])

def fitodr(frame, beta0=[1000], errorfunc=log_error, log=True):
	fit_model = scipy.odr.Model(reciprocal_odr)
	fit_model_log = scipy.odr.Model(reciprocal_odr_log)

	fitdata = scipy.odr.RealData(frame['tau_w_ms'],
								 frame['slope_over_nuobs'],
								 sx=rangeerror_odr(frame)[0],
								 sy=rangeerror_odr(frame)[1])
	fitdata_log = scipy.odr.RealData(np.log(frame['tau_w_ms']),
									 np.log(frame['slope_over_nuobs']),
									 sx=errorfunc(frame)[0],
									 sy=errorfunc(frame)[1])

	odrfitter_log = scipy.odr.ODR(fitdata_log, fit_model_log, beta0=beta0)
	odrfitter_log.set_job(fit_type=0)

	odrfitter = scipy.odr.ODR(fitdata, fit_model, beta0=beta0)
	odrfitter.set_job(fit_type=0)

	if log:
		# print('log odr')
		return odrfitter_log.run()
	else:
		# print('linear odr')
		return odrfitter.run()

def sloperanges(source):
	"""
	Given all burst and model data at different trial DMs,
	computes the range of slopes durations across the range of trial DMs
	"""

	yaxis = 'slope_over_nuobs'
	xaxis ='tau_w_ms'
	for burst in source.index.unique():
		burstdf = source.loc[burst]
		eduration   = np.sqrt(burstdf['red_chisq'])*burstdf['tau_w_error']
		eslopenuobs = np.sqrt(burstdf['red_chisq'])*burstdf['slope error (mhz/ms)']/burstdf['center_f']

		dmax, dmin = np.max(burstdf[yaxis] + eslopenuobs), np.min(burstdf[yaxis] - eslopenuobs)
		tmax, tmin = np.max(burstdf[xaxis] + eduration)  , np.min(burstdf[xaxis] - eduration)

		source.loc[burst, 'slope_nu_max'] = dmax
		source.loc[burst, 'slope_nu_min'] = dmin
		source.loc[burst, 'slope_max'] = dmax*burstdf['center_f']
		source.loc[burst, 'slope_min'] = dmin*burstdf['center_f']
		source.loc[burst, 'tw_max']    = tmax
		source.loc[burst, 'tw_min']    = tmin
		source.loc[burst, 'slope_nu_plus']  = burstdf[yaxis] + eslopenuobs
		source.loc[burst, 'slope_nu_minus'] = burstdf[yaxis] - eslopenuobs
		source.loc[burst, 'tw_plus']  = burstdf[xaxis] + eduration
		source.loc[burst, 'tw_minus'] = burstdf[xaxis] - eduration

		# print(f'burst: {burst},\t\tsloperange = ({dmin}, {dmax}),\t\ttwrange = ({tmin}, {tmax})')

	return source

def bakeMeasurements(sources, names, exclusions, targetDMs, logging=True,
					 tagColors=['r', 'r', 'b', 'g', 'y', 'c', 'tomato', 'c'],
					 showDMtraces=False):
	"""
	1. filter
	2. color
	3. tags data that is at the target dm and splits up measurements by dm
	4. finds range of fits by source
	5. do some extra computations (?)
	"""
	if not (len(sources) == len(names) == len(exclusions) == len(targetDMs)):
		raise ValueError(
			"Sources, names, exclusions, and targetDMs are not the same length",
			len(sources), len(names), len(exclusions), len(targetDMs)
			)
	def exclusionLogger(exclusionstr, frame):
		if not logging:
			return
		print(exclusionstr)
		count = 0
		if frame.empty:
			print(">> bursts excluded: \n\tnone")
		else:
			print('>> bursts excluded:',
				  *[f"\n\t{exname} DM = {row['DM']}" for exname, row in frame.iterrows()])
			count += len(frame)
		print(f'>> # of measurements excluded: {count}')
		print(f'>> # of bursts limited: {len(frame.index.unique())}')
		return count

	def DMrangeLogger(sname, frame, burstwise=False):
		if not logging:
			return
		if not burstwise:
			if type(frame) == pd.core.series.Series:
				print(f"{sname}\tDM Range: {frame[0]} pc/cm3")
			else:
				print(f"{sname}\t DM Range: {min(frame['DM'])} - {max(frame['DM'])} pc/cm3")
		else:
			print(f'>>> {sname} DM Range by burst:')
			for bname in frame.index.unique():
				DMrangeLogger(f'\tburst {bname}:', frame.loc[bname])
		return

	bakeddata = { 'frames': [], 'labels': [], 'colors': itertools.cycle(tagColors) }
	fitdata = pd.DataFrame()
	extradata = {}
	extradata['Bconstants'] = []
	for source, name, exclude, targetDM in zip(sources, names, exclusions, targetDMs):
		fitframes, fitlabels = [], []
		print('\n#', name + ':')
		DMrangeLogger('', source)
		if logging: print('> manually excluded:', exclude)
		source = source.drop(exclude)

		## Burst Exclusions
		# exclude bursts that where a fit was not found
		exclusionLogger('> exclusion rule: no fit found', source[~(source.amplitude > 0)])
		source = source[source.amplitude > 0]

		# exclude positive drifts, we assume they are non-physical and an artifact of over dedispersion
		exclusionLogger('> exclusion rule: slope_abs < 0', source[~(source.slope_abs > 0)])
		source = source[source.slope_abs > 0]

		# exclude slopes with large relative errors
		relerrthreshold = 0.4 # (40%)
		exclusionLogger(f'> exclusion rule: rel slope error > {relerrthreshold}', source[~(abs(source['slope error (mhz/ms)']/source['slope (mhz/ms)']) < relerrthreshold)])
		source = source[abs(source['slope error (mhz/ms)']/source['slope (mhz/ms)']) < relerrthreshold]

		# exclude slopes with huge errors, they are vertical and poorly measured
		errorthreshold = 1e8
		exclusionLogger('> exclusion rule: slope error > {}'.format(errorthreshold), source[~(source['slope error (mhz/ms)'] < errorthreshold)])
		source = source[source['slope error (mhz/ms)'] < errorthreshold]

		# exclude slope ranges that flip signs over the dm range (this is analgous to the large error exclusion)
		exclusionLogger('> exclusion rule: negative slope ranges', source[~(source['slope_over_nuobs'] > source['slope_over_nuobs_err'])])
		source = source[source['slope_over_nuobs'] > source['slope_over_nuobs_err']]

		DMrangeLogger(name, source, burstwise=True) # log the dm range after we're done excluding
		source = sloperanges(source) # compute drift ranges after burst exclusions

		for dm, color in zip(source.DM.unique(), itertools.cycle(['r', 'y', 'b', 'k', 'g', 'c', 'm'])):
			df = source[source.DM == dm]
			df['color'] = color
			fitframes.append(df)
			fitlabels.append(dm)

			# Figure 1
			if np.isclose(dm, targetDM):
				print(f'>> num bursts remaining = {len(df.index.unique())}')
				df['color'] = next(bakeddata['colors'])
				bakeddata['frames'].append(df)
				bakeddata['labels'].append(name)

			# For Ref figure A, turn off otherwise
			# if dm == 565:
			#     df['color'] = 'b'
			# if dm == 348.82:
			#     df['color'] = 'r'

			# Figure 5: compute drift/nu^2 vs nu range of fits
			nus = np.linspace(0, 20000, num=2000)
			popt, pcov = scipy.optimize.curve_fit(lambda x, c: c, nus, df['slope_abs_nuobssq'])
			extradata['Bconstants'].append((popt[0], np.sqrt(pcov[0])))

		markcolors = []
		for f in fitframes:
			markcolors.append(f['color'].iloc[0])
		markcolors = itertools.cycle(markcolors)
		# manylines = ['r--', 'y-', 'b-.', 'k--', 'g--', 'c-.', 'r-', 'b-', 'k-', 'g-', 'c-']
		linestyles = ['--', '-', '-.', '--', '--', '-.', '-', '-', '-', '-', '-']
		manylines = [next(markcolors)+lst for lst in linestyles]
		labels = [round(lbl, 2) for lbl in fitlabels]
		tempax, fits = plotSlopeVsDuration(fitframes, labels, annotatei=[], fitlines=manylines,
			logscale=True, fiterrorfunc=log_error, dmtrace=True, hidefitlabel=True)

		for fit in fits:
			fitdata = fitdata.append({'name': name.split(' DM')[0],
										'DM': fit[0],
										'param': fit[1],
										'err': fit[2] ,
										'red_chisq': fit[3],
										'numbursts': fit[5]},
										ignore_index=True)
	if not showDMtraces: plt.close('all')
	return bakeddata, fitdata, sources, extradata

def listDMFitResults(fitdata):
	for name in fitdata.name.unique():
		df = fitdata.loc[(fitdata.name == name)]# & (fitdata.numbursts > 1.0)]
		df = df.set_index('DM')
		display(df)
		print('minimal red_chisq:')
		display(df[df.red_chisq == df.red_chisq.min()])

def plotSlopeVsDuration(frames=[], labels=[], title=None, logscale=True, annotatei=0,
						markers=['o', 'p', 'X', 'd', 's'], hidefit=[], hidefitlabel=False,
						fitlines=['r-', 'b--', 'g-.'], fitextents=None,
						errorfunc=modelerror, fiterrorfunc=rangelog_error, dmtrace=False):
	""" wip """
	plt.rcParams["errorbar.capsize"] = 4
	plt.rcParams["font.family"] = "serif"

	markersize = 125#100
	fontsize = 25 #18
	annotsize = 14
	filename = 'log_slope_over_nu_obsvsduration' if logscale else 'slope_over_nu_obsvsduration'
	figsize = (17, 8)
	figsize = (17, 9)
	# figsize = (14, 10)

	yaxis = 'slope_over_nuobs'
	yaxis_lbl = 'Sub-burst Slope $\\,\\left|\\frac{d\\nu_\\mathrm{obs}}{dt_\\mathrm{D}}\\right|(1/\\nu_{\\mathrm{obs}})$ (ms$^{-1}$)'
	# yaxis = 'recip_slope_over_nuobs'
	# yaxis_lbl = 'nu_obs / slope'

	if type(markers) == list:
		markers = itertools.cycle(markers)
	if type(fitlines) == list:
		fitlines = itertools.cycle(fitlines)

	ax = frames[0].plot.scatter(x='tau_w_ms', y=yaxis,
			xerr=errorfunc(frames[0])[0],
			yerr=errorfunc(frames[0])[1],
			figsize=figsize, s=markersize, c='color', colorbar=False, fontsize=fontsize,
			logy=logscale, logx=logscale, marker=next(markers), edgecolors='k',
			label=labels[0])

	for frame, lbl in zip(frames[1:], labels[1:]):
		frame.plot.scatter(ax=ax, x='tau_w_ms', y=yaxis,
			xerr=errorfunc(frame)[0],
			yerr=errorfunc(frame)[1],
			figsize=figsize, s=markersize, c='color', colorbar=False, fontsize=fontsize,
			logy=logscale, logx=logscale, marker=next(markers), edgecolors='k',
			label=lbl)

	if type(annotatei) == int:
		annotatei =[annotatei]
	for ai in annotatei:
		if ai < len(frames):
			for k, v in frames[ai].iterrows():
				if v[yaxis] > 0 or not logscale:
					ax.annotate(k, (v['tau_w_ms'], v[yaxis]), xytext=(-3,5),
						textcoords='offset points', weight='bold', size=annotsize)

	alldata = pd.concat([f for f in frames])
	if not fitextents:
		fitextents = min(alldata['tau_w_ms'])*0.9, max(alldata['tau_w_ms'])*1.1

	logfit = True
	if type(hidefit) == int:
		hidefit = [hidefit]

	fits = []
	for fi, (frame, label, line) in enumerate(zip(frames, labels, fitlines)):
		x = np.linspace(fitextents[0], fitextents[1], num=1200)
		if logfit:
			fit = fitodr(frame, errorfunc=fiterrorfunc)
			param, err = np.exp(fit.beta[0]), np.exp(fit.beta[0])*(np.exp(fit.sd_beta[0])-1)
		else:
			fit = fitodr(frame, log=logfit)
			param, err = fit.beta[0], fit.sd_beta[0]

		## compute reduced chisq
		# parameter error
		ex = frame['tau_w_error']*np.sqrt(frame['red_chisq'])
		ey = frame['slope error (mhz/ms)']/frame['center_f']*np.sqrt(frame['red_chisq'])
		data_err = np.sqrt(ex**2 + ey**2)
		residuals = frame['slope_over_nuobs'] - param/frame['tau_w_ms']
		chisq = np.sum((residuals / data_err) ** 2)
		red_chisq = chisq / (len(frame) - 1)
		# print(residuals)
		fits.append([label, param, err, red_chisq, residuals, len(frame)])

		lstr = ''
		if not hidefitlabel:
			lstr = '{} fit ({:.3f} $\\pm$ {:.3f}) $t_w^{{-1}}$'.format(label, param, err)
		if fi not in hidefit:
			color = re.match(r"\w+", line)[0]
			lstyle = line.split(color)[-1]
			plt.plot(x, param/x, color=color, ls=lstyle, label=lstr)

	if title:
		ax.set_title(title, size=fontsize)

	if dmtrace:
		sorteddata = pd.concat([frames[dmi] for dmi in np.argsort(labels)])
		for bid in sorteddata.index.unique():
			plt.plot(sorteddata.loc[bid]['tau_w_ms'], sorteddata.loc[bid]['slope_over_nuobs'])

	ax.set_xlabel('Sub-burst Duration $t_\\mathrm{w}$ (ms)', size=fontsize)
	ax.set_ylabel(yaxis_lbl, size=fontsize)
	plt.legend(fontsize='xx-large')
	# plt.legend()
	plt.tight_layout()

	return ax, fits

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
	ax.set_ylabel('Sub-burst slope Angle $\\theta$ (rad)', size=fontsize)

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
