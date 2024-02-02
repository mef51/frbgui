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
	sigma = np.sqrt(abs(sum(data*(x-xo)**2)/sum(data)))
	guess = [np.max(data), xo, sigma]
	popt, pcov = scipy.optimize.curve_fit(gauss_model, x, data, p0=guess)
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

def plotburst(data, retfig=False, extent=None):
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
	axs['B'].plot(np.nanmean(data, axis=1), np.linspace(*extent[2:], num=data.shape[0]))
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

		# 1 burst:
		pktime = np.nanargmax(np.nanmean(wfall, axis=0))*res_time*1000
		t_popt, _ = fitgauss(tseries, duration)
		window = 2*abs(t_popt[2]) # 2*stddev of a guassian fit to the integrated time series

		##### 4 component model: use 4 1d gaussians to make 4 windows
		xos = [9.26, 15.75, 22.61, 27.94] # ms
		xos_chans = np.floor(np.array(xos)/res_time_ms)
		t4_popt, _ = fit4gauss(tseries, duration, xos=xos)
		windows4 = 2*np.abs(t4_popt[8:]) # 2*sigma
		subfalls = []
		for xoi, s in zip(xos_chans, t4_popt[8:]):
			s = np.floor(4*np.abs(s)/res_time_ms)
			subfall = wfall[..., int(xoi-s):int(xoi+s)] # 4sigma window around burst
			subfall = driftrate.subtractbg(subfall, 0, int(subfall.shape[1]*0.1)) # again
			subfalls.append(subfall)

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

		#### Fitting (4 subburst model)
		dtdnus, subdfs = [], []
		subpeaks = []
		for subfall, xosi, win4i in zip(subfalls, xos, windows4):
			subdf = fitrows(subfall)
			subdf = subdf[(subdf.amp > 0)]
			subpktime = np.nanargmax(np.nanmean(subfall, axis=0))*res_time*1000
			subdf = subdf[(subpktime-win4i < subdf[tpoint]) & (subdf[tpoint] < subpktime+win4i)]

			dtdnu, dtdnu_err = measuredtdnu(subdf, subpktime)

			# Sub-burst plot
			subfig, subaxs = plotburst(
				subfall,
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
			AB
			AC
			AD
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

		### Individual Profiles --> Replace with spectrum summed from 2sigma timeseries profile
		profiles = arrtimesdf.sort_values('amp', ascending=False).head(3)
		chans = profiles.index
		lines = chans*res_freq + freqs_bin0
		for line in lines:
			ax_wfall.axhline(
				y=line,
				c='w',
				linewidth=0.75,
				ls='-.'
			)

		for irow, panel in zip(profiles.iterrows(), ['B', 'C', 'D']):
			chan, row = irow
			x = np.linspace(0, duration, num=wfall.shape[1])*1000 # times in ms
			xm = np.linspace(0, duration, num=2000)*1000 # times in ms
			axs[panel].plot(x-xobar, wfall[chan,:]/np.max(wfall)) # match fit normalization
			axs[panel].axvline(x=row['xo']-xobar, ls='--', c='k')
			axs[panel].set_title(f'{row["freqs"]:.1f} MHz Time series', fontsize='medium')

		axs['T'].sharex(axs['A'])
		axs['A'].sharex(axs['B'])
		axs['B'].sharex(axs['C'])
		axs['C'].sharex(axs['D'])
		axs['D'].set_xlabel('Time (ms)')
		axs['B'].set_xlim(row['xo']-xobar - 0.015, row['xo']-xobar + 0.015)
		ax_wfall.set_xlim(extent[0], extent[1])
		ax_wfall.set_ylim(extent[2], extent[3])
		plt.setp(ax_tseries.get_xticklabels(), visible=False)
		plt.setp(axs['B'].get_xticklabels(), visible=False)
		plt.setp(axs['C'].get_xticklabels(), visible=False)

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

