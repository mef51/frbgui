from frbrepeaters import driftlaw
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle
import warnings, importlib
import scipy
warnings.filterwarnings('ignore')
importlib.reload(driftlaw)

pd.options.mode.chained_assignment = None  # default='warn'
pd.set_option('display.max_rows', None)
%matplotlib inline

sources_filenames = [
	"results/frb121102_michilli_dmvariations.csv",
	'results/gajjar11a_dmvariations_oct15_2021.csv',
	'results/dmvars_180916.csv',
	'results/dmvars_180814.csv',
	'results/FRB180301_DM516-518-5_luo_272rows_Aug21.csv',
	'results/FRB121102_oostrum_525rows_Aug28.csv',
	'results/FRB121102_aggarwal_results_1078rows_Dec08.csv',
	# 'results/FRB121102_aggarwal_results_1122rows_Jan14.csv',
	'results/GajjarFRB121102_results_506rows_Oct19.csv',
	# 'results/simulated_nov15.csv'
]

sources = []
for sfile in sources_filenames:
	df = pd.read_csv(sfile)
	df = df.set_index("name")
	print(sfile)
	df = driftlaw.computeModelDetails(df)
	sources.append(df)

names = [
	'FRB20121102A Michilli et al. (2018) DM Variations',
	'FRB20121102A 11A Gajjar et al. (2018) DM Variations',
	'FRB20180916B bursts DM Variations',
	'FRB20180814A bursts DM Variations',
	'FRB20180301A bursts DM Variations',
	'FRB20121102A Oostrum et al. (2020) DM Variations',
	'FRB20121102A Aggarwal et al. (2021) DM Variations',
	'FRB20121102A Gajjar et al. (2018) DM Variations',
	# 'Simulated DM Variations'
]
annotations = [
	[0, 1, 2, 3, 4, 5],
	[],
	[3],
	[0, 1, 2, 3, 4, 5],
	[0, 1, 2, 3, 4, 5],
	[],
	# [],
	[],
	[]
]
dmlimited = [
	[3, 5, 6, 14],
	[],
	[15.5, 16, 23, 24, 29, 31, 32, 33],
	[180919.0],
	[],
	[],
	[],
	[],
	[]
]
exclusions = [
	[],
	[],
	[],
	[180917.0], # first is a multi burst, second has negative drift (which should better be corrected by limiting 180814's DM range)
	[],
	[],
	['B023', 'B016', 'B028', 'B056', 'B109'],
	[],
	[]
]

targetDMs   = [
	558.816263157895,
	568.333333333333,
	348.82,
	188.8,
	516.3571428571429,
	560.526315789474,
	560,
	556,
	565
]

def getFRB20121102AFit(pltdata):
	""" merge michilli and gajjar and fit together"""
	frb121102frame = pd.concat([pltdata['frames'][0],pltdata['frames'][1]])
	_, fits = driftlaw.plotSlopeVsDuration([frb121102frame], ['merged'], logscale=True, fiterrorfunc=driftlaw.log_error)
	return fits[0][1], fits[0][2]



########### pasted original below

fig1pltdata = { 'frames': [], 'labels': [], 'colors': cycle(['r', 'r', 'b', 'g', 'y', 'c', 'tomato', 'c']) }
fitdata = pd.DataFrame()
Bconstants = []
for source, name, annotate, exclude, targetDM in zip(sources, names, annotations, exclusions, targetDMs):
	frames = []
	labels = []

	print('\n#', name.split(' DM')[0]+':')
	DMrangeLogger('', source)
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

	DMrangeLogger(name.split(' DM')[0], source, burstwise=True) # log the dm range after we're done excluding
	source = driftlaw.sloperanges(source) # compute drift ranges after burst exclusions

	# special case: exclude 180911 at 190 pc/cm3 cuz drift range switches sign
	# TODO remove in favor of relative error (see above)
	# or exclude if 'slope_nu_minus' is negative but the other way seems more flexible
	if name == names[3]:
		source = source[~((source.DM == 190) & (source.slope_nu_min < 0))]
	source = driftlaw.sloperanges(source) # recompute range

	for dm, color in zip(source.DM.unique(), cycle(['r', 'y', 'b', 'k', 'g', 'c', 'm'])):
		df = source[source.DM == dm]
		df['color'] = color
		frames.append(df)
		labels.append(dm)

		# Figure 1
		if np.isclose(dm, targetDM):
			print(f'>> num bursts remaining = {len(df.index.unique())}')
			df['color'] = next(fig1pltdata['colors'])
			fig1pltdata['frames'].append(df)
			fig1pltdata['labels'].append(name.split(' DM')[0])

		# For Ref figure A, turn off otherwise
		# if dm == 565:
		#     df['color'] = 'b'
		# if dm == 348.82:
		#     df['color'] = 'r'

		# Figure 5: compute drift/nu^2 vs nu range of fits
		nus = np.linspace(0, 20000, num=2000)
		popt, pcov = scipy.optimize.curve_fit(lambda x, c: c, nus, df['slope_abs_nuobssq'])
		Bconstants.append((popt[0], np.sqrt(pcov[0])))

	# driftlaw.plotSlopeVsDuration(frames, labels, title=name, annotatei=annotate, logscale=True)
	markcolors = []
	for f in frames:
		markcolors.append(f['color'].iloc[0])
	markcolors = cycle(markcolors)
	# manylines = ['r--', 'y-', 'b-.', 'k--', 'g--', 'c-.', 'r-', 'b-', 'k-', 'g-', 'c-']
	linestyles = ['--', '-', '-.', '--', '--', '-.', '-', '-', '-', '-', '-']
	manylines = [next(markcolors)+lst for lst in linestyles]
	labels = [round(lbl, 2) for lbl in labels]
	tempax, fits = driftlaw.plotSlopeVsDuration(frames, labels, annotatei=[], fitlines=manylines, logscale=True, fiterrorfunc=driftlaw.log_error, dmtrace=True, hidefitlabel=True)
	if name.split(' DM')[0] == 'FRB20180916B bursts':
		hdls, lbls = tempax.get_legend_handles_labels()
		order = [2, 1, 0, 3, 4]
		tempax.legend([hdls[idx] for idx in order], [lbls[idx] for idx in order], fontsize='x-large')
	if name.split(' DM')[0] == 'FRB20121102A Gajjar et al. (2018)':
		tempax.set_xlim(0.073, 0.389)
		tempax.legend(fontsize='x-large', loc=1)

	for fit in fits:
		fitdata = fitdata.append({'name': name.split(' DM')[0],
								  'DM': fit[0],
								  'param': fit[1],
								  'err': fit[2] ,
								  'red_chisq': fit[3],
								  'numbursts': fit[5]},
								  ignore_index=True)

plt.close('all')

# fig1pltdata['frames'][0] = pd.concat([fig1pltdata['frames'][0], fig1pltdata['frames'][1]])
# del fig1pltdata['frames'][1]
# del fig1pltdata['labels'][1]
annotatei = [6]
ax, _ = driftlaw.plotSlopeVsDuration(fig1pltdata['frames'], fig1pltdata['labels'], fitlines=['r-', 'r--', 'b--', 'g-.'], markers=['o', 'd', 'p', 'X'], title=None,
									 annotatei=annotatei, fitextents=[0.0728, 15], hidefit=[0, 1],
									 logscale=True, errorfunc=driftlaw.rangeerror, fiterrorfunc=driftlaw.log_error)

## Mark all bursts whose DM ranges were limited by the exclusion process
print(f"{len(fig1pltdata['frames']) = }, {len(sources) = }")
for i, limitedpoints in enumerate(dmlimited):
	for lp in limitedpoints:
		# print('YOO: ', i, limitedpoints, lp)
		# print(fig1pltdata['frames'][i].head())
		v = fig1pltdata['frames'][i].loc[lp]
		plt.scatter(v['tau_w_ms'], v['slope_over_nuobs'], facecolors='none', edgecolors='k', linewidth=1.0, s=400)
		# ax.annotate('o', (v['tau_w_ms'], v['slope_over_nuobs']), xytext=(-10,-18), textcoords='offset pixels', weight='bold', size=25)

# see Josephy et al. 2019
C1data = ['C1', 1.9077, 184, 63, -3.92, 0.01, 630, 0.003, 3, 11.4, 11.4/1000, 0.12, 1]
josephyburst = pd.DataFrame([C1data], columns=['name', 'angle', 'sigmax', 'sigmay', 'slope (mhz/ms)', 'slope error (mhz/ms)', 'center_f', 'time_res', 'freq_res', 'tau_w_ms', 'tau_w', 'tau_w_error', 'red_chisq']).set_index('name')
josephyburst['slope_over_nuobs'] = -1*josephyburst['slope (mhz/ms)']/josephyburst['center_f']
xerr, yerr = driftlaw.modelerror(josephyburst)
ax.scatter(josephyburst['tau_w_ms'], josephyburst['slope_over_nuobs'], c='r', s=125, edgecolors='k', marker='s', label='FRB20121102A Josephy et al. (2019)')
ax.errorbar(josephyburst['tau_w_ms'], josephyburst['slope_over_nuobs'], xerr=xerr[0], yerr=yerr[0], ecolor='r')

# Plot shaded regions
alldata = pd.concat([f for f in fig1pltdata['frames']+[josephyburst]])
x = np.linspace(min(alldata['tau_w_ms'])*0.9, max(alldata['tau_w_ms'])*1.1, num=1200)
# plt.rcParams['hatch.linewidth'] = 1
for name, color, hatch in zip(fitdata.name.unique(), fig1pltdata['colors'], ['', '/', '\\', '-', '', '', '', '']):
	params = fitdata.loc[(fitdata.name == name) & (fitdata.numbursts > 1.0)]['param']
	print(f"{name} ({min(params):.3f} - {max(params):.3f}) $t_w^{{-1}}$")
	if 'Michilli' not in name: # michilli region is contained in gajjar region so no need to plot it
		lstr = f"{name.split(' ')[0]} range of fits"
		lstr = name
		ax.fill_between(x, min(params)/x, max(params)/x, alpha=0.1, color=color, edgecolor=color, label=lstr, zorder=-10)
		# ax.fill_between(x, min(params)/x, max(params)/x, alpha=0.1, color='None', edgecolor='w', hatch=hatch, zorder=-9)

## Combine michilli and gajjar fits. See `getFRB20121102AFit()`. the result of the fit when merging frames is 0.083 +- 0.005
# frb121102_A = getFRB20121102AFit(fig1pltdata) # opens a figure
frb121102_A = (0.07816130353202802, 0.005971946986254264)
plt.plot(x, frb121102_A[0]/x, 'r-', label='{} fit ({:.3f} $\\pm$ {:.3f}) $t_w^{{-1}}$'.format('FRB20121102A', frb121102_A[0], frb121102_A[1]))

# Legend stuff
hdls, lbls = ax.get_legend_handles_labels()
# order = [3, 4, 7, 5, 6, 2, 0, 1, 8, 9, 10]
# ax.legend([hdls[idx] for idx in order], [lbls[idx] for idx in order], fontsize='x-large')
# hdls[0], hdls[1], hdls[2], hdls[5], hdls[6], hdls[7] = hdls[2], hdls[0], hdls[1], hdls[7], hdls[5], hdls[6]
# lbls[0], lbls[1], lbls[2], lbls[5], lbls[6], lbls[7] = lbls[2], lbls[0], lbls[1], lbls[7], lbls[5], lbls[6]
# ax.legend(fontsize='medium')
ax.legend(hdls[8:18], lbls[8:18], fontsize='medium')
# ax.get_legend().remove()

# ax.set_xlim(0.0737, 12.5)
# ax.set_ylim(0.0013395279509532563, 508.64764364084175)
# plt.savefig('figures/driftvnuobs_dmvars.pdf')
