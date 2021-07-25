import dpg
import frbrepeaters
import driftrate, driftlaw
import os, glob, itertools
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import datetime
import warnings
warnings.filterwarnings("ignore")

# subfall, pkidx = frbrepeaters.loadpsrfits('data/oostrum2020/R1_frb121102/R1_B07.rf')
# subfall, pkidx, wfall = frbrepeaters.loadpsrfits('data/oostrum2020/R1_frb121102/R1_B07.rf')

# width = 150
# subfall = subfall[:, pkidx-width:pkidx+width]
# plt.imshow(subfall[:, pkidx-width:pkidx+width], origin='lower',interpolation='none', aspect='auto')
# plt.show()
# dpg.show_demo()

# dpg.show_debug()
dpg.show_logger()

# GUI data is stored in this object. Defaults initialized here
gdata = {
	'globfilter'     : '*.npz',
	'masks'          : {}, # will store lists of masked channel numbers
	'datadir'        : 'B:\\dev\\frbrepeaters\\data\\luo2020\\180813_ar_file\\ar_file\\converted',
	'results'        : [], # avoid using, might remove in favor of resultsdf
	'resultsdf'      : None,
	'popt'           : None, # currently displayed 2d fit
	'p0'             : None # currently displayed initial fit guess
}

def getscale(m, M=-1): # used for dynamically rounding numbers for display
	if M == -1: M = m
	ret = 1
	c = abs((m+M)/2)
	while c < 1:
		c *= 10; ret += 1
	return ret

def applyMasks(wfall):
	for mask in gdata['masks'][gdata['currfile']]:
		if mask < len(wfall):
			wfall[mask] = 0
	return wfall

def makeburstname(filename):
	return filename[-12:]


def log_cb(sender, data):
	dpg.log_debug(f"{sender}, {data}")
def error_log_cb(sender, data):
	dpg.log_error(f"{sender}, {data}")

def loaddata_cb(sender, data):
	wfall = None
	if 'fileidx' in data.keys():
		filename = gdata['files'][data['fileidx']]
		gdata['currfile'] = filename
		burstname = makeburstname(filename)
		dpg.set_value('burstname', burstname)
		loaded = np.load(filename)

		if type(loaded) == np.ndarray:
			wfall = loaded
		elif type(loaded) == np.lib.npyio.NpzFile:
			wfall = loaded['wfall']
			gdata['burstmeta'] = {}
			for key in loaded.files:
				if key != 'wfall':
					gdata['burstmeta'][key] = loaded[key]
					if key == 'dfs':
						dfs = loaded[key]
						downf = loaded['raw_shape'][0] / wfall.shape[0]
						df = (dfs[-1] - dfs[0])/len(dfs) * downf
						gdata['burstmeta']['fres'] = df
						dpg.set_value('df', df)
						dpg.configure_item('df', format='%.{}f'.format(getscale(df)+1))
					elif key == 'dt':
						downt = loaded['raw_shape'][1] / wfall.shape[1]
						dt = loaded[key][0] / loaded['raw_shape'][1] * downt * 1000
						gdata['burstmeta']['tres'] = dt
						dpg.set_value('dt', dt)
						dpg.configure_item(key, format='%.{}f'.format(getscale(dt)+1))
					else:
						dpg.set_value(key, loaded[key]) # this line sets all the burst fields

			# initialize DM range elements
			gdata['displayedDM'] = loaded['DM']
			gdata['burstDM']     = loaded['DM']
			dpg.set_value('dmdisplayed', str(gdata['displayedDM']))
			if dpg.get_value('dmrange')[0] == 0:
				dmrange = [gdata['burstDM']*0.99, gdata['burstDM']*1.01]
				dpg.set_value('dmrange', dmrange)
				dpg.configure_item('dmrange', speed=0.1)
				dmrange_cb(sender, None)

		gdata['wfall']          = wfall          # wfall at burstdm
		gdata['wfall_original'] = np.copy(wfall) # wfall at burstdm without additional subsampling

		# update subsample controls
		dpg.set_value('Wfallshapelbl', 'Maximum Size: {}'.format(np.shape(wfall)))
		dpg.set_value('Subfallshapelbl', 'Current Size: {}'.format(np.shape(wfall)))
		dpg.configure_item('dfreqinput', enabled=True, min_value=0, max_value=wfall.shape[0])
		dpg.configure_item('dtimeinput', enabled=True, min_value=0, max_value=wfall.shape[1])
		dpg.set_value('dfreqinput', wfall.shape[0])
		dpg.set_value('dtimeinput', wfall.shape[1])

		# update result controls
		if gdata['resultsdf'] is not None:
			hasResults = burstname in gdata['resultsdf'].index.unique()
			dpg.configure_item('PrevDM', enabled=hasResults)
			dpg.configure_item('NextDM', enabled=hasResults)
			dpg.configure_item('DeleteBurstBtn', enabled=False) # unimplemented
			dpg.configure_item('DeleteAllBtn', enabled=False) # unimplemented
			if hasResults:
				gdata['burstdf'] = gdata['resultsdf'].loc[burstname]
				updateResultTable(gdata['burstdf'])
				gdata['popt'] = gdata['burstdf'].set_index('DM').loc[gdata['displayedDM']][5:11]
			else:
				updateResultTable(pd.DataFrame())
				gdata['popt'] = None
				gdata['p0'] = None
				dpg.set_value('EnableP0Box', False)
				items = ['P0AllDMsBox','AmplitudeDrag','AngleDrag','x0y0Drag','SigmaXYDrag', 'RedoBtn']
				toggle_config('EnableP0Box', {'kwargs': ['enabled'], 'items': items})

	elif sender == 'subsample_cb' and data['subsample']: # ie. sender == 'subsample_cb' dpg.get_value('DM')
		wfall = gdata['wfall']
		dpg.set_value('Subfallshapelbl', 'Current Size: {}'.format(np.shape(wfall)))
	else:
		wfall = gdata['wfall']

	if gdata['currfile'] not in gdata['masks'].keys():
		gdata['masks'][gdata['currfile']] = [] # initialize list

	if wfall.shape == gdata['wfall_original'].shape:
		wfall = applyMasks(np.copy(gdata['wfall_original']))

	gdata['wfall'] = wfall
	# gdata['ts']    = np.nanmean(wfall, axis=0) # time series at burstDM
	# gdata['pkidx'] = np.nanargmax(gdata['ts']) # pkidx at burstDM, for displaying across DMs

	plotdata_cb(sender, data)

def cropwfall(wfall, twidth=150, pkidx=None):
	wfall = wfall.copy()
	ts    = np.nanmean(wfall, axis=0)
	if not pkidx:
		pkidx = np.nanargmax(ts)
	return wfall[..., pkidx-twidth:pkidx+twidth]

def getcorr2dtexture(corr, popt=None, p0=None):
	plt.figure(figsize=(5, 5))

	plt.imshow(corr, origin='lower', interpolation='none', aspect='auto', cmap='gray')
	plt.clim(0, np.max(corr)/20)

	if popt is not None and popt[0] > 0:
		fitmap = driftrate.makeFitmap(popt, corr)
		plt.contour(fitmap, [popt[0]/4, popt[0]*0.9], colors='b', alpha=0.75, origin='lower')

	if p0 is not None and p0[0] > 0:
		fitmap = driftrate.makeFitmap(p0, corr)
		plt.contour(fitmap, [p0[0]/4, p0[0]*0.9], colors='g', alpha=0.75, origin='lower')

	# remove axes, whitespace
	plt.gca().set_axis_off()
	plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
	plt.margins(0,0)
	plt.gca().xaxis.set_major_locator(plt.NullLocator())
	plt.gca().yaxis.set_major_locator(plt.NullLocator())

	fig = plt.gcf()
	fig.canvas.draw()
	texture = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
	w, h = fig.get_size_inches()*fig.dpi
	return texture, int(w), int(h)

def plotdata_cb(sender, data):
	if not data:
		data = {}

	wfall       = gdata['wfall']
	df, dt      = gdata['burstmeta']['fres'], gdata['burstmeta']['tres']
	lowest_freq = gdata['burstmeta']['dfs'][0] # mhz

	ddm = gdata['displayedDM'] - gdata['burstDM']
	wfall_dd = driftrate.dedisperse(wfall, ddm, lowest_freq, df, dt)
	wfall_dd_cr = cropwfall(wfall_dd)
	tseries = np.nanmean(wfall_dd_cr, axis=0)

	extents, correxts = driftrate.getExtents(wfall_dd_cr, df=df, dt=dt, lowest_freq=lowest_freq)
	gdata['extents'], gdata['correxts'] = extents, correxts

	corr = driftrate.autocorr2d(wfall_dd_cr)

	## enable scale sliders
	mostmin, mostmax = np.min(wfall_dd_cr), np.max(wfall_dd_cr)
	dpg.configure_item('wfallscale', enabled=True, min_value=mostmin, max_value=mostmax,
						format='%.{}f'.format(getscale(mostmin, mostmax)+1))
	mmincorr, mmaxcorr = np.min(corr), np.max(corr)
	dpg.configure_item('corrscale', enabled=True, min_value=mmincorr, max_value=mmaxcorr,
						format='%.{}f'.format(getscale(mmincorr, mmaxcorr)+1))

	if 'scale' in data.keys() and data['scale'] != None:
		smin, smax = data['scale']
	else:
		smin, smax = mostmin, mostmax
	dpg.set_value('wfallscale', [smin, smax])

	if 'cscale' in data.keys() and data['cscale'] != None:
		scmin, scmax = data['cscale']
	else:
		scmin, scmax = mmincorr, mmaxcorr
	dpg.set_value('corrscale', [scmin, scmax])

	wx, wy = dpg.get_plot_xlimits('WaterfallPlot'), dpg.get_plot_ylimits('WaterfallPlot')

	dpg.add_heat_series("WaterfallPlot", "Waterfall",
		values=list(np.flipud(wfall_dd_cr).flatten()),
		rows=wfall_dd_cr.shape[0], columns=wfall_dd_cr.shape[1],
		scale_min=smin, scale_max=smax,
		bounds_min=(0,0), bounds_max=(wfall_dd_cr.shape[1], wfall_dd_cr.shape[0]), format='')
		# bounds_min=(extents[0],extents[2]), bounds_max=(extents[1], extents[3]), format='')

	dpg.set_plot_xlimits_auto('WaterfallPlot')
	dpg.set_plot_ylimits_auto('WaterfallPlot')
	if 'keepview' in data.keys() and data['keepview']:
		## BROKEN
		# print('keepving view', wx, wy)
		# dpg.set_plot_xlimits('WaterfallPlot', wx[0], wx[1])
		# dpg.set_plot_ylimits('WaterfallPlot', wy[0], wy[1])
		pass

	p0 = gdata['p0'] if dpg.get_value('EnableP0Box') else None
	corr2dtexture, txwidth, txheight = getcorr2dtexture(corr, gdata['popt'], p0)
	dpg.add_texture("corr2dtexture", corr2dtexture, txwidth, txheight, format=dpg.mvTEX_RGB_INT)
	dpg.add_image_series("Corr2dPlot", "Corr2d", "corr2dtexture",
		bounds_min=[correxts[0],correxts[2]], bounds_max=[correxts[1], correxts[3]]
	)
	# dpg.add_heat_series("Corr2dPlot", "Corr2d",
	# 	values=list(np.flipud(corr).flatten()),
	# 	rows=corr.shape[0], columns=corr.shape[1],
	# 	scale_min=scmin, scale_max=scmax,
	# 	# bounds_min=(0,0), bounds_max=(corr.shape[1], corr.shape[0]), format='')
	# 	bounds_min=(correxts[0],correxts[2]), bounds_max=(correxts[1], correxts[3]), format='')

	dpg.add_line_series("TimeSeriesPlot", "TimeSeries", list(range(0, len(tseries))), tseries)

def subsample_cb(sender, data):
	df, dt = dpg.get_value("dfreqinput"), dpg.get_value("dtimeinput")
	print(df, dt)

	try:
		# Make a copy of the original fall, apply the masks, then downsample
		wfall = applyMasks(np.copy(gdata['wfall_original']))
		subfall = driftrate.subsample(wfall, df, dt)
		gdata['wfall'] = subfall
		log_cb('subsample_cb', (df, dt))
		loaddata_cb('subsample_cb', {'subsample': True})
	except (ValueError, ZeroDivisionError) as e:
		error_log_cb('subsample_cb', (df, dt, e))

def directory_cb(sender, data):
	dpg.set_value('Dirtext', 'Selected: {}'.format(data[0]))
	dpg.configure_item('Filter', enabled=True)
	dpg.configure_item('Clear filter', enabled=True)
	files = glob.glob(data[0]+'/{}'.format(gdata['globfilter']))
	dpg.configure_item('burstselect', items=[os.path.basename(x) for x in files])
	gdata['datadir'] = data[0]
	gdata['files']   = files
	log_cb(sender, data[0])

def filter_cb(sender, data):
	globfilter = dpg.get_value('Filter')
	if globfilter == '':
		globfilter = '*'
	gdata['globfilter'] = globfilter
	directory_cb(sender, [gdata['datadir']])

def clearfilter_cb(s, d):
	dpg.set_value('Filter', '')
	filter_cb(s, d)

def burstselect_cb(sender, data):
	fileidx = dpg.get_value('burstselect')
	loaddata_cb(sender, {'fileidx': fileidx})
	log_cb(sender, 'Opening file {}'.format(gdata['files'][fileidx]))

def exportmask_cb(sender, data):
	np.save('masks_{}.npy'.format('test'), [gdata['masks']])
	print(gdata['masks'])

def importmask_cb(sender, data):
	filename = data[1]
	log_cb(sender, 'mask selected: {}'.format(data))
	if filename.split('.')[-1] == 'npy':
		masks = np.load(filename, allow_pickle=True)[0]
		if type(masks) == dict:
			gdata['masks'] = masks
			loaddata_cb(sender, {'keepview': True})
			masktable_cb(sender, None)
		else:
			error_log_cb(sender, 'invalid mask dictionary selected.')
	else:
		error_log_cb(sender, 'invalid mask file selected.')

def removemask_cb(sender, data):
	coords = dpg.get_table_selections('Masktable') # should be length 1?
	coord = coords[0]
	mask = int(dpg.get_table_item('Masktable', coord[0], coord[1]))
	# print(type(mask), type(gdata['masks'][gdata['currfile']]), type(gdata['masks'][gdata['currfile']][0]))
	# print(mask, gdata['masks'][gdata['currfile']], mask in gdata['masks'][gdata['currfile']])
	if mask in gdata['masks'][gdata['currfile']]:
		gdata['masks'][gdata['currfile']].remove(mask)
		dpg.log_debug('removing {} from {} mask'.format(mask, gdata['currfile']))
		loaddata_cb(sender, {'keepview': True})
		masktable_cb(sender, None)

def masktable_cb(sender, data):
	# dpg makes working with tables impossible so we will delete the table and re-add it every time
	dpg.delete_item('Masktable')

	tableheight = round(min(25*len(gdata['masks'][list(gdata['masks'].keys())[0]]), 250))
	dpg.add_table('Masktable', [], height=tableheight, parent='Masking', callback=removemask_cb)

	columns = [s.split('.')[0][-8:] for s in gdata['masks'].keys()]
	for key, col in zip(gdata['masks'].keys(), columns):
		dpg.add_column('Masktable', col, gdata['masks'][key])

def resulttable_cb(sender, data):
	coords = dpg.get_table_selections('Resulttable')
	newsel = coords[0].copy()
	for coord in coords:
		dpg.set_table_selection('Resulttable', coord[0], coord[1], False)
	# dpg.set_table_selection('Resulttable', newsel[0], newsel[1], True)

	displaydm_cb('User', {'newdmidx': newsel[0]})

def updateResultTable(resultsdf):
	dpg.delete_item('Resulttable')
	dpg.add_table('Resulttable', [], height=225, parent='ResultsGroup', callback=resulttable_cb)

	# subset of driftrate.columns:
	columns = ['name', 'DM', 'amplitude', 'slope (mhz/ms)', 'theta', 'center_f']
	dpg.set_headers('Resulttable', columns)

	# [burstname, trialDM, center_f, slope, slope_err, theta, red_chisq], popt, perr, [fres_MHz, tres_ms/1000]
	for burstname, rowdata in resultsdf.iterrows():
		newrow = [burstname] + [rowdata[col] for col in columns[1:]]
		dpg.add_row('Resulttable', newrow)

def mousepos_cb(sender, data):
	isOnWaterfall = dpg.is_item_hovered('WaterfallPlot')
	if isOnWaterfall:
		tchan, fchan = dpg.get_plot_mouse_pos()
		mask = round(fchan)
		if mask not in gdata['masks'][gdata['currfile']]:
			gdata['masks'][gdata['currfile']].append(mask)

		loaddata_cb(sender, {'keepview': True})
		masktable_cb(sender, None)
		log_cb('mousepos_cb ', [[tchan, fchan], isOnWaterfall])
	else:
		return

def dmrange_cb(sender, data):
	dmrange   = dpg.get_value('dmrange')
	dmrange   = np.round(dmrange, 4) # dpg.get_value introduces rounding errors
	numtrials = dpg.get_value('numtrials')
	burstDM = gdata['burstDM']
	if dmrange[1] < dmrange[0]:
		dmrange.sort()
		dpg.set_value('dmrange', dmrange)
	if not (dmrange[0] < burstDM < dmrange[1]):
		dpg.configure_item('DMWarning', show=True)
	else:
		dpg.configure_item('DMWarning', show=False)
	# trialDMs = np.append(np.linspace(dmrange[0], dmrange[1], num=numtrials), burstDM)
	gdata['trialDMs'] = np.linspace(dmrange[0], dmrange[1], num=numtrials)

def getCurrentBurst():
	burstname = dpg.get_value('burstname').replace(',', '')
	df, dt = gdata['burstmeta']['fres'], gdata['burstmeta']['tres']
	lowest_freq = gdata['burstmeta']['dfs'][0] # mhz
	wfall = gdata['wfall'].copy()
	wfall_cr = cropwfall(wfall)
	burstDM = gdata['burstDM']
	return burstname, df, dt, lowest_freq, wfall_cr, burstDM

def slope_cb(sender, data):
	dpg.set_value('SlopeStatus', 'Status: Calculating...')

	burstname, df, dt, lowest_freq, wfall_cr, burstDM = getCurrentBurst()

	trialDMs = np.append(gdata['trialDMs'], burstDM)
	results, burstdf = driftrate.processDMRange(burstname, wfall_cr, burstDM, trialDMs, df, dt, lowest_freq)
	gdata['results'] += results
	fileidx = dpg.get_value('burstselect')

	# burstdf = burstdf.loc[burstdf.index.unique()[0]]
	burstdf = burstdf.loc[burstname]
	gdata['popt'] = burstdf.set_index('DM').loc[gdata['displayedDM']][5:11]

	gdata['burstdf'] = burstdf
	if gdata['resultsdf'] is None:
		gdata['resultsdf'] = burstdf
	else:
		gdata['resultsdf'] = gdata['resultsdf'].append(burstdf)

	dpg.set_value('SlopeStatus', 'Status: Done.')
	dpg.set_value('NumMeasurementsText', "# of Measurements: {}".format(len(burstdf)))
	dpg.set_value('TotalNumMeasurementsText', "Total # of Measurements: {}".format(len(gdata['resultsdf'])))
	dpg.configure_item('PrevDM', enabled=True)
	dpg.configure_item('NextDM', enabled=True)
	dpg.configure_item('ExportResultsBtn', enabled=True)
	dpg.configure_item('EnableP0Box', enabled=True)
	initializeP0Group()
	updateResultTable(burstdf)
	# print(driftlaw.computeModelDetails(burstdf))
	plotdata_cb(sender, data)

def redodm_cb(sender, data):
	p0 = gdata['p0']
	burstname, df, dt, lowest_freq, wfall_cr, burstDM = getCurrentBurst()
	result, burstdf = driftrate.processDMRange(
		burstname, wfall_cr, burstDM, [float(gdata['displayedDM'])],
		df, dt, lowest_freq, p0=p0
	)

	df = gdata['resultsdf']
	df[(df.index == burstname) & (df['DM'] == gdata['displayedDM'])] = burstdf

	gdata['resultsdf'] = df
	gdata['burstdf'] = gdata['resultsdf'].loc[burstname]
	gdata['popt'] = burstdf.set_index('DM').loc[gdata['displayedDM']][5:11]

	updateResultTable(gdata['burstdf'])
	plotdata_cb(sender, data)

def exportresults_cb(sender, data):
	results = gdata['results']
	df = driftrate.exportresults(results)
	datestr = datetime.now().strftime('%b%d')
	prefix = dpg.get_value('ExportPrefix')
	filename = '{}_results_{}rows_{}.csv'.format(prefix, len(df.index), datestr)
	df.to_csv(filename)
	dpg.set_value('ExportResultText', 'Saved to {}'.format(filename))
	dpg.configure_item('ExportResultText', show=True)

def displaydm_cb(sender, data):
	burstDM = gdata['burstDM']
	trialDMs = np.append(gdata['trialDMs'], burstDM)
	dmlist = list(trialDMs)
	dmidx = dmlist.index(gdata['displayedDM'])

	if sender == 'User':
		newdmidx = data['newdmidx']
	elif sender == 'NextDM':
		newdmidx = dmidx + 1
		if not (newdmidx < len(dmlist)):
			newdmidx = 0
	elif sender == 'PrevDM':
		newdmidx = dmidx - 1
	gdata['displayedDM'] = dmlist[newdmidx]
	dpg.set_value('dmdisplayed', str(round(gdata['displayedDM'], getscale(gdata['displayedDM']))))

	df = gdata['burstdf'].set_index('DM')
	df = df.set_index(df.index.astype('float'))
	gdata['popt'] = df.loc[gdata['displayedDM']][5:11]
	plotdata_cb(sender, data)

def confirmpopup(data, cb):
	with popup("main", "ConfirmDelete", modal=True):
		add_text("All those beautiful results will be deleted.\nThis operation cannot be undone!")
		add_button("OK", width=75,
			callback=lambda s, d: cb('ConfirmDelete', data)
		)
		add_same_line()
		add_button("Cancel", width=75,
			callback=lambda s, d: cb("ConfirmDelete", data)
		)

def deleteresults_cb(sender, data):
	""" TODO: implement """
	burstname = dpg.get_value('burstname').replace(',', '')
	if sender == 'ConfirmDelete':
		if data == 'burst':
			df = gdata['resultsdf'].loc[~burstname] # nope
		elif data == 'all':
			gdata['resultsdf'] = None
	else:
		confirmpopup(data, deleteresults_cb)

	print(sender, data)

def initializeP0Group():
	params = ['amplitude', 'xo', 'yo', 'sigmax', 'sigmay', 'angle']
	p0 = []
	for burst, row in gdata['burstdf'].iterrows():
		if row['amplitude'] > 0:
			p0 = [row[param] for param in params]
	gdata['p0'] = p0
	dpg.set_value("AmplitudeDrag", p0[0])
	dpg.set_value("AngleDrag", p0[5])
	dpg.set_value("x0y0Drag", [p0[1], p0[2]])
	dpg.set_value("SigmaXYDrag", [p0[3], p0[4]])
	for item in ['AmplitudeDrag', 'AngleDrag', 'x0y0Drag', 'SigmaXYDrag']:
		val = dpg.get_value(item)
		if type(val) == list:
			val = val[0]
		dpg.configure_item(item, speed=1/10**getscale(val), format='%.{}f'.format(getscale(val)+1))

def enablep0_cb(sender, data):
	toggle_config(sender, data)
	updatep0_cb(sender, data)

def updatep0_cb(sender, data):
	p0 = []
	for item in ['AmplitudeDrag', 'x0y0Drag', 'SigmaXYDrag', 'AngleDrag']:
		val = dpg.get_value(item)
		if type(val) == list:
			p0 = p0 + val
		else:
			p0.append(val)
	gdata['p0'] = p0
	plotdata_cb(sender, data)

def helpmarker(message):
	dpg.add_same_line()
	dpg.add_text("(?)", color=[150, 150, 150], tip=message)

def toggle_config(sender, data):
	config_dict = {}
	for kwarg in data['kwargs']:
		config_dict[kwarg] = dpg.get_value(sender)
	for item in data['items']:
		dpg.configure_item(item, **config_dict)

### Analysis window
with dpg.window('FRB Analysis', width=560, height=745, x_pos=10, y_pos=30):
	with dpg.collapsing_header("1. Data", default_open=True):
		dpg.add_button("Select Directory...", callback=lambda s, d: dpg.select_directory_dialog(directory_cb))
		dpg.add_text("Dirtext", default_value="Selected: (no directory selected)")
		dpg.add_text("Filter:"); dpg.add_same_line()
		dpg.add_input_text("Filter", label='', hint="eg. *.npy", callback=filter_cb, enabled=False)
		dpg.add_same_line();
		dpg.add_button('Clear filter', callback=clearfilter_cb, enabled=False)

		dpg.add_text("Files found:")
		dpg.add_listbox("burstselect", label='', items=[], num_items=10, width=520,
						callback=burstselect_cb, tip="Select to load burst...")

		dpg.add_text("Burst Metadata:")
		dpg.add_input_text('burstname', label='Burst Name')
		dpg.add_input_float('DM', label='DM (pc/cm^3)')
		dpg.add_input_float('dt', label='Time Resolution (ms)')
		dpg.add_input_float('df', label='Freq Resolution (MHz)')
		dpg.add_input_float('center_f', label='Center Frequency (MHz)')
		dpg.add_input_float('bandwidth', label='Bandwidth (MHz)')
		dpg.add_input_float('duration', label='Data Duration')
		dpg.add_input_float('burstSN', label='Burst SNR')
		dpg.add_input_text('telescope', label='Telescope')
		## TODO: Use these units to populate the resolution inputs
		dpg.add_input_text('freq_unit', label='Frequency Unit')
		dpg.add_input_text('time_unit', label='Time Unit')
		dpg.add_input_text('int_unit', label='Intensity Unit')

		dpg.add_text('Dedispersion Range for all Bursts: ')
		dpg.add_text('DMWarning', default_value='Warning: Range chosen does not include burst DM',
			color=[255, 0, 0], show=False)
		dpg.add_drag_float2('dmrange', label='DM range (pc/cm^3)', callback=dmrange_cb,
			min_value=0, max_value=0)
		helpmarker('Double click to edit')
		dpg.add_input_int('numtrials', label='# of Trial DMs', default_value=10, callback=dmrange_cb)

	with dpg.collapsing_header("2. Burst Cleanup", default_open=False):
		with dpg.tree_node('Masking', default_open=True):
			dpg.add_text("Click on the waterfall plot to begin masking frequency channels.")
			dpg.add_text("NOTE: only mask on the original waterfall (todo: add a 'mask' button)")

			dpg.add_button('Export Masks', callback=exportmask_cb, enabled=True)
			dpg.add_same_line()
			dpg.add_button('Import Masks',
				callback=lambda s, d: dpg.open_file_dialog(importmask_cb, extensions='.npy'),
				enabled=True)

			dpg.add_table('Masktable', [], height=50)

		with dpg.tree_node('Downsampling', default_open=True):
			dpg.add_text("Wfallshapelbl", default_value="Maximum Size: (no burst selected)")
			dpg.add_text("Subfallshapelbl", default_value="Current Size: (no burst selected)")
			dpg.add_input_int("dfreqinput", width=100, label="df", callback=subsample_cb, enabled=False)
			dpg.add_same_line()
			dpg.add_input_int("dtimeinput", width=100, label="dt", callback=subsample_cb, enabled=False)

	with dpg.collapsing_header('SlopeSection', label="3. Sub-burst Slope Measurements", default_open=True):
		dpg.add_button("Measure Slope over DM Range", callback=slope_cb)
		dpg.add_same_line()
		dpg.add_text("SlopeStatus", default_value="Status: (click 'Measure Slope' to calculate)")
		dpg.add_button('DeleteBurstBtn',
			label="Delete this bursts's results",
			callback=deleteresults_cb,
			callback_data='burst',
			enabled=False
		)
		dpg.add_same_line()
		dpg.add_button('DeleteAllBtn',
			label="Delete ALL results",
			callback=deleteresults_cb,
			callback_data='all',
			enabled=False
		)

		dpg.add_text('NumMeasurementsText', default_value="# of Measurements for this burst: (none)")
		dpg.add_text('TotalNumMeasurementsText', default_value="Total # of Measurements: (none)")
		with dpg.group("DMselector", horizontal=True):
			dpg.add_text("DM Displayed (pc/cm^3): ")
			dpg.add_text("dmdisplayed", default_value=str(dpg.get_value('DM')))
			dpg.add_button("PrevDM", arrow=True, direction=dpg.mvDir_Left, enabled=False,
				callback=displaydm_cb)
			dpg.add_button("NextDM", arrow=True, direction=dpg.mvDir_Right, enabled=False,
				callback=displaydm_cb)
			helpmarker('Selecting a row also displays that DM')

		with dpg.group('P0Group'):
			with dpg.tree_node('Fit Initial Guess', default_open=False):
				items = ['P0AllDMsBox','AmplitudeDrag','AngleDrag','x0y0Drag','SigmaXYDrag', 'RedoBtn']
				dpg.add_checkbox("EnableP0Box", label='Use Initial Guess', default_value=False,
					enabled=False,
					callback=enablep0_cb,
					callback_data={'kwargs': ['enabled'], 'items': items}
				)
				dpg.add_same_line()
				dpg.add_checkbox('P0AllDMsBox', label='Use for all DMs', default_value=False, enabled=False)
				dpg.add_drag_float("AmplitudeDrag", label="Amplitude", enabled=False,
					min_value=0, max_value=0, callback=updatep0_cb)
				dpg.add_drag_float("AngleDrag", label="Angle", enabled=False,
					min_value=0, max_value=0, callback=updatep0_cb)
				dpg.add_drag_float2("x0y0Drag", label="x0, y0", enabled=False,
					min_value=0, max_value=0, callback=updatep0_cb)
				dpg.add_drag_float2("SigmaXYDrag", label="sigmax, sigmay", enabled=False,
					min_value=0, max_value=0, callback=updatep0_cb)
				dpg.add_button("RedoBtn", label="Redo this DM", callback=redodm_cb, enabled=False)

		with dpg.group('ResultsGroup'):
			dpg.add_table('Resulttable', [], height=10, parent='ResultsGroup', callback=resulttable_cb)

		dpg.add_input_text('ExportPrefix', label='Filename Prefix', default_value="FRBName")
		dpg.add_button('ExportResultsBtn', label="Export Full Results", callback=exportresults_cb, enabled=False)
		dpg.add_same_line()
		dpg.add_text('ExportResultText', default_value='', show=False, color=(0, 255, 0))

### Plotting window
with dpg.window("FRB Plots", width=1035, height=745, x_pos=600, y_pos=30):
	dpg.set_mouse_click_callback(mousepos_cb)
	dpg.add_slider_float2("wfallscale", label='Wfall Min/Max', enabled=False,
						  width=400, callback=plotdata_cb,
						  callback_data=lambda: {'scale': dpg.get_value('wfallscale')})
	dpg.add_same_line()
	dpg.add_slider_float2("corrscale", label='Corr Min/Max', enabled=False,
						  width=400, callback=plotdata_cb,
						  callback_data=lambda: {'cscale': dpg.get_value('corrscale')})

	# Colors: From implot.cpp: {"Default","Deep","Dark","Pastel","Paired","Viridis","Plasma","Hot","Cool","Pink","Jet"};
	dpg.add_plot("WaterfallPlot", no_legend=True, height=480, width=500,
				x_axis_name='Time (ms)', y_axis_name='Frequency (MHz)')
	dpg.set_color_map("WaterfallPlot", 5) # Viridis
	dpg.add_same_line()
	dpg.add_plot("Corr2dPlot", no_legend=True, height=480, width=500, x_axis_name='Time lag (ms)',
				 y_axis_name='Freq lag (MHz)')
	dpg.set_color_map("Corr2dPlot", 9) # Pink. "Hot" is good too

	dpg.add_plot("TimeSeriesPlot", x_axis_name="Time", y_axis_name="Intensity",
				height=200, width=500, no_legend=True)


### Main Menu Bar
with dpg.window("main"):
	with dpg.menu_bar("MenuBar##frbrs"):
		with dpg.menu("Menu##frbrs"):
			pass
		with dpg.menu("Themes##frbrs"):
			dpg.add_menu_item("Dark", callback = lambda sender, data: dpg.set_theme(sender), check=True)
			dpg.add_menu_item("Light", callback = lambda sender, data: dpg.set_theme(sender), check=True)
			dpg.add_menu_item("Classic", callback = lambda sender, data: dpg.set_theme(sender), check=True, shortcut="Ctrl+Shift+T")
			dpg.add_menu_item("Dark 2", callback = lambda sender, data: dpg.set_theme(sender), check=True)
			dpg.add_menu_item("Grey", callback = lambda sender, data: dpg.set_theme(sender), check=True)
			dpg.add_menu_item("Dark Grey", callback = lambda sender, data: dpg.set_theme(sender), check=True)
			dpg.add_menu_item("Cherry", callback = lambda sender, data: dpg.set_theme(sender), check=True)
			dpg.add_menu_item("Purple", callback = lambda sender, data: dpg.set_theme(sender), check=True)
			dpg.add_menu_item("Gold", callback = lambda sender, data: dpg.set_theme(sender), check=True, shortcut="Ctrl+Shift+P")
			dpg.add_menu_item("Red", callback = lambda sender, data: dpg.set_theme(sender), check=True, shortcut="Ctrl+Shift+Y")
		with dpg.menu("Tools##frbrs"):
			dpg.add_menu_item("Show Logger##frbrs", callback=dpg.show_logger)
			dpg.add_menu_item("Show About##frbrs", callback=dpg.show_about)
			dpg.add_menu_item("Show Metrics##frbrs", callback=dpg.show_metrics)
			dpg.add_menu_item("Show Documentation##frbrs", callback=dpg.show_documentation)
			dpg.add_menu_item("Show Debug##frbrs", callback=dpg.show_debug)
			dpg.add_menu_item("Show Style Editor##frbrs", callback=dpg.show_style_editor)
			# dpg.add_menu_item("Show Demo##frbrs", callback=dpg.show_demo)


# dpg.show_documentation()
dpg.set_main_window_size(1700, 850)
dpg.set_main_window_title("FRB Repeaters")
# dpg.start_dearpygui()

# Load defaults
dpg.set_value('Filter', gdata['globfilter'])
directory_cb('user', [gdata['datadir']]) # default directory
burstselect_cb('burstselect', None)
importmask_cb('user', ['B:\\dev\\frbrepeaters', 'luomasks.npy'])

## dm range defaults
dpg.set_value('dmrange', [515.2, 525.5])
dpg.set_value('numtrials', 2)
dmrange_cb('user', None)

dpg.start_dearpygui(primary_window='main')
