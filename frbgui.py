import dpg
import frbrepeaters
import driftrate, driftlaw
import os, glob, itertools
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import datetime
import warnings
from your.utils.rfi import sk_sg_filter
warnings.filterwarnings("ignore")

dpg.show_logger()

twidth_default = 150

# GUI data is stored in this object. Defaults initialized here and at the bottom
gdata = {
	'globfilter'     : '*.npz',
	'masks'          : {},                   # will store lists of masked channel numbers
	'datadir'        : '',
	'multiburst'     : {                     # store all multiburst metadata
						'numregions': 1,
						'enabled'   : False,
						'regions'   : {}
						},
	'resultsdf'      : None,
	'p0'             : None,                  # currently displayed initial fit guess
	'displayedBurst' : '',                    # name of displayed burst
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
	if 'sksgmask' in gdata and dpg.get_value('EnableSKSGMaskBox'):
		wfall[gdata['sksgmask'], :] = 0

	return wfall

def makeburstname(filename):
	return filename.split(os.sep)[-1].split('.')[0]

def log_cb(sender, data):
	dpg.log_debug(f"{sender}, {data}")
def error_log_cb(sender, data):
	dpg.log_error(f"{sender}, {data}")

def updatedata_cb(sender, data):
	wfall = None
	if 'fileidx' in data.keys(): # load burst from disk
		filename = gdata['files'][data['fileidx']]
		gdata['currfile'] = filename
		burstname = makeburstname(filename)
		dpg.set_value('burstname', burstname)
		loaded = np.load(filename)

		if type(loaded) == np.ndarray:
			wfall = loaded.astype(np.float64)
		elif type(loaded) == np.lib.npyio.NpzFile:
			wfall = loaded['wfall'].astype(np.float64)
			storedshape = wfall.shape
			gdata['burstmeta'] = {}
			for key in loaded.files:
				if key != 'wfall':
					gdata['burstmeta'][key] = loaded[key]
					if key == 'dfs':
						# dfs = loaded[key]
						gdata['burstmeta']['bandwidth'] = abs(loaded['bandwidth'])
						df = gdata['burstmeta']['bandwidth'] / wfall.shape[0]
						gdata['burstmeta']['fres_original'] = df
						gdata['burstmeta']['fres'] = df
						dpg.set_value('df', df)
						dpg.configure_item('df', format='%.{}f'.format(getscale(df)+1))
					elif key == 'dt':
						dt = loaded['duration'] / storedshape[1]*1000
						gdata['burstmeta']['duration'] = loaded['duration']
						gdata['burstmeta']['tres_original'] = dt
						gdata['burstmeta']['tres'] = dt
						dpg.set_value('dt', dt)
						dpg.configure_item(key, format='%.{}f'.format(getscale(dt)+3))
					else:
						dpg.set_value(key, loaded[key]) # this line sets all the burst fields

				# post loading tweaks
				if key == 'time_unit':
					dpg.configure_item('duration', label=f"Data Duration ({loaded[key]})")

			# initialize DM range elements
			gdata['displayedDM'] = loaded['DM']
			gdata['burstDM']     = loaded['DM']
			dpg.set_value('dmdisplayed', str(gdata['displayedDM']))
			dpg.set_value('burstdisplayed', burstname)
			if dpg.get_value('dmrange')[0] == 0:
				dmrange = [gdata['burstDM']*0.99, gdata['burstDM']*1.01]
				dpg.set_value('dmrange', dmrange)
				dpg.configure_item('dmrange', speed=0.1)
				dmrange_cb(sender, None)

		gdata['wfall']          = wfall          # wfall at burstdm
		gdata['wfall_original'] = np.copy(wfall) # wfall at burstdm without additional subsampling

		# update subsample controls
		dpg.set_value('Wfallshapelbl', 'Original Size: {}'.format(np.shape(wfall)))
		dpg.set_value('Subfallshapelbl', 'Current Size: {}'.format(np.shape(wfall)))
		dpg.configure_item('numfreqinput', enabled=True, min_value=0, max_value=wfall.shape[0])
		dpg.configure_item('numtimeinput', enabled=True, min_value=0, max_value=wfall.shape[1])
		dpg.configure_item('ResetSamplingBtn', enabled=True)
		dpg.set_value('numfreqinput', wfall.shape[0])
		dpg.set_value('numtimeinput', wfall.shape[1])

		# update meta controls
		dpg.set_value('twidth', twidth_default)
		dpg.configure_item('twidth', max_value=round(wfall.shape[1]/2))

		# update result controls and reproduce measurement state
		if gdata['resultsdf'] is not None:
			hasResults = burstname in gdata['resultsdf'].index.unique()
			dpg.configure_item('PrevDM', enabled=hasResults)
			dpg.configure_item('NextDM', enabled=hasResults)
			dpg.configure_item('ExportCSVBtn', enabled=hasResults)
			dpg.configure_item('ExportPDFBtn', enabled=hasResults)
			dpg.configure_item('DeleteBurstBtn', enabled=False) # unimplemented
			dpg.configure_item('DeleteAllBtn', enabled=False) # unimplemented
			dpg.configure_item('EnableP0Box', enabled=hasResults)
			if hasResults:
				if gdata['multiburst']['enabled']:
					gdata['burstdf'] = gdata['resultsdf'][gdata['resultsdf'].index.str.startswith(burstname)]
				else:
					gdata['burstdf'] = gdata['resultsdf'].loc[burstname]

				# reload waterfall cleanup
				cols = ['tsamp_width', 'subbg_start (ms)', 'subbg_end (ms)', 'sksigma', 'skwindow']
				twidth, subbgstart, subbgend, sksigma, skwindow = gdata['burstdf'][cols].iloc[0]
				if sender == 'burstselect':
					dpg.set_value('twidth', int(twidth))
					if not pd.isnull(subbgstart):
						dpg.set_value('EnableSubBGBox', True)
						toggle_config('EnableSubBGBox', {'kwargs': ['enabled'], 'items': ['SubtractBGRegion']})
						dpg.set_value('SubtractBGRegion', [subbgstart, subbgend])
					else:
						dpg.set_value('EnableSubBGBox', False)
						toggle_config('EnableSubBGBox', {'kwargs': ['enabled'], 'items': ['SubtractBGRegion']})

					skitems = ['SKSGSigmaInput','SKSGWindowInput','SKSGStatus']
					if not pd.isnull(sksigma):
						dpg.set_value("EnableSKSGMaskBox", True)
						dpg.set_value('SKSGSigmaInput', int(sksigma))
						dpg.set_value('SKSGWindowInput', int(skwindow))
						toggle_config('EnableSKSGMaskBox', {'kwargs': ['enabled'], 'items': skitems})
					else:
						dpg.set_value('EnableSKSGMaskBox', False)
						toggle_config('EnableSKSGMaskBox', {'kwargs': ['enabled'], 'items': skitems})

				# check if subsampling is needed (for eg. if results have been loaded)
				fres, tres = gdata['burstdf'][['f_res (mhz)', 'time_res (s)']].iloc[0]
				tres = tres*1000
				if fres != gdata['burstmeta']['fres']:
					subfac = fres/gdata['burstmeta']['fres']
					if round(subfac) == subfac:
						dpg.set_value('numfreqinput', int(wfall.shape[0]/subfac))
						wfall = subsample_cb('loadresults_cb', None)
				if tres != gdata['burstmeta']['tres']:
					subfac = tres/gdata['burstmeta']['tres']
					if round(subfac) == subfac:
						dpg.set_value('numtimeinput', int(wfall.shape[1]/subfac))
						wfall = subsample_cb('loadresults_cb', None)

				updateResultTable(gdata['burstdf'])
				initializeP0Group()
			else:
				updateResultTable(pd.DataFrame())
				gdata['burstdf'] = pd.DataFrame()
				gdata['p0'] = None
				dpg.set_value('EnableP0Box', False)
				items = ['P0AllDMsBox','AmplitudeDrag','AngleDrag','x0y0Drag','SigmaXYDrag', 'RedoBtn']
				toggle_config('EnableP0Box', {'kwargs': ['enabled'], 'items': items})

			dpg.set_value('NumMeasurementsText', "# of Measurements for this burst: {}".format(len(gdata['burstdf'])))
			dpg.set_value('TotalNumMeasurementsText', "Total # of Measurements: {}".format(len(gdata['resultsdf'])))

		# setup burst splitting
		for regid in range(1, gdata['multiburst']['numregions']):
			if dpg.does_item_exist('RegionSelector{}'.format(regid)):
				maxval = driftrate.cropwfall(wfall, twidth=dpg.get_value('twidth')).shape[1]*gdata['burstmeta']['tres']
				dpg.configure_item('Region{}'.format(regid), max_value=maxval, speed=maxval*0.005)
		if burstname in gdata['multiburst']['regions']:
			if not gdata['multiburst']['enabled']:
				dpg.set_value('MultiBurstBox', True)
				enablesplitting_cb('MultiBurstBox', {'kwargs': ['enabled']})
			regions = gdata['multiburst']['regions'][burstname]
			numregions = len(regions.keys())
			while numregions > gdata['multiburst']['numregions']-1:
				addregion_cb(sender, data)
			while numregions < gdata['multiburst']['numregions']-1:
				removeregion_cb([gdata['multiburst']['numregions']-1], data)
			# set the elements based on gdata['regions']
			for regid, name in enumerate(regions):
				regid += 1 # off by one
				regiontype = 0 if 'background' in name else 1 # 0 is background, 1 is burst
				dpg.set_value('Region{}'.format(regid), regions[name])
				dpg.set_value('RegionType{}'.format(regid), regiontype)
				drawregion_cb([regid], None)

	elif sender == 'subsample_cb' and data['subsample']: # ie. sender == 'subsample_cb' dpg.get_value('DM')
		wfall = gdata['wfall']

		disp_bandwidth = gdata['extents'][3] - gdata['extents'][2] # == gdata['burstmeta']['bandwidth']
		disp_duration = gdata['extents'][1] - gdata['extents'][0]
		twidth = disp_duration/gdata['burstmeta']['tres_original']/2
		dpg.set_value('twidth', round(twidth * (wfall.shape[1]/gdata['wfall_original'].shape[1])))

		wfall_cr = driftrate.cropwfall(wfall, twidth=dpg.get_value('twidth'))
		gdata['burstmeta']['fres'] = gdata['burstmeta']['bandwidth'] / wfall.shape[0]
		gdata['burstmeta']['tres'] = gdata['burstmeta']['duration']*1000 / wfall.shape[1]

		dpg.set_value('Subfallshapelbl', 'Current Size: {}'.format(np.shape(wfall)))
	else:
		wfall = gdata['wfall']

	if gdata['currfile'] not in gdata['masks'].keys():
		gdata['masks'][gdata['currfile']] = [] # initialize list

	if wfall.shape == gdata['wfall_original'].shape:
		wfall = applyMasks(np.copy(gdata['wfall_original']))
		if dpg.get_value('EnableSubBGBox'):
			tleft, tright = dpg.get_value('SubtractBGRegion')
			timerange = [round(gdata['extents'][0]), round(gdata['extents'][1])]
			tleft  = round(np.interp(tleft, timerange, [0, wfall.shape[1]]))
			tright = round(np.interp(tright, timerange, [0, wfall.shape[1]]))
			wfall = driftrate.subtractbg(wfall, tleft, tright)


	gdata['wfall'] = wfall
	# gdata['ts']    = np.nanmean(wfall, axis=0) # time series at burstDM
	# gdata['pkidx'] = np.nanargmax(gdata['ts']) # pkidx at burstDM, for displaying across DMs

	plotdata_cb(sender, data)

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

def plotdata_cb(_, data):
	if not data:
		data = {}

	df, dt      = gdata['burstmeta']['fres'], gdata['burstmeta']['tres']
	lowest_freq = min(gdata['burstmeta']['dfs']) # mhz
	ddm = gdata['displayedDM'] - gdata['burstDM']
	burstname = dpg.get_value('burstname').replace(',', '')

	popt = None
	if gdata['resultsdf'] is not None:
		subburstdf = gdata['resultsdf'][gdata['resultsdf'].index.str.startswith(burstname)]
		if not subburstdf.empty:
			dmframe = subburstdf.loc[burstname].set_index('DM')
			popt = dmframe.loc[np.isclose(dmframe.index, gdata['displayedDM'])].iloc[0][5:11]

	subname = None if 'resultidx' not in data else subburstdf.index[data['resultidx']]
	if ('resultidx' not in data) or (subname == burstname):
		gdata['displayedBurst'] = burstname
		wfall = gdata['wfall'].copy()
		wfall_dd = driftrate.dedisperse(wfall, ddm, lowest_freq, df, dt)
		wfall_dd_cr = driftrate.cropwfall(wfall_dd, twidth=dpg.get_value('twidth'))
	elif ('resultidx' in data) and (subname != burstname):
		gdata['displayedBurst'] = subname
		subbursts = getSubbursts()
		subburst = subbursts[subname]
		wfall_cr = subburst
		wfall_dd_cr = driftrate.dedisperse(wfall_cr, ddm, lowest_freq, df, dt)
		dmframe = subburstdf.loc[subname].set_index('DM')
		popt = dmframe.loc[np.isclose(dmframe.index, gdata['displayedDM'])].iloc[0][5:11]

	tseries = np.nanmean(wfall_dd_cr, axis=0)

	extents, correxts = driftrate.getExtents(wfall_dd_cr, df=df, dt=dt, lowest_freq=lowest_freq)
	gdata['extents'], gdata['correxts'] = extents, correxts
	dpg.set_value('twidth_ms', 2*dpg.get_value('twidth')*dt)
	# print(extents, df, dt, lowest_freq)

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
		# bounds_min=(0,0), bounds_max=(wfall_dd_cr.shape[1], wfall_dd_cr.shape[0]), format='')
		bounds_min=(extents[0],extents[2]), bounds_max=(extents[1], extents[3]), format='')

	dpg.set_plot_xlimits_auto('WaterfallPlot')
	dpg.set_plot_ylimits_auto('WaterfallPlot')
	if 'keepview' in data.keys() and data['keepview']:
		## BROKEN
		# print('keepving view', wx, wy)
		# dpg.set_plot_xlimits('WaterfallPlot', wx[0], wx[1])
		# dpg.set_plot_ylimits('WaterfallPlot', wy[0], wy[1])
		pass

	p0 = gdata['p0'] if dpg.get_value('EnableP0Box') else None
	corr2dtexture, txwidth, txheight = getcorr2dtexture(corr, popt, p0)
	dpg.add_texture("corr2dtexture", corr2dtexture, txwidth, txheight, format=dpg.mvTEX_RGB_INT)
	dpg.add_image_series("Corr2dPlot", "Corr2d", "corr2dtexture",
		bounds_min=[correxts[0],correxts[2]], bounds_max=[correxts[1], correxts[3]]
	)

	tx = np.linspace(extents[0], extents[1], num=len(tseries))
	dpg.add_line_series("TimeSeriesPlot", "TimeSeries", tx, tseries)

def twidth_cb(_, data):
	wfall_cr = getCurrentBurst()[4]
	for regid in range(1, gdata['multiburst']['numregions']):
		if dpg.does_item_exist('RegionSelector{}'.format(regid)):
			maxval = wfall_cr.shape[1]*gdata['burstmeta']['tres']
			dpg.configure_item('Region{}'.format(regid), max_value=maxval, speed=maxval*0.005)
	plotdata_cb(_, data)

def subsample_cb(sender, data):
	if sender == 'ResetSamplingBtn':
		dpg.set_value('numfreqinput', gdata['wfall_original'].shape[0])
		dpg.set_value('numtimeinput', gdata['wfall_original'].shape[1])

	numf, numt = dpg.get_value("numfreqinput"), dpg.get_value("numtimeinput")

	try:
		# Make a copy of the original fall, apply the masks, then downsample
		wfall = applyMasks(np.copy(gdata['wfall_original']))
		if dpg.get_value('EnableSubBGBox'):
			tleft, tright = dpg.get_value('SubtractBGRegion')
			timerange = [round(gdata['extents'][0]), round(gdata['extents'][1])]
			tleft  = round(np.interp(tleft, timerange, [0, wfall.shape[1]]))
			tright = round(np.interp(tright, timerange, [0, wfall.shape[1]]))
			wfall = driftrate.subtractbg(wfall, tleft, tright)
		subfall = driftrate.subsample(wfall, numf, numt)
	except (ValueError, ZeroDivisionError) as e:
		error_log_cb('subsample_cb', (numf, numt, e))
	else:
		gdata['wfall'] = subfall
		log_cb('subsample_cb', (numf, numt))
		updatedata_cb('subsample_cb', {'subsample': True})
		return subfall

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
	updatedata_cb(sender, {'fileidx': fileidx})
	log_cb(sender, 'Opening file {}'.format(gdata['files'][fileidx]))

def exportmask_cb(sender, data):
	datestr = datetime.now().strftime('%b%d')
	np.save('masks_{}.npy'.format(datestr), [gdata['masks']])
	print(gdata['masks'])

def importmask_cb(sender, data):
	if data is None:
		return
	if type(data) == list:
		data = data[-1]
	filename = data
	log_cb(sender, 'mask selected: {}'.format(data))
	if filename.split('.')[-1] == 'npy':
		masks = np.load(filename, allow_pickle=True)[0]
		if type(masks) == dict:
			gdata['masks'] = masks
			updatedata_cb(sender, {'keepview': True})
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
		updatedata_cb(sender, {'keepview': True})
		masktable_cb(sender, None)

def masktable_cb(sender, data):
	# dpg makes working with tables impossible so we will delete the table and re-add it every time
	dpg.delete_item('Masktable')

	mostmasks = 1
	for burst, masks in gdata['masks'].items():
		mostmasks = len(masks) if len(masks) > mostmasks else mostmasks
	tableheight = round(min(25*mostmasks, 250))
	dpg.add_table('Masktable', [], height=tableheight, parent='Masking', callback=removemask_cb)

	columns = [s.split('.')[0][-8:] for s in gdata['masks'].keys()]
	for key, col in zip(gdata['masks'].keys(), columns):
		dpg.add_column('Masktable', col, gdata['masks'][key])

def sksgmask_cb(sender, data):
	items = ['SKSGSigmaInput','SKSGWindowInput','SKSGStatus']
	if sender == 'EnableSKSGMaskBox':
		toggle_config(sender, {'kwargs': ['enabled'], 'items': items})

	sigma, window = dpg.get_value('SKSGSigmaInput'), dpg.get_value('SKSGWindowInput')
	df, dt      = gdata['burstmeta']['fres'], gdata['burstmeta']['tres']
	sksgmask = sk_sg_filter(
		data=gdata['wfall_original'].copy().T,
		foff=df,
		tsamp=dt/1000,
		spectral_kurtosis_sigma=sigma,
		savgol_frequency_window=window,
		savgol_sigma=sigma,
	)
	gdata['sksgmask'] = sksgmask
	updatedata_cb(sender, {'keepview': True})

def resulttable_cb(sender, data):
	coords = dpg.get_table_selections('Resulttable')
	newsel = coords[0].copy()
	for coord in coords:
		dpg.set_table_selection('Resulttable', coord[0], coord[1], False)
	print("resultidx = ", newsel[0])
	displayresult_cb('User', {'resultidx': newsel[0]})

def updateResultTable(resultsdf):
	dpg.delete_item('Resulttable')
	dpg.add_table('Resulttable', [], height=225, parent='ResultsGroup', callback=resulttable_cb)

	# subset of driftrate.columns:
	columns = ['name', 'DM', 'amplitude', 'slope (mhz/ms)', 'theta', 'center_f']
	columns = ['name', 'DM', 'amplitude', 'tsamp_width','subbg_start (ms)', 'subbg_end (ms)']
	dpg.set_headers('Resulttable', columns)

	# [burstname, trialDM, center_f, slope, slope_err, theta, red_chisq], popt, perr, [fres_MHz, tres_ms/1000]
	for burstname, rowdata in resultsdf.iterrows():
		newrow = [burstname] + [rowdata[col] for col in columns[1:]]
		dpg.add_row('Resulttable', newrow)

def mousemask_cb(sender, data):
	isOnWaterfall = dpg.is_item_hovered('WaterfallPlot')
	if isOnWaterfall:
		tchan, fchan = dpg.get_plot_mouse_pos()
		rawmask = round(fchan)

		# map frequency (rawmask) to channel number
		spectralrange = [round(gdata['extents'][2]), round(gdata['extents'][3])]
		mask = np.interp(rawmask, spectralrange, [0, gdata['wfall_original'].shape[0]])
		mask = int(mask)

		if mask not in gdata['masks'][gdata['currfile']]:
			gdata['masks'][gdata['currfile']].append(mask)

		updatedata_cb(sender, {'keepview': True})
		masktable_cb(sender, None)
		log_cb('mousemask_cb ', [[tchan, fchan], isOnWaterfall])
	else:
		return

def dmrange_cb(sender, data):
	dmrange   = dpg.get_value('dmrange')
	dmrange   = np.round(dmrange, 4) # dpg.get_value introduces rounding errors
	numtrials = dpg.get_value('numtrials')
	dmstep = round(dpg.get_value('dmstep'), 3)
	burstDM = gdata['burstDM']
	if dmrange[1] < dmrange[0]:
		dmrange.sort()
		dpg.set_value('dmrange', dmrange)
	if not (dmrange[0] < burstDM < dmrange[1]):
		dpg.configure_item('DMWarning', show=True)
	else:
		dpg.configure_item('DMWarning', show=False)

	if sender == 'dmstep' or sender == 'user':
		numtrials = round((dmrange[1] - dmrange[0])/dmstep)+1
		dpg.set_value('numtrials', numtrials)
	elif sender == 'numtrials':
		dmstep = round((dmrange[1] - dmrange[0])/numtrials, 3)
		dpg.set_value('dmstep', dmstep)
	gdata['trialDMs'] = np.arange(dmrange[0], dmrange[1]+dmstep, step=dmstep)

def getCurrentBurst():
	""" Return the currently loaded waterfall at its burst DM """
	burstname = dpg.get_value('burstname').replace(',', '')
	df, dt = gdata['burstmeta']['fres'], gdata['burstmeta']['tres']
	lowest_freq = min(gdata['burstmeta']['dfs']) # mhz
	wfall = gdata['wfall'].copy()
	wfall_cr = driftrate.cropwfall(wfall, twidth=dpg.get_value('twidth'))
	burstDM = gdata['burstDM']
	return burstname, df, dt, lowest_freq, wfall_cr, burstDM

def getMeasurementInfo(wfall_cr):
	# TODO: region info, raw shape
	tsamp_width = dpg.get_value('twidth')
	fchans, tchans = wfall_cr.shape
	subbgstart, subbgend, sksigma, skwindow = None, None, None, None
	if dpg.get_value('EnableSubBGBox'):
		subbgstart, subbgend = dpg.get_value('SubtractBGRegion')
	if dpg.get_value('EnableSKSGMaskBox'):
		sksigma, skwindow = dpg.get_value('SKSGSigmaInput'), dpg.get_value('SKSGWindowInput')
	cols = ['fchans', 'tchans', 'tsamp_width','subbg_start (ms)', 'subbg_end (ms)','sksigma','skwindow']
	row = [fchans, tchans, tsamp_width, subbgstart, subbgend, sksigma, skwindow]
	return cols, row

def slope_cb(sender, data):
	burstname, df, dt, lowest_freq, wfall_cr, burstDM = getCurrentBurst()

	if data is not None:
		dpg.set_value('SlopeStatus', 'Status: Doing second pass...')
		p0 = data['p0']
		trialDMs = data['badfitDMs']
	else:
		dpg.set_value('SlopeStatus', 'Status: Calculating...')
		p0 = []
		trialDMs = np.unique(np.append(gdata['trialDMs'], burstDM))

	results, burstdf = driftrate.processDMRange(burstname, wfall_cr, burstDM, trialDMs, df, dt, lowest_freq, p0)

	if gdata['multiburst']['enabled']:
		subbursts = getSubbursts()
		subresults, subdf = [], pd.DataFrame()
		for subname, subburst in subbursts.items():
			print('processing {}'.format(subname))
			ret, retdf = driftrate.processDMRange(subname, subburst, burstDM, trialDMs, df, dt, lowest_freq)
			subresults.append(ret)
			subdf = subdf.append(retdf)

		burstdf = burstdf.append(subdf)

	# Do a second pass using the best p0 just found
	# TODO: only repeat bad fits
	p0 = getOptimalFit(burstdf)
	if data is None:
		# ensure p0 is in center of autocorr
		if not (wfall_cr.shape[1]*0.9 <= float(p0[1]) <= wfall_cr.shape[1]*1.1):
			p0[1] = wfall_cr.shape[1]
		if not (wfall_cr.shape[0]*0.9 <= float(p0[2]) <= wfall_cr.shape[0]*1.1):
			p0[2] = wfall_cr.shape[0]
		# todo: find bad dms
		return slope_cb(sender, {'p0' : p0, 'badfitDMs': trialDMs})

	# Add measurement info to row (things needed to reproduce/reload the measurement)
	cols, row = getMeasurementInfo(wfall_cr)
	burstdf[cols] = row

	# Save results
	if gdata['resultsdf'] is None:
		gdata['resultsdf'] = burstdf
	else:
		# overwrite if there are already results
		if gdata['resultsdf'].index.str.startswith(burstname).any():
			df = gdata['resultsdf']
			gdata['resultsdf'] = df.drop(df[df.index.str.startswith(burstname)].index)
		gdata['resultsdf'] = gdata['resultsdf'].append(burstdf)
	backupresults()
	print(gdata['resultsdf'])

	# burstdf = burstdf.loc[burstname] # remove??
	gdata['burstdf'] = burstdf

	dpg.set_value('SlopeStatus', 'Status: Done.')
	dpg.set_value('NumMeasurementsText', "# of Measurements for this burst: {}".format(len(burstdf)))
	dpg.set_value('TotalNumMeasurementsText', "Total # of Measurements: {}".format(len(gdata['resultsdf'])))
	dpg.configure_item('PrevDM', enabled=True)
	dpg.configure_item('NextDM', enabled=True)
	dpg.configure_item('ExportCSVBtn', enabled=True)
	dpg.configure_item('ExportPDFBtn', enabled=True)
	dpg.configure_item('EnableP0Box', enabled=True)
	initializeP0Group()
	updateResultTable(burstdf)
	plotdata_cb(sender, data)

def redodm_cb(sender, data):
	p0 = gdata['p0']
	useForAll = dpg.get_value('P0AllDMsBox')

	burstname, df, dt, lowest_freq, wfall_cr, burstDM = getCurrentBurst()
	displayedname = gdata['displayedBurst']
	if gdata['multiburst']['enabled'] and displayedname != burstname:
		subbursts = getSubbursts()
		wfall_cr = subbursts[displayedname]

	print('redoing ', displayedname, gdata['displayedDM'])
	result, burstdf = driftrate.processDMRange(
		displayedname, wfall_cr, burstDM, [float(gdata['displayedDM'])],
		df, dt, lowest_freq, p0=p0
	)
	cols, row = getMeasurementInfo(wfall_cr)
	burstdf[cols] = row
	df = gdata['resultsdf']
	df[(df.index == displayedname) & (df['DM'] == gdata['displayedDM'])] = burstdf

	gdata['resultsdf'] = df
	gdata['burstdf'] = gdata['resultsdf'].loc[burstname]

	if gdata['multiburst']['enabled']:
		subburstdf = gdata['resultsdf'][gdata['resultsdf'].index.str.startswith(burstname)]
		updateResultTable(subburstdf)
		subburstdf = subburstdf.reset_index()
		dispdm = gdata['displayedDM']
		resultidx = subburstdf[(subburstdf.name == displayedname) & (subburstdf.DM == dispdm)].index[0]
		if data is None:
			data = {}
		data['resultidx'] = resultidx
	else:
		updateResultTable(gdata['burstdf'])
	backupresults()
	plotdata_cb(sender, data)

def getAllRegions():
	ret = {}
	suffixes = itertools.cycle(subburst_suffixes)
	if gdata['multiburst']['enabled']:
		for regid in range(1, gdata['multiburst']['numregions']):
			region = dpg.get_value('Region{}'.format(regid))
			regiontype = dpg.get_value('RegionType{}'.format(regid))
			if regiontype == 1: # burst
				suffix = next(suffixes)
				ret[suffix] = region
			else:
				ret['background'] = region
	return ret

def loadresults_cb(sender, data):
	resultsfile = os.path.join(*data)
	resultsdf = pd.read_csv(resultsfile).set_index('name')
	gdata['resultsdf'] = resultsdf
	burstselect_cb(sender, data)

def exportresults_cb(sender, data):
	resultsdf = gdata['resultsdf']

	df = driftlaw.computeModelDetails(resultsdf)

	datestr = datetime.now().strftime('%b%d')
	prefix = dpg.get_value('ExportPrefix')
	filename = '{}_results_{}rows_{}.csv'.format(prefix, len(df.index), datestr)
	try:
		df.to_csv(filename)
		dpg.set_value('ExportCSVText', 'Saved to {}'.format(filename))
	except PermissionError as e:
		dpg.set_value('ExportCSVText', 'Permission Denied')
	dpg.configure_item('ExportCSVText', show=True)

	if sender == 'ExportPDFBtn':
		dpg.configure_item('ExportPDFText', show=True)
		dpg.set_value('ExportPDFText', 'Saving PDF...')
		success = driftrate.plotResults(filename, datafiles=gdata['files'], masks=gdata['masks'], regionsfile=gdata['multiburst']['regions'])
		if success:
			dpg.set_value('ExportPDFText', f'Saved to {filename.split(".")[0]+".pdf"}')
		else:
			dpg.set_value('ExportPDFText', f'Permission Denied')

def backupresults():
	resultsdf = gdata['resultsdf']
	df = driftlaw.computeModelDetails(resultsdf)
	datestr = datetime.now().strftime('%b%d')
	prefix = 'backup'
	filename = '{}_results_{}.csv'.format(prefix, datestr)
	df.to_csv(filename)

def displayresult_cb(sender, data):
	burstDM = gdata['burstDM']
	trialDMs = np.unique(np.append(gdata['trialDMs'], burstDM))
	dmlist = list(trialDMs)

	burstname = dpg.get_value('burstname').replace(',', '')
	subname = gdata['displayedBurst']
	dispdm = gdata['displayedDM']
	subburstdf = gdata['resultsdf'][gdata['resultsdf'].index.str.startswith(burstname)].reset_index()
	resultidx = subburstdf[(subburstdf.name == subname) & (np.isclose(subburstdf.DM, dispdm))].index[0]
	subburstdf = subburstdf.set_index('name')

	if sender == 'User':
		resultidx = data['resultidx']
	elif sender == 'NextDM':
		resultidx = resultidx + 1
		if not (resultidx < len(subburstdf)):
			resultidx = 0
	elif sender == 'PrevDM':
		resultidx = resultidx - 1

	if data is None:
		data = {}
	data['resultidx'] = resultidx
	subname = subburstdf.index[resultidx]
	gdata['displayedDM'] = dmlist[resultidx % len(dmlist)]
	dpg.set_value('dmdisplayed', str(round(gdata['displayedDM'], getscale(gdata['displayedDM']))))
	dpg.set_value('burstdisplayed', subname)
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
			df = gdata['resultsdf'].loc[~burstname] # use df.drop see slope_cb
		elif data == 'all':
			gdata['resultsdf'] = None
	else:
		confirmpopup(data, deleteresults_cb)

	print(sender, data)

def getOptimalFit(df):
	params = ['amplitude', 'xo', 'yo', 'sigmax', 'sigmay', 'angle']
	bestrow = df[df.red_chisq == df.red_chisq.min()]
	if len(bestrow) > 1:
		bestrow = bestrow.head(1)
	p0 = [bestrow[param] for param in params]
	return p0

def initializeP0Group():
	p0 = []
	burstname = dpg.get_value('burstname').replace(',', '')
	subburstdf = gdata['resultsdf'][gdata['resultsdf'].index.str.startswith(burstname)]
	wfall_cr = getCurrentBurst()[-2]
	p0 = getOptimalFit(subburstdf)
	p0f = [p0i.iloc[0] for p0i in p0]
	if p0f[0] < 0:
		p0f = [-p0fi if p0fi < 0 else p0fi for p0fi in p0f]

	gdata['p0'] = p0
	dpg.set_value("AmplitudeDrag", p0f[0])
	dpg.set_value("AngleDrag", p0f[5])
	# solution will always be near the center
	dpg.set_value("x0y0Drag", [wfall_cr.shape[1], wfall_cr.shape[0]])
	dpg.set_value("SigmaXYDrag", [p0f[3], p0f[4]])
	for item in ['AmplitudeDrag', 'AngleDrag', 'x0y0Drag', 'SigmaXYDrag']:
		val = dpg.get_value(item)
		if type(val) == list:
			val = val[0]
		dpg.configure_item(item, speed=1/10**getscale(val), format='%.{}f'.format(getscale(val)+1))
	return p0

def enablep0_cb(sender, data):
	toggle_config(sender, data)
	updatep0_cb(sender, data)

def updatep0_cb(sender, data):
	p0 = []

	# get resultidx
	subname, dispdm = gdata['displayedBurst'], gdata['displayedDM']
	burstname = dpg.get_value('burstname').replace(',', '')
	subburstdf = gdata['resultsdf'][gdata['resultsdf'].index.str.startswith(burstname)].reset_index()
	resultidx = subburstdf[(subburstdf.name == subname) & (subburstdf.DM == dispdm)].index[0]
	if gdata['multiburst']['enabled']:
		if data is None:
			data = {}
		data['resultidx'] = resultidx

	for item in ['AmplitudeDrag', 'x0y0Drag', 'SigmaXYDrag', 'AngleDrag']:
		val = dpg.get_value(item)
		if type(val) == list:
			p0 = p0 + val
		else:
			p0.append(val)
	gdata['p0'] = p0
	plotdata_cb(sender, data)

def enablesplitting_cb(sender, data):
	items = ['Region', 'RegionType', 'RemoveRegionBtn', 'AddRegionBtn', 'ExportRegionBtn']
	gdata['multiburst']['enabled'] = dpg.get_value(sender)
	items = items.copy()
	items.remove('AddRegionBtn')
	items.remove('ExportRegionBtn')
	toggle_config(sender, {'kwargs': ['enabled'], 'items': ['AddRegionBtn']})
	toggle_config(sender, {'kwargs': ['enabled'], 'items': ['ExportRegionBtn']})
	for regid in range(1, gdata['multiburst']['numregions']):
		if dpg.does_item_exist('RegionSelector{}'.format(regid)):
			itemsid = list(map(lambda item: item+'{}'.format(regid), items))
			toggle_config(sender, {'kwargs': ['enabled'], 'items': itemsid})

			if gdata['multiburst']['enabled']:
				drawregion_cb([regid], None)
			else:
				seriesname = 'Region{}Series'.format(regid)
				[dpg.delete_series(plot, seriesname) for plot in ['WaterfallPlot', 'TimeSeriesPlot']]

def exportregions_cb(sender, data):
	regions = getAllRegions()
	saveobj = {gdata['displayedBurst']: regions}
	datestr = datetime.now().strftime('%b%d')
	filename = 'burstregions_{}.npy'.format(datestr)
	if os.path.exists(filename):
		loadobj = np.load(filename, allow_pickle=True)[0]
		loadobj[gdata['displayedBurst']] = regions
		saveobj = loadobj
	np.save(filename, [saveobj])
	print(saveobj)
	print('Saved', filename)

def importregions_cb(sender, data):
	if data is None:
		return
	filename = data
	if os.path.exists(filename):
		loadobj = np.load(filename, allow_pickle=True)[0]
		gdata['multiburst']['regions'] = loadobj

def addregion_cb(sender, data):
	regionSelector()
	drawregion_cb([gdata['multiburst']['numregions']-1], None)

def removeregion_cb(sender, data):
	regid = sender[-1]
	dpg.delete_item('RegionSelector{}'.format(regid))
	seriesname = 'Region{}Series'.format(regid)
	[dpg.delete_series(plot, seriesname) for plot in ['WaterfallPlot', 'TimeSeriesPlot']]
	gdata['multiburst']['numregions'] -= 1

colors = [
	(0, 255, 0),
	(255, 0, 0),
	(0, 0, 255),
	(255, 255, 0),
	(255, 0, 255),
	(0, 255, 255)
]
def drawregion_cb(sender, data):
	regid = sender[-1]
	region = dpg.get_value('Region{}'.format(regid))
	regiontype = dpg.get_value('RegionType{}'.format(regid))

	if region[1] < region[0]:
		region.sort()
		dpg.set_value('Region{}'.format(regid), region)

	seriesname = 'Region{}Series'.format(regid)
	for plot in ['WaterfallPlot', 'TimeSeriesPlot']:
		dpg.delete_series(plot, seriesname)
		dpg.add_vline_series(plot, seriesname, region,
			update_bounds=False,
			weight=1.5,
			color=colors[(int(regid)-1) % len(colors)]
		)

def regionSelector():
	regid = gdata['multiburst']['numregions']
	before = "AddRegionBtn" if regid > 1 else ""
	enabled = gdata['multiburst']['enabled']
	if 'wfall' in gdata:
		maxval =  gdata['extents'][1]
	else:
		maxval = 100

	with dpg.group('RegionSelector{}'.format(regid), horizontal=True, parent='SplittingSection',
		before=before):
		dpg.add_drag_float2('Region{}'.format(regid),
			label='(ms)',
			width=280,
			enabled=enabled,
			max_value=maxval,
			speed=maxval*0.005,
			default_value=[0, maxval],
			callback=drawregion_cb
		)
		helpmarker("double click to edit")
		dpg.add_radio_button('RegionType{}'.format(regid), items=["Background", "Burst"],
			horizontal=True,
			callback=drawregion_cb,
			enabled=enabled
		)
		dpg.add_button('RemoveRegionBtn{}'.format(regid), label='X', enabled=enabled, callback=removeregion_cb)
	gdata['multiburst']['numregions'] += 1

subburst_suffixes = driftrate.subburst_suffixes
def getSubbursts():
	if not gdata['multiburst']['enabled']:
		return []

	burstname, df, dt, lowest_freq, wfall_cr, burstDM = getCurrentBurst()
	background = None
	subbursts  = []
	for regid in range(1, gdata['multiburst']['numregions']):
		if dpg.does_item_exist('RegionSelector{}'.format(regid)):
			regionname = 'Region{}'.format(regid)
			typename = 'RegionType{}'.format(regid)

			region = dpg.get_value(regionname)
			trange = driftrate.getExtents(wfall_cr, df=df, dt=dt, lowest_freq=lowest_freq)[0][:2]
			for i, edge in enumerate(region):
				region[i] = round(np.interp(edge, trange, [0, wfall_cr.shape[1]]))

			regiontype =  dpg.get_value(typename)
			if regiontype == 0:   # Background
				background = wfall_cr[:, region[0]:region[1]]
			elif regiontype == 1: # Burst
				subburst = wfall_cr[:, region[0]:region[1]]
				subbursts.append(subburst)

	if background is None:
		error_log_cb("getSubbursts", "Please specify a background region")
		raise "Please specify a background region"

	# pad with background
	# for subburst in subbursts:
	subburstsobj = {}
	for subburst, suffix in zip(subbursts, subburst_suffixes):
		subburst = np.concatenate((background, subburst, background), axis=1)
		subname = burstname +'_'+ suffix # '_' is used in plotResults to split names
		subburstsobj[subname] = subburst
	return subburstsobj

def subtractbg_cb(sender, data):
	if sender == 'EnableSubBGBox':
		toggle_config(sender, {'kwargs': ['enabled'], 'items': ['SubtractBGRegion']})
		maxval = gdata['extents'][1]
		dpg.configure_item('SubtractBGRegion', max_value=maxval)
		if dpg.get_value('SubtractBGRegion')[1] > maxval:
			dpg.set_value('SubtractBGRegion', [0, maxval*0.2])

	subsample_cb(sender, data)

def helpmarker(message):
	dpg.add_same_line()
	dpg.add_text("(?)", color=[150, 150, 150], tip=message)

def toggle_config(sender, data):
	config_dict = {}
	for kwarg in data['kwargs']:
		config_dict[kwarg] = dpg.get_value(sender)
	for item in data['items']:
		dpg.configure_item(item, **config_dict)

def frbgui(filefilter=gdata['globfilter'],
		datadir=None,
		maskfile=None,
		regionfile=None,
		dmrange=None,
		numtrials=10,
		dmstep=0.1,
		winwidth=1700,
		winheight=850,
	):
	gdata['datadir'] = datadir
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
			dpg.add_input_int('twidth', label='Display width (# chans)',
				default_value=twidth_default,
				step=10,
				callback=twidth_cb
			)
			dpg.add_input_float('twidth_ms', label='Display width (ms)', enabled=False)

			dpg.add_text('Dedispersion Range for all Bursts: ')
			dpg.add_text('DMWarning', default_value='Warning: Range chosen does not include burst DM',
				color=[255, 0, 0], show=False)
			dpg.add_drag_float2('dmrange', label='DM range (pc/cm^3)', callback=dmrange_cb,
				min_value=0, max_value=0)
			helpmarker('Double click to edit')
			dpg.add_input_int('numtrials', label='# of Trial DMs', default_value=10, callback=dmrange_cb)
			dpg.add_text(' or ')
			dpg.add_input_float('dmstep', label='DM Step (pc/cm^3)', default_value=0.1,
				min_value=0.001,
				min_clamped=True,
				callback=dmrange_cb)

		with dpg.collapsing_header("2. Waterfall Cleanup", default_open=True):
			with dpg.tree_node('Subtract Background', default_open=True):
				dpg.add_checkbox('EnableSubBGBox', label='Subtract background sample', callback=subtractbg_cb)
				dpg.add_drag_float2('SubtractBGRegion',
					label='t_start (ms), t_end (ms)',
					width=280,
					enabled=False,
					max_value=100,
					speed=0.5,
					default_value=[0, 25],
					callback=subtractbg_cb
				)

			with dpg.tree_node('Masking', default_open=True):
				dpg.add_text("Click on the waterfall plot to begin masking frequency channels.")
				dpg.add_text("NOTE: only mask on the original waterfall (todo: add a 'mask' button)")

				dpg.add_button('Export Masks', callback=exportmask_cb, enabled=True)
				dpg.add_same_line()
				dpg.add_button('Import Masks',
					callback=lambda s, d: dpg.open_file_dialog(importmask_cb, extensions='.npy'),
					enabled=True)

				dpg.add_table('Masktable', [], height=50)

			with dpg.tree_node('Auto Masking', default_open=True):
				dpg.add_checkbox('EnableSKSGMaskBox', label='Enable SK-SG Filter', default_value=False,
					callback=sksgmask_cb)
				dpg.add_input_int("SKSGSigmaInput", width=100, default_value=3, label="sigma",
					callback=sksgmask_cb, enabled=False, min_value=1)
				dpg.add_same_line()
				dpg.add_input_int("SKSGWindowInput", width=100, default_value=15, label="window",
					callback=sksgmask_cb, enabled=False, min_value=1)
				dpg.add_text('SKSGStatus', default_value=' ')

			with dpg.tree_node('Downsampling', default_open=True):
				dpg.add_text("Wfallshapelbl", default_value="Original Size: (no burst selected)")
				dpg.add_text("Subfallshapelbl", default_value="Current Size: (no burst selected)")
				dpg.add_input_int("numfreqinput", width=100, label="numf", callback=subsample_cb, enabled=False)
				dpg.add_same_line()
				dpg.add_input_int("numtimeinput", width=100, label="numt", callback=subsample_cb, enabled=False)
				dpg.add_button('ResetSamplingBtn', label='Reset', callback=subsample_cb, enabled=False)

		with dpg.collapsing_header("SplittingSection", label="3. Burst Splitting", default_open=True):
			dpg.add_checkbox('MultiBurstBox', label='Are there multiple bursts in this waterfall?',
				default_value=False, enabled=True,
				callback=enablesplitting_cb,
				callback_data={'kwargs': ['enabled']}
			)
			regionSelector()
			dpg.add_button('AddRegionBtn', label="Add Region", callback=addregion_cb, enabled=False)
			dpg.add_same_line()
			dpg.add_button('ExportRegionBtn', label="Export Regions", callback=exportregions_cb, enabled=False)

		with dpg.collapsing_header('SlopeSection', label="4. Sub-burst Slope Measurements", default_open=True):
			dpg.add_button("Measure Slope over DM Range", callback=slope_cb)
			dpg.add_same_line()
			dpg.add_text("SlopeStatus", default_value="Status: (click 'Measure Slope' to calculate)")
			dpg.add_button("Load Results",
				callback=lambda s, d:dpg.open_file_dialog(loadresults_cb, extensions='.csv'))
			dpg.add_same_line()
			dpg.add_text("LoadResultsWarning", color=(255, 0, 0),
				default_value="Warning: will overwrite unsaved results")
			dpg.add_button('DeleteBurstBtn',
				label="Delete this burst's results",
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
				dpg.add_button("PrevDM", arrow=True, direction=dpg.mvDir_Left, enabled=False,
					callback=displayresult_cb)
				dpg.add_button("NextDM", arrow=True, direction=dpg.mvDir_Right, enabled=False,
					callback=displayresult_cb)
				dpg.add_text("Burst Displayed: ")
				dpg.add_text("burstdisplayed", default_value=str(dpg.get_value('burstname')))
				dpg.add_text("DM Displayed (pc/cm^3): ")
				dpg.add_text("dmdisplayed", default_value=str(dpg.get_value('DM')))
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

			dpg.add_input_text('ExportPrefix', label='Filename Prefix', default_value="Test")
			dpg.add_button('ExportCSVBtn', label="Export Results CSV", callback=exportresults_cb, enabled=False)
			dpg.add_same_line()
			dpg.add_text('ExportCSVText', default_value=' ', show=True, color=(0, 255, 0))
			dpg.add_button('ExportPDFBtn', label="Export Results PDF", callback=exportresults_cb, enabled=False)
			dpg.add_same_line()
			dpg.add_text('ExportPDFText', default_value=' ', show=True, color=(0, 255, 0))


	### Plotting window
	with dpg.window("FRB Plots", width=1035, height=745, x_pos=600, y_pos=30):
		dpg.set_mouse_click_callback(mousemask_cb)
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
	dpg.set_main_window_size(winwidth, winheight)
	dpg.set_main_window_title("FRB Repeaters")
	# dpg.start_dearpygui()

	# Load defaults
	dpg.set_value('Filter', filefilter)
	directory_cb('user', [datadir])
	burstselect_cb('burstselect', None)
	importmask_cb('user', maskfile)
	importregions_cb('user', regionfile)

	## dm range defaults
	dpg.set_value('dmrange', dmrange)
	dpg.set_value('numtrials', numtrials)
	dpg.set_value('dmstep', dmstep)
	dmrange_cb('user', None)

	# loadresults_cb('user2', ['B:\\dev\\frbrepeaters\\results', 'FRB121102_oostrum_525rows_Aug28.csv'])

	dpg.start_dearpygui(primary_window='main') # blocks til closed
	return gdata

if __name__ == '__main__':
	frbgui(
		datadir='B:\\dev\\frbrepeaters\\data\\gajjar2018',
		# datadir='B:\\dev\\frbrepeaters\\data\\aggarwal2021',
		# datadir='B:\\dev\\frbrepeaters\\data\\oostrum2020\\R1_frb121102',
		# maskfile='aggarwalmasks_sept12.npy',
		regionfile='burstregions_gajjaroct1.npy',
		dmstep=5,
		dmrange=[555, 575]
	)
