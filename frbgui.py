import dearpygui.demo as demo
import dearpygui.dearpygui as dpg
from dearpygui_ext import logger, themes
import driftrate, driftlaw
import os, glob, itertools, io
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
from datetime import datetime
import warnings
from your.utils.rfi import sk_sg_filter
import time
warnings.filterwarnings("ignore")

defaultmpl_backend = matplotlib.get_backend()
twidth_default = 150
logwin = None
# flag for switching the main thread from DPG to Matplotlib. Matplotlib is not threadsafe
exportPDF = False

# GUI data is stored in this object. Defaults initialized here and at the bottom
gdata = {
	'globfilter'     : '*.npz',
	'masks'          : {},                   # will store masks, either lists of channels or ranges
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

def compressArr(arr):
	""" adapted from https://www.geeksforgeeks.org/compress-the-array-into-ranges/ """
	i, j, n = 0, 0, len(arr)
	arr.sort()
	while i < n:
		j = i
		while (j + 1 < n) and (arr[j + 1] == arr[j] + 1):
			j += 1

		if i == j:
			print(arr[i], end=" ")
			i += 1
		else:
			print(arr[i], "-", arr[j], end=" ")
			i = j + 1

def applyMasks(wfall):
	for mask in gdata['masks'][gdata['currfile']]['chans']:
		if mask < len(wfall):
			wfall[mask] = 0

	for rangeid, maskrange in gdata['masks'][gdata['currfile']]['ranges'].items():
		masks = np.arange(maskrange[0], maskrange[1]+1)
		for mask in masks:
			if mask < len(wfall):
				wfall[mask] = 0

	if 'sksgmask' in gdata and dpg.get_value('EnableSKSGMaskBox'):
		wfall[gdata['sksgmask'], :] = 0

	return wfall

def makeburstname(filename):
	return filename.split(os.sep)[-1].split('.')[0]

def log_cb(sender, data):
	print(f"{sender}, {data}")
	logwin.log_debug(f"{sender}, {data}")
def error_log_cb(sender, data):
	print(f"{sender}, {data}")
	logwin.log_error(f"{sender}, {data}")

def updatedata_cb(sender, data, udata):
	if not data: data = {}
	if not udata: udata = {}
	wfall = None
	if 'filename' in udata.keys(): # load burst from disk
		filename = udata['filename']
		gdata['currfile'] = filename
		if gdata['currfile'] not in gdata['masks'].keys():
			gdata['masks'][gdata['currfile']] = {'chans':[], 'ranges':{}}
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
						# print(df, loaded['bandwidth'], storedshape[1]*1000)
						gdata['burstmeta']['fres_original'] = df
						gdata['burstmeta']['fres'] = df
						dpg.set_value('df', df)
						dpg.configure_item('df', format='%.{}f'.format(getscale(df)+1))
					elif key == 'duration' or key == 'dt':
						dt = loaded['duration'] / storedshape[1]*1000
						# print(dt, loaded['duration'], storedshape[1]*1000)
						gdata['burstmeta']['duration'] = loaded['duration']
						gdata['burstmeta']['tres_original'] = dt
						gdata['burstmeta']['tres'] = dt
						dpg.set_value('dt', dt)
						dpg.set_value('duration', loaded['duration'])
						dpg.configure_item(key, format='%.{}f'.format(getscale(dt)+3))
					else:
						if dpg.does_alias_exist(key):
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
		dpg.configure_item('SaveDMButton', enabled=False)

		# update mask range controls
		gdata['maskrangeid'] = 0
		dpg.delete_item("MaskRangeGroup", children_only=True)
		if gdata['masks'][gdata['currfile']]['ranges']:
			rangeids = []
			for rangeid, maskrange in gdata['masks'][gdata['currfile']]['ranges'].items():
				rangeids.append(rangeid)
				addMaskRange_cb(sender, rangeid)
				dpg.set_value(f'rangeslider##{rangeid}', maskrange)
			gdata['maskrangeid'] = max(rangeids) # guarantee no collisions

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
				gdata['burstdf'] = gdata['resultsdf'][gdata['resultsdf'].index.str.startswith(burstname)]

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

					skitems = ['SKSGSigmaInput','SKSGWindowInput']
					if not pd.isnull(sksigma):
						dpg.set_value("EnableSKSGMaskBox", True)
						dpg.set_value('SKSGSigmaInput', int(sksigma))
						dpg.set_value('SKSGWindowInput', int(skwindow))
						toggle_config('EnableSKSGMaskBox', {'kwargs': ['enabled'], 'items': skitems})
					else:
						dpg.set_value('EnableSKSGMaskBox', False)
						toggle_config('EnableSKSGMaskBox', {'kwargs': ['enabled'], 'items': skitems})

				# check if subsampling is needed (for eg. if results have been loaded)
				dmframe = gdata['burstdf'].set_index('DM')
				downf, downt = map(int, dmframe.loc[np.isclose(dmframe.index, gdata['displayedDM'])].iloc[0][['downf', 'downt']])
				dpg.set_value('numfreqinput', int(wfall.shape[0]/downf))
				dpg.set_value('numtimeinput', int(wfall.shape[1]/downt))
				if downf != 1 or downt != 1:
					wfall = subsample_cb('loadresults_cb', None) # this updates gdata['burstmeta']
				dpg.set_value('twidth', int(twidth)) # force twidth to match results row after subsampling

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
				addregion_cb(sender, None)
			while numregions < gdata['multiburst']['numregions']-1:
				removeregion_cb([gdata['multiburst']['numregions']-1], None)
			# set the elements based on gdata['regions']
			for regid, name in enumerate(regions):
				regid += 1 # off by one
				regiontype = 0 if 'background' in name else 1 # 0 is background, 1 is burst
				dpg.set_value('Region{}'.format(regid), regions[name])
				dpg.set_value('RegionType{}'.format(regid), regiontype)
				drawregion_cb([regid], None)

	elif sender == 'subsample_cb' and udata['subsample']: # ie. sender == 'subsample_cb' dpg.get_value('DM')
		wfall = gdata['wfall']

		disp_bandwidth = gdata['extents'][3] - gdata['extents'][2] # == gdata['burstmeta']['bandwidth']
		disp_duration = gdata['extents'][1] - gdata['extents'][0]
		twidth = disp_duration/gdata['burstmeta']['tres_original']/2
		dpg.set_value('twidth', round(twidth * (wfall.shape[1]/gdata['wfall_original'].shape[1])))
		dpg.configure_item('twidth', max_value=round(wfall.shape[1]/2))

		# wfall_cr = driftrate.cropwfall(wfall, twidth=dpg.get_value('twidth'))
		gdata['burstmeta']['fres'] = gdata['burstmeta']['bandwidth'] / wfall.shape[0]
		gdata['burstmeta']['tres'] = gdata['burstmeta']['duration']*1000 / wfall.shape[1]

		dpg.set_value('Subfallshapelbl', 'Current Size: {}'.format(np.shape(wfall)))
	else:
		wfall = gdata['wfall']

	if wfall.shape == gdata['wfall_original'].shape:
		wfall = applyMasks(np.copy(gdata['wfall_original']))
		if dpg.get_value('EnableSubBGBox'):
			tleft, tright, _, _ = dpg.get_value('SubtractBGRegion')
			timerange = [round(gdata['extents'][0]), round(gdata['extents'][1])]
			tleft  = round(np.interp(tleft, timerange, [0, wfall.shape[1]]))
			tright = round(np.interp(tright, timerange, [0, wfall.shape[1]]))
			wfall = driftrate.subtractbg(wfall, tleft, tright)


	gdata['wfall'] = wfall
	# gdata['ts']    = np.nanmean(wfall, axis=0) # time series at burstDM
	# gdata['pkidx'] = np.nanargmax(gdata['ts']) # pkidx at burstDM, for displaying across DMs

	plotdata_cb(sender, data, udata)

# plt.figure will fail unless it is on the main thread, and DPG runs callbacks on a seperate thread
matplotlib.use('agg') # this is a workaround to the thread issue
def getcorr2dtexture(corr, popt=None, p0=None, extents=None, slope=None, clim=None):
	plt.figure(figsize=(5, 5))

	plt.imshow(corr, origin='lower', interpolation='none', aspect='auto', cmap='gray',
				extent=extents)
	plt.clim(0, np.max(corr)/20)
	if clim and clim > 0:
		plt.clim(0, clim)
	if popt is not None and popt[0] > 0:
		fitmap = driftrate.makeDataFitmap(popt, corr, extents)
		if 1 not in fitmap.shape: # subsample ui can result in fitmaps unsuitable for countour plotting
			plt.contour(fitmap, [popt[0]/4, popt[0]*0.9], colors='b', alpha=0.75, extent=extents)
			if slope is not None:
				print(f"{slope = } {list(popt) = }")
				xo = popt[1]
				x = np.array([extents[2] / slope + xo, extents[3] / slope + xo])
				plt.plot(x, slope*(x-xo), 'g--')

	if p0 is not None and p0[0] > 0:
		fitmap = driftrate.makeDataFitmap(p0, corr, extents)
		plt.contour(fitmap, [p0[0]/4, p0[0]*0.9], colors='g', alpha=0.75, origin='lower',
					extent=extents)
	plt.xlim(extents[:2])
	plt.ylim(extents[2:])

	# remove axes, whitespace
	plt.gca().set_axis_off()
	plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
	plt.margins(0,0)
	plt.gca().xaxis.set_major_locator(plt.NullLocator())
	plt.gca().yaxis.set_major_locator(plt.NullLocator())

	fig = plt.gcf()
	fig.canvas.draw()
	texture = np.frombuffer(fig.canvas.tostring_argb(), dtype='uint8')
	w, h = fig.get_size_inches()*fig.dpi
	texture = texture.reshape(int(w*h), 4)[:, [1, 2, 3, 0]].flatten()/255 # reshape to rgba, and normalize
	plt.close('all')
	return texture, int(w), int(h)

def plotdata_cb(sender, data, userdata):
	if not data: data = {}
	if not userdata: userdata = {}

	df, dt      = gdata['burstmeta']['fres'], gdata['burstmeta']['tres']
	lowest_freq = min(gdata['burstmeta']['dfs']) # mhz
	ddm = gdata['displayedDM'] - gdata['burstDM']
	burstname = dpg.get_value('burstname').replace(',', '')

	popt, slope = None, None
	poptcols = ['amplitude', 'xo', 'yo', 'sigmax', 'sigmay', 'angle']
	if gdata['resultsdf'] is not None:
		subburstdf = gdata['resultsdf'][gdata['resultsdf'].index.str.startswith(burstname)]
		if not subburstdf.empty:
			dmframe = subburstdf.loc[burstname].set_index('DM')
			popt = dmframe.loc[np.isclose(dmframe.index, gdata['displayedDM'])].iloc[0][poptcols]
			slope = dmframe.loc[np.isclose(dmframe.index, gdata['displayedDM'])]['slope (mhz/ms)'].iloc[0]

	subname = None if 'resultidx' not in userdata else subburstdf.index[userdata['resultidx']]
	if ('resultidx' not in userdata) or (subname == burstname):
		gdata['displayedBurst'] = burstname
		wfall = gdata['wfall'].copy()
		wfall_cr = driftrate.cropwfall(wfall, twidth=dpg.get_value('twidth'))
		wfall_dd_cr = driftrate.dedisperse(wfall_cr, ddm, lowest_freq, df, dt)
	elif ('resultidx' in userdata) and (subname != burstname):
		gdata['displayedBurst'] = subname
		subbursts = getSubbursts()
		subburst = subbursts[subname]
		wfall_cr = subburst
		wfall_dd_cr = driftrate.dedisperse(wfall_cr, ddm, lowest_freq, df, dt)
		dmframe = subburstdf.loc[subname].set_index('DM')
		popt = dmframe.loc[np.isclose(dmframe.index, gdata['displayedDM'])].iloc[0][poptcols]
		slope = dmframe.loc[np.isclose(dmframe.index, gdata['displayedDM'])]['slope (mhz/ms)'].iloc[0]

	tseries = np.nanmean(wfall_dd_cr, axis=0)

	extents, correxts = driftrate.getExtents(wfall_dd_cr, df=df, dt=dt, lowest_freq=lowest_freq)
	gdata['extents'], gdata['correxts'] = extents, correxts
	dpg.set_value('twidth_ms', wfall_dd_cr.shape[1]*dt)
	dpg.set_value('pkidx', int(np.nanargmax(np.nanmean(wfall_dd_cr, axis=0))))

	corr = driftrate.autocorr2d(wfall_dd_cr)

	## enable scale sliders
	mostmin, mostmax = np.min(wfall_dd_cr), np.max(wfall_dd_cr)
	mmincorr, mmaxcorr = np.min(corr), np.max(corr)
	if sender not in ['wfallscale', 'corrscale']:
		dpg.configure_item('wfallscale', enabled=True, min_value=mostmin, max_value=mostmax,
							format='%.{}f'.format(getscale(mostmin, mostmax)+1))
		dpg.configure_item('corrscale', enabled=True, min_value=mmincorr, max_value=mmaxcorr,
							format='%.{}f'.format(getscale(mmincorr, mmaxcorr)+1))

	smin, smax = mostmin, mostmax
	scmin, scmax = mmincorr, mmaxcorr/20
	if sender == 'wfallscale':
		smin, smax = dpg.get_value('wfallscale')[:2]
	if sender == 'corrscale':
		scmin, scmax = dpg.get_value('corrscale')[:2]
	dpg.set_value('wfallscale', [smin, smax])
	dpg.set_value('corrscale', [scmin, scmax])

	flatwfall = list(np.flipud(wfall_dd_cr).flatten())
	if dpg.does_item_exist("Waterfall"):
		dpg.delete_item('Waterfall')
	dpg.add_heat_series(flatwfall,
		wfall_dd_cr.shape[0], wfall_dd_cr.shape[1],
		tag="Waterfall", parent="freq_axis", scale_min=smin, scale_max=smax,
		bounds_min=(extents[0],extents[2]), bounds_max=(extents[1], extents[3]), format='')
	for axis in ['time_axis', 'freq_axis']: dpg.fit_axis_data(axis)

	p0 = gdata['p0'] if dpg.get_value('EnableP0Box') else None
	corr2dtexture, txwidth, txheight = getcorr2dtexture(corr, popt, p0, correxts, slope, clim=scmax)
	if dpg.does_alias_exist('corr2dtexture'):
		dpg.delete_item('Corr2d') # texture won't delete unless items using it are deleted too
		time.sleep(0.017) # hack for strange race bug when deleting textures. wait one frame
		dpg.delete_item('corr2dtexture')
	dpg.add_static_texture(txwidth, txheight, corr2dtexture, tag='corr2dtexture', parent="TextureRegistry")
	dpg.add_image_series("corr2dtexture", parent="freq_lag_axis",
		bounds_min=[correxts[0],correxts[2]], bounds_max=[correxts[1], correxts[3]],
		tag='Corr2d')
	for axis in ['time_lag_axis', 'freq_lag_axis']: dpg.fit_axis_data(axis)

	tx = np.linspace(extents[0], extents[1], num=len(tseries))
	if dpg.does_alias_exist('TimeSeries'): dpg.delete_item('TimeSeries')
	with dpg.theme() as theme:
		with dpg.theme_component(dpg.mvLineSeries):
			dpg.add_theme_color(dpg.mvPlotCol_Line, (60+19, 103+10, 164+10, 255), category=dpg.mvThemeCat_Plots)
	dpg.add_line_series(tx, tseries, tag="TimeSeries", parent="tseries_int_axis")
	dpg.bind_item_theme('TimeSeries', theme)
	for axis in ['tseries_t_axis', 'tseries_int_axis']: dpg.fit_axis_data(axis)

def twidth_cb(_, data):
	twidth = data
	wfall_cr = getCurrentBurst()[4]
	if round(wfall_cr.shape[1]/2) != twidth:
		print(f"twidth_cb: {wfall_cr.shape = } {twidth = }")
	for regid in range(1, gdata['multiburst']['numregions']):
		if dpg.does_item_exist('RegionSelector{}'.format(regid)):
			maxval = wfall_cr.shape[1]*gdata['burstmeta']['tres']
			dpg.configure_item('Region{}'.format(regid), max_value=maxval, speed=maxval*0.005)
	plotdata_cb(_, None, None)

def pkidx_cb(_, data):
	enabled = dpg.get_value('pkidxbool')
	if dpg.get_value('pkidxbool'):
		pkidx = dpg.get_value('pkidx')
		print(pkidx)
		twidth_cb(_, data)

def subsample_cb(sender, data):
	if sender == 'ResetSamplingBtn':
		dpg.set_value('numfreqinput', gdata['wfall_original'].shape[0])
		dpg.set_value('numtimeinput', gdata['wfall_original'].shape[1])

	numf, numt = dpg.get_value("numfreqinput"), dpg.get_value("numtimeinput")

	try:
		# Make a copy of the original fall, apply the masks, then downsample
		wfall = applyMasks(np.copy(gdata['wfall_original']))
		if dpg.get_value('EnableSubBGBox'):
			tleft, tright, _, _ = dpg.get_value('SubtractBGRegion')
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
		updatedata_cb('subsample_cb', data, {'subsample': True})
		return subfall

def directory_cb(sender, data):
	path = data['file_path_name']
	dpg.set_value('Dirtext', 'Selected: {}'.format(path))
	dpg.configure_item('Filter', enabled=True)
	dpg.configure_item('clearfilter', enabled=True)
	files = sorted(glob.glob(path+'/{}'.format(gdata['globfilter'])))
	dpg.configure_item('burstselect', items=[os.path.basename(x) for x in files])
	gdata['datadir'] = path
	gdata['files']   = files
	log_cb(sender, path)

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
	filename = dpg.get_value('burstselect')
	filename = gdata['files'][[os.path.basename(f) for f in gdata['files']].index(filename)]
	updatedata_cb(sender, data, {'filename': filename})
	log_cb(sender, f'Opening file {filename}')

def exportmask_cb(sender, data):
	datestr = datetime.now().strftime('%b%d')
	filename = f'masks_{datestr}.npy'
	np.save(filename, [gdata['masks']])
	dpg.set_value('maskstatus', f'saved: {filename}')
	print(f'Saved {filename}')

def importmask_cb(sender, data):
	if data is None:
		return
	if type(data) == str:
		filename = data
	filename = data['file_path_name']
	log_cb(sender, 'mask selected: {}'.format(data))
	if filename.split('.')[-1] == 'npy':
		masks = np.load(filename, allow_pickle=True)[0]
		if type(masks) == dict:
			gdata['masks'].update(masks)
			updatedata_cb(sender, {'filename': gdata['currfile']}, None)
			masktable_cb(sender, None)
		else:
			error_log_cb(sender, 'invalid mask dictionary selected.')
	else:
		error_log_cb(sender, 'invalid mask file selected.')

def removemask_cb(sender, data, userdata):
	mask = userdata
	if mask in gdata['masks'][gdata['currfile']]['chans']:
		gdata['masks'][gdata['currfile']]['chans'].remove(mask)
		logwin.log_debug('removing {} from {} mask'.format(mask, gdata['currfile']))
		updatedata_cb(sender, {'keepview': True}, None)
		masktable_cb(sender, None)

def masktable_cb(sender, data):
	dpg.delete_item('Masktable')

	tableheight = 125
	with dpg.table(tag='Masktable', header_row=True, height=tableheight,
		borders_innerV=True, borders_outerH=True, borders_outerV=True,
		policy=dpg.mvTable_SizingFixedFit,
		scrollX=True, parent='Masking', before='MaskRangeGroup'):
		shortnames = [s.split('.')[0][-8:] for s in gdata['masks'].keys()]

		numcols = 0
		for key, masks in gdata['masks'].items():
			numcols = max(numcols, len(masks['chans']))
		for col in range(0, 1+numcols): # enough cols for the burst name and masks
			if col == 0: dpg.add_table_column(label="Burst")
			if col == 1: dpg.add_table_column(label="chan#:")
			if col > 1: dpg.add_table_column()

		for (key, masks), name in zip(gdata['masks'].items(), shortnames):
			with dpg.table_row():
				dpg.add_text(name)
				for chanmask in masks['chans']:
					dpg.add_selectable(label=chanmask, callback=removemask_cb, user_data=chanmask)


def sksgmask_cb(sender, data):
	items = ['SKSGSigmaInput','SKSGWindowInput']
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
	updatedata_cb(sender, {'keepview': True}, None)

def resulttable_cb(sender, data, userdata):
	displayresult_cb('User', data, {'trialDM': userdata})

def updateResultTable(resultsdf):
	if 'slope (mhz/ms)' in resultsdf.columns:
		resultsdf = driftlaw.computeModelDetails(resultsdf)
	dpg.delete_item('Resulttable')

	# subset of driftrate.columns:
	columns = ['name', 'DM', 'amplitude', 'slope (mhz/ms)', 'tau_w_ms', 'angle', 'center_f']
	# columns = ['name', 'DM', 'amplitude', 'tsamp_width','subbg_start (ms)', 'subbg_end (ms)']
	with dpg.table(tag='Resulttable', header_row=True, height=225, parent='ResultsGroup',
		borders_innerV=True, borders_outerH=True, borders_outerV=True,
		scrollY=True,
		resizable=True,
		policy=dpg.mvTable_SizingStretchSame,
		row_background=True,
		reorderable=True):
		for col in columns:
			dpg.add_table_column(label=col)

		# [burstname, trialDM, center_f, slope, slope_err, theta, red_chisq], popt, perr, [fres_MHz, tres_ms/1000]
		for burstname, rowdata in resultsdf.iterrows():
			trialDM = rowdata['DM']
			with dpg.table_row():
				# adding a component automatically moves to the next cell
				dpg.add_selectable(label=burstname, height=20, span_columns=True,
					callback=resulttable_cb, user_data=trialDM)
				for col in columns[1:]:
					dpg.add_text(f"{rowdata[col]}")

def mousemask_cb(sender, data):
	isOnWaterfall = dpg.is_item_hovered('WaterfallPlot')
	if isOnWaterfall:
		tchan, fchan = dpg.get_plot_mouse_pos()
		rawmask = round(fchan)

		# map frequency (rawmask) to channel number
		spectralrange = [round(gdata['extents'][2]), round(gdata['extents'][3])]
		mask = np.interp(rawmask, spectralrange, [0, gdata['wfall_original'].shape[0]])
		mask = int(mask)

		if mask not in gdata['masks'][gdata['currfile']]['chans']:
			gdata['masks'][gdata['currfile']]['chans'].append(mask)

		updatedata_cb(sender, {'keepview': True}, None)
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
	# pkidx = dpg.get_value('pkidx') if dpg.get_value('pkidxbool') else None
	wfall_cr = driftrate.cropwfall(wfall, twidth=dpg.get_value('twidth'), pkidx=None)
	burstDM = gdata['burstDM']
	return burstname, df, dt, lowest_freq, wfall_cr, burstDM

def getMeasurementInfo(wfall_cr):
	tsamp_width = dpg.get_value('twidth')
	downf = gdata['wfall_original'].shape[0] / gdata['wfall'].shape[0]
	downt = gdata['wfall_original'].shape[1] / gdata['wfall'].shape[1]
	fchans, tchans = wfall_cr.shape
	subbgstart, subbgend, sksigma, skwindow = None, None, None, None
	if dpg.get_value('EnableSubBGBox'):
		subbgstart, subbgend, _, _= dpg.get_value('SubtractBGRegion')
	if dpg.get_value('EnableSKSGMaskBox'):
		sksigma, skwindow = dpg.get_value('SKSGSigmaInput'), dpg.get_value('SKSGWindowInput')
	cols = ['downf', 'downt', 'fchans', 'tchans', 'tsamp_width','subbg_start (ms)', 'subbg_end (ms)','sksigma','skwindow']
	row = [downf, downt, fchans, tchans, tsamp_width, subbgstart, subbgend, sksigma, skwindow]

	# TODO: region info, raw shape
	regions = getAllRegions() # is empty when regions is disabled
	for regname, region in regions.items():
		if regname == 'background':
			cols.append('background')
			row.append(region[1])
		else:
			cols.append(f'regstart_{regname}')
			row.append(region[0])
			cols.append(f'regend_{regname}')
			row.append(region[1])

	return cols, row

def progress_cb(val, data):
	dpg.set_value('SlopeStatus', val)
	overlay = dpg.get_item_configuration("SlopeStatus")['overlay'].split(' (')[0]
	dpg.configure_item('SlopeStatus', overlay=f"{overlay} ({data})")

def slope_cb(sender, data, userdata):
	burstname, df, dt, lowest_freq, wfall_cr, burstDM = getCurrentBurst()

	if userdata is not None:
		dpg.configure_item('SlopeStatus', overlay='Status: Doing second pass...')
		p0 = userdata['p0']
		trialDMs = userdata['badfitDMs']
	else:
		dpg.configure_item('SlopeStatus', overlay='Status: Calculating...')
		p0 = [] if not gdata['p0'] else gdata['p0']
		trialDMs = np.unique(np.round(np.append(gdata['trialDMs'], burstDM), decimals=5))
		if not dpg.get_value('RepeatBox') and gdata['resultsdf'] is not None:
			trialDMs = list(set(gdata['resultsdf'].loc[burstname]['DM']) ^ set(trialDMs))
			# remove DMs that are in gdata['resultsdf'] from trialDMs

	progress = io.StringIO()
	results, burstdf = driftrate.processDMRange(burstname, wfall_cr, burstDM, trialDMs, df, dt,
												lowest_freq, p0,
												tqdmout=None, progress_cb=progress_cb)

	if gdata['multiburst']['enabled']:
		subbursts, corrsigma, wfallsigma = getSubbursts(getsigmas=True)
		subresults, subdf = [], pd.DataFrame()
		for subname, subburst in subbursts.items():
			print('processing {}'.format(subname))
			ret, retdf = driftrate.processDMRange(subname, subburst, burstDM, trialDMs, df, dt,
												  lowest_freq, corrsigma=corrsigma, wfallsigma=wfallsigma)
			subresults.append(ret)
			subdf = subdf.append(retdf)

		burstdf = burstdf.append(subdf)

	# Do a second pass using the best p0 just found
	p0 = getOptimalFit(burstdf)
	if userdata is None:
		return slope_cb(sender, data, {'p0' : p0, 'badfitDMs': trialDMs})

	# Add measurement info to row (things needed to reproduce/reload the measurement)
	cols, row = getMeasurementInfo(wfall_cr)
	print('measurement info >>', cols, row)
	burstdf[cols] = row

	# Save results
	if gdata['resultsdf'] is None:
		gdata['resultsdf'] = burstdf
	else:
		# overwrite if there are already results while adding new results
		gdata['resultsdf'] = (burstdf.append(gdata['resultsdf']) # add new and updated measurements
									.reset_index().drop_duplicates(['name', 'DM']) # drop old measurements
									.set_index('name').sort_values(['name', 'DM'])) # sort table
	backupresults()
	print(gdata['resultsdf'])

	burstdf = gdata['resultsdf'][gdata['resultsdf'].index.str.startswith(burstname)]
	gdata['burstdf'] = burstdf

	dpg.configure_item('SlopeStatus', overlay='Status: Done.')
	dpg.set_value('NumMeasurementsText', "# of Measurements for this burst: {}".format(len(burstdf)))
	dpg.set_value('TotalNumMeasurementsText', "Total # of Measurements: {}".format(len(gdata['resultsdf'])))
	dpg.configure_item('PrevDM', enabled=True)
	dpg.configure_item('NextDM', enabled=True)
	dpg.configure_item('ExportCSVBtn', enabled=True)
	dpg.configure_item('ExportPDFBtn', enabled=True)
	dpg.configure_item('EnableP0Box', enabled=True)
	initializeP0Group()
	updateResultTable(burstdf)
	plotdata_cb(sender, data, userdata)

def redodm_cb(sender, data, userdata):
	p0 = gdata['p0']
	useForAll = dpg.get_value('P0AllDMsBox')

	burstname, df, dt, lowest_freq, wfall_cr, burstDM = getCurrentBurst()
	displayedname = gdata['displayedBurst']
	corrsigma, wfallsigma = None, None
	if gdata['multiburst']['enabled'] and displayedname != burstname:
		subbursts, corrsigma, wfallsigma = getSubbursts(getsigmas=True)
		wfall_cr = subbursts[displayedname]

	print('redoing ', displayedname, gdata['displayedDM'], f'with {df = } {dt =} {lowest_freq = }')
	result, burstdf = driftrate.processDMRange(
		displayedname, wfall_cr, burstDM, [float(gdata['displayedDM'])],
		df, dt, lowest_freq, p0=p0, corrsigma=corrsigma, wfallsigma=wfallsigma
	)
	print(f'{result = }')
	cols, row = getMeasurementInfo(wfall_cr)
	burstdf[cols] = row
	df = gdata['resultsdf']
	df[(df.index == displayedname) & (np.isclose(df['DM'], gdata['displayedDM']))] = burstdf

	gdata['resultsdf'] = df
	gdata['burstdf'] = gdata['resultsdf'].loc[burstname]
	dispdm = gdata['displayedDM']

	if gdata['multiburst']['enabled']:
		subburstdf = gdata['resultsdf'][gdata['resultsdf'].index.str.startswith(burstname)]
		updateResultTable(subburstdf)
		subburstdf = subburstdf.reset_index()
		resultidx = subburstdf[(subburstdf.name == displayedname) & (np.isclose(subburstdf.DM, dispdm))].index[0]
		if userdata is None:
			userdata = {}
		userdata['resultidx'] = resultidx
	else:
		updateResultTable(gdata['burstdf'])
	backupresults()
	plotdata_cb(sender, data, userdata)

	# if useForAll programatically click "next DM" then repeat this block
	# that will preserve the behaviour where regions only need to be read at the start
	if useForAll:
		tabledf = gdata['resultsdf'][gdata['resultsdf'].index.str.startswith(burstname)].reset_index()
		redoidx = tabledf[(tabledf.name == gdata['displayedBurst']) & (np.isclose(tabledf.DM, dispdm))].index[0]
		maxidx = len(tabledf)
		if redoidx+1 < maxidx:
			displayresult_cb('NextDM', None, None)
			redodm_cb(sender, None, None)

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
	if data is None:
		return
	if type(data) == str:
		resultfile = data
	resultsfile = data['file_path_name']
	resultsdf = pd.read_csv(resultsfile).set_index('name')
	gdata['resultsdf'] = resultsdf
	regionsobj = driftrate.readRegions(resultsdf)
	gdata['multiburst']['regions'] = regionsobj
	burstselect_cb(sender, data)

def exportresults_cb(sender, data):
	global exportPDF
	resultsdf = gdata['resultsdf']
	df = driftlaw.computeModelDetails(resultsdf)

	datestr = datetime.now().strftime('%b%d')
	prefix = dpg.get_value('ExportPrefix')
	filename = f'{prefix}_{len(df.index)}rows_{datestr}.csv'

	if sender == 'ExportCSVBtn' or sender == 'ExportPDFBtn':
		dpg.configure_item('ExportCSVText', show=True)
		try:
			df.to_csv(filename)
			dpg.set_value('ExportCSVText', 'Saved to {}'.format(filename))
		except PermissionError as e:
			dpg.set_value('ExportCSVText', 'Permission Denied')
		if sender == "ExportPDFBtn":
			exportPDF = True # callback thread lets main thread know to run matplotlib on next frame
	elif sender == 'MainThread':
		success = driftrate.plotResults(filename, datafiles=gdata['files'], masks=gdata['masks'])
		if success:
			dpg.set_value('ExportPDFText', f'Saved to {filename.split(".")[0]+".pdf"}')
		else:
			dpg.set_value('ExportPDFText', f'Permission Denied')

def backupresults():
	resultsdf = gdata['resultsdf']
	df = driftlaw.computeModelDetails(resultsdf)
	datestr = datetime.now().strftime('%b%d')
	prefix = 'backup'
	if not os.path.isdir('backups'):
		os.makedirs('backups', exist_ok=True)
	filename = 'backups/{}_results_{}.csv'.format(prefix, datestr)
	df.to_csv(filename)

def displayresult_cb(sender, data, userdata):
	burstDM = gdata['burstDM']
	dmlist = list(gdata['burstdf'].loc[gdata['burstdf'].index[0]]['DM'])

	burstname = dpg.get_value('burstname').replace(',', '')
	subname = gdata['displayedBurst']
	dispdm = gdata['displayedDM']
	subburstdf = gdata['resultsdf'][gdata['resultsdf'].index.str.startswith(burstname)].reset_index()
	resultidx = subburstdf[(subburstdf.name == subname) & (np.isclose(subburstdf.DM, dispdm))].index[0]
	subburstdf = subburstdf.set_index('name')

	if sender == 'User':
		resultidx = dmlist.index(userdata['trialDM'])
	elif sender == 'NextDM':
		resultidx = resultidx + 1
		if not (resultidx < len(subburstdf)):
			resultidx = 0
	elif sender == 'PrevDM':
		resultidx = resultidx - 1

	if userdata is None:
		userdata = {}
	userdata['resultidx'] = resultidx
	subname = subburstdf.index[resultidx]
	gdata['displayedDM'] = dmlist[resultidx % len(dmlist)]
	dpg.set_value('dmdisplayed', str(round(gdata['displayedDM'], getscale(gdata['displayedDM']))))
	dpg.set_value('burstdisplayed', subname)
	plotdata_cb(sender, data, userdata)

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
	p0 = [float(bestrow[param]) for param in params]
	return p0

def initializeP0Group():
	p0 = []
	burstname = dpg.get_value('burstname').replace(',', '')
	subburstdf = gdata['resultsdf'][gdata['resultsdf'].index.str.startswith(burstname)]
	wfall_cr = getCurrentBurst()[-2]
	p0 = getOptimalFit(subburstdf)
	p0f = [p0i for p0i in p0] # why?
	if p0f[0] < 0:
		p0f = [-p0fi if p0fi < 0 else p0fi for p0fi in p0f]

	gdata['p0'] = p0
	dpg.set_value("AmplitudeDrag", p0f[0])
	dpg.set_value("AngleDrag", p0f[5])
	# solution will always be near the center
	dpg.set_value("x0y0Drag", [0.01, 0.01]) # dpg glitches if you use 0,0
	dpg.set_value("SigmaXYDrag", [p0f[3], p0f[4]])
	for item in ['AmplitudeDrag', 'AngleDrag', 'x0y0Drag', 'SigmaXYDrag']:
		val = dpg.get_value(item)
		if type(val) == list:
			val = val[0]
		dpg.configure_item(item, speed=1/10**getscale(val), format='%.{}f'.format(getscale(val)+1))
	dpg.configure_item('AngleDrag', speed=0.01, format='%.8f')
	return p0

def enablep0_cb(sender, data, userdata):
	toggle_config(sender, userdata)
	updatep0_cb(sender, data, userdata)

def updatep0_cb(sender, data, userdata):
	if not userdata: userdata = {}
	p0 = []

	# get resultidx
	subname, dispdm = gdata['displayedBurst'], gdata['displayedDM']
	burstname = dpg.get_value('burstname').replace(',', '')
	subburstdf = gdata['resultsdf'][gdata['resultsdf'].index.str.startswith(burstname)].reset_index()
	resultidx = subburstdf[(subburstdf.name == subname) & (np.isclose(subburstdf.DM, dispdm))].index[0]
	if gdata['multiburst']['enabled']:
		userdata['resultidx'] = resultidx

	for item in ['AmplitudeDrag', 'x0y0Drag', 'SigmaXYDrag', 'AngleDrag']:
		val = dpg.get_value(item)
		if type(val) == list: # x0y0Drag or SigmaXYDrag
			p0 = p0 + val[:2]
		else:
			p0.append(val)
	gdata['p0'] = p0
	print(f"{sender = }, {data = }, {userdata = }")
	plotdata_cb(sender, data, userdata)

def enablesplitting_cb(sender, data):
	items = ['Region', 'RegionType', 'RemoveRegionBtn', 'AddRegionBtn']
	gdata['multiburst']['enabled'] = dpg.get_value(sender)
	items = items.copy()
	items.remove('AddRegionBtn')
	toggle_config(sender, {'kwargs': ['enabled'], 'items': ['AddRegionBtn']})
	for regid in range(1, gdata['multiburst']['numregions']):
		if dpg.does_item_exist('RegionSelector{}'.format(regid)):
			itemsid = list(map(lambda item: item+'{}'.format(regid), items))
			toggle_config(sender, {'kwargs': ['enabled'], 'items': itemsid})

			if gdata['multiburst']['enabled']:
				drawregion_cb([regid], None)
			else:
				seriesnames = [f'Region{regid}Series', f'Region{regid}TSeries']
				for axis, seriesname in zip(['freq_axis', 'tseries_int_axis'], seriesnames):
					dpg.delete_item(seriesname)

def addregion_cb(sender, _):
	regionSelector()
	drawregion_cb([gdata['multiburst']['numregions']-1], None)

def removeregion_cb(sender, _):
	regid = sender[-1]
	dpg.delete_item('RegionSelector{}'.format(regid))
	seriesname = 'Region{}Series'.format(regid)
	seriesnames = [f'Region{regid}Series', f'Region{regid}TSeries']
	for seriesname in seriesnames:
		dpg.delete_item(seriesname)
	# dpg.delete_series(seriesname) # needs to be removed from both plots, tags are unique
	# [dpg.delete_series(plot, seriesname) for plot in ['WaterfallPlot', 'TimeSeriesPlot']]
	gdata['multiburst']['numregions'] -= 1

colors = [
	(0, 255, 0),
	(255, 0, 0),
	(0, 0, 255),
	(255, 255, 0),
	(255, 0, 255),
	(0, 255, 255)
]
def drawregion_cb(sender, _):
	regid = sender[-1]
	region = dpg.get_value('Region{}'.format(regid))
	regiontype = dpg.get_value('RegionType{}'.format(regid))

	if region[1] < region[0]:
		region.sort()
		dpg.set_value('Region{}'.format(regid), region)

	seriesnames = [f'Region{regid}Series', f'Region{regid}TSeries']
	for axis, seriesname in zip(['freq_axis', 'tseries_int_axis'], seriesnames):
		dpg.delete_item(seriesname)
		with dpg.theme() as linetheme:
			with dpg.theme_component(dpg.mvVLineSeries):
				dpg.add_theme_color(dpg.mvPlotCol_Line, colors[(int(regid)-1) % len(colors)],
					category=dpg.mvThemeCat_Plots)
				dpg.add_theme_style(dpg.mvPlotStyleVar_LineWeight, 1.5,
					category=dpg.mvThemeCat_Plots)
		dpg.add_vline_series(region, tag=seriesname, parent=axis)
		dpg.bind_item_theme(seriesname, linetheme)

def addMaskRange_cb(sender, rangeid):
	maxchan = gdata['wfall_original'].shape[0]-1
	if not rangeid:
		gdata["maskrangeid"] += 1
		rangeid = gdata["maskrangeid"]
	with dpg.group(tag=f'MaskRange{rangeid}', horizontal=True, parent='MaskRangeGroup'):
		dpg.add_drag_intx(tag=f'rangeslider##{rangeid}', label='Mask Range', width=200,
			callback=maskrange_cb, size=2,
			min_value=0, max_value=maxchan)
		dpg.add_button(tag=f'removerange##{rangeid}', label='X', enabled=True,
			callback=removemaskrange_cb)
	if rangeid not in gdata['masks'][gdata['currfile']]['ranges']:
		gdata['masks'][gdata['currfile']]['ranges'][rangeid] = [0, 0]

def removemaskrange_cb(sender, data):
	rangeid = int(sender[-1])
	if rangeid in gdata['masks'][gdata['currfile']]['ranges']:
		del gdata['masks'][gdata['currfile']]['ranges'][rangeid]
	dpg.delete_item(f'MaskRange{rangeid}')
	updatedata_cb(sender, {'keepview': True}, None)

def maskrange_cb(sender, data):
	rangeid = int(sender[-1])
	*maskrange, _, _ = dpg.get_value(sender)
	# print(gdata['masks'][gdata['currfile']], '\n----')
	gdata['masks'][gdata['currfile']]['ranges'][rangeid] = maskrange
	updatedata_cb(sender, {'keepview': True}, None)
	# masktable_cb(sender, None)

def regionSelector():
	regid = gdata['multiburst']['numregions']
	enabled = gdata['multiburst']['enabled']
	if 'wfall' in gdata:
		maxval =  gdata['extents'][1]
	else:
		maxval = 100

	with dpg.group(tag='RegionSelector{}'.format(regid), horizontal=True, parent='SplittingSection',
		before="AddRegionBtn"):
		dpg.add_drag_floatx(tag='Region{}'.format(regid),
			size=2,
			label='(ms)',
			width=280,
			enabled=enabled,
			max_value=maxval,
			speed=maxval*0.005,
			default_value=[0, maxval],
			callback=drawregion_cb
		)
		dpg.add_text(label="(?)", color=[150, 150, 150])
		with dpg.tooltip(dpg.last_item()):
			dpg.add_text("double click to edit")
		dpg.add_radio_button(tag='RegionType{}'.format(regid), items=["Background", "Burst"],
			horizontal=True,
			callback=drawregion_cb,
			enabled=enabled,
			default_value="Background"
		)
		dpg.add_button(tag='RemoveRegionBtn{}'.format(regid), label='X', enabled=enabled, callback=removeregion_cb)
	gdata['multiburst']['numregions'] += 1

subburst_suffixes = driftrate.subburst_suffixes
def getSubbursts(getsigmas=False):
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
			if regiontype == "Background":
				background = wfall_cr[:, region[0]:region[1]]
			elif regiontype == "Burst":
				subburst = wfall_cr[:, region[0]:region[1]]
				subbursts.append(subburst)

	if background is None:
		error_log_cb("getSubbursts", "Please specify a background region")
		raise ValueError("Please specify a background region")

	subburstsobj = {}
	# pad with background
	for subburst, suffix in zip(subbursts, subburst_suffixes):
		subburst = np.concatenate((0*background, subburst, 0*background), axis=1)
		subname = burstname +'_'+ suffix # '_' is used in plotResults to split names
		subburstsobj[subname] = subburst

	# zero padding bursts ruins the sampling of sigma so we will sample from the full waterfall
	if getsigmas:
		corr = driftrate.autocorr2d(wfall_cr)
		corrsigma = np.std( corr[:, 0:50] )
		wfallsigma = np.std( wfall_cr[:, 0:wfall_cr.shape[1]//20] ) # use the first 5% of channels
		return subburstsobj, corrsigma, wfallsigma
	else:
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

def dmchange_cb(sender, data):
	newDM = round(data, 3)

	# update
	gdata['displayedDM'] = newDM
	dpg.set_value('dmdisplayed', str(gdata['displayedDM']))
	if gdata['wfall'].shape == gdata['wfall_original'].shape:
		dpg.configure_item('SaveDMButton',
			enabled=bool(gdata['burstmeta']['DM'] != gdata['displayedDM']))

	plotdata_cb(sender, None, None)

def savedm_cb(sender, data):
	newDM = round(dpg.get_value('DM'), 3)
	npz = gdata['currfile']
	if gdata['wfall'].shape == gdata['wfall_original'].shape:
		if os.path.isfile(npz):
			# write new DM and waterfall to disk, then trigger a reload
			df, dt      = gdata['burstmeta']['fres'], gdata['burstmeta']['tres']
			lowest_freq = min(gdata['burstmeta']['dfs']) # mhz
			ddm = newDM - gdata['burstDM']
			wfall = np.copy(gdata['wfall_original'])
			wfall = driftrate.dedisperse(wfall, ddm, lowest_freq, df, dt)

			driftrate.updatenpz(npz, 'wfall', wfall)
			driftrate.updatenpz(npz, 'DM', newDM)
			burstselect_cb(sender, None)
			print(f'updated DM of {npz} to {newDM}')
			dpg.configure_item('SaveDMButton', enabled=False)

def frbgui(filefilter=gdata['globfilter'],
		datadir='',
		maskfile=None,
		regionfile=None,
		dmrange=None,
		numtrials=10,
		dmstep=0.1,
		winwidth=1700,
		winheight=850,
	):
	global logwin, exportPDF
	dpg.create_context()

	#### Themeing
	darktheme = themes.create_theme_imgui_dark()
	lighttheme = themes.create_theme_imgui_light()

	with dpg.theme() as global_theme:
		with dpg.theme_component(dpg.mvAll):
			dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 4, category=dpg.mvThemeCat_Core)
			dpg.add_theme_style(dpg.mvStyleVar_WindowRounding, 4)

		# fix for disabled theme see DearPyGui/issues/2068
		disabledcolor = (45, 45, 48)
		comps = [dpg.mvInputText, dpg.mvButton, dpg.mvRadioButton, dpg.mvTabBar, dpg.mvTab, dpg.mvImage, dpg.mvMenuBar, dpg.mvViewportMenuBar, dpg.mvMenu, dpg.mvMenuItem, dpg.mvChildWindow, dpg.mvGroup, dpg.mvDragFloatMulti, dpg.mvSliderFloat, dpg.mvSliderInt, dpg.mvFilterSet, dpg.mvDragFloat, dpg.mvDragInt, dpg.mvInputFloat, dpg.mvInputInt, dpg.mvColorEdit, dpg.mvClipper, dpg.mvColorPicker, dpg.mvTooltip, dpg.mvCollapsingHeader, dpg.mvSeparator, dpg.mvCheckbox, dpg.mvListbox, dpg.mvText, dpg.mvCombo, dpg.mvPlot, dpg.mvSimplePlot, dpg.mvDrawlist, dpg.mvWindowAppItem, dpg.mvSelectable, dpg.mvTreeNode, dpg.mvProgressBar, dpg.mvSpacer, dpg.mvImageButton, dpg.mvTimePicker, dpg.mvDatePicker, dpg.mvColorButton, dpg.mvFileDialog, dpg.mvTabButton, dpg.mvDrawNode, dpg.mvNodeEditor, dpg.mvNode, dpg.mvNodeAttribute, dpg.mvTable, dpg.mvTableColumn, dpg.mvTableRow]
		themecolors = [dpg.mvThemeCol_Button, dpg.mvThemeCol_ButtonActive, dpg.mvThemeCol_ButtonHovered, dpg.mvThemeCol_FrameBg, dpg.mvThemeCol_FrameBgActive, dpg.mvThemeCol_FrameBgHovered]
		for comp_type in comps:
			with dpg.theme_component(comp_type, enabled_state=False):
				dpg.add_theme_color(dpg.mvThemeCol_Text, (0.50 * 255, 0.50 * 255, 0.50 * 255, 1.00 * 255))
				for color in themecolors:
					dpg.add_theme_color(color, disabledcolor)

		dpg.bind_theme(global_theme)
	####

	dpg.create_viewport(title='frbgui', width=winwidth, height=winheight,
						# x_pos=1600, y_pos=-550, # with two monitors
						x_pos=0, y_pos=0, # one monitor
						small_icon='frbgui.ico', large_icon='frbgui.ico')
	dpg.add_texture_registry(tag="TextureRegistry")
	logwin = logger.mvLogger() # todo: replace with something more usable

	gdata['datadir'] = datadir
	with dpg.window(label='FRB Analysis', width=560, height=745, pos=[10, 30]):
		with dpg.collapsing_header(label="1. Data", default_open=True):
			dpg.add_file_dialog(show=False, directory_selector=True,
				callback=directory_cb,
				tag='dirselector', label="Select Directory", modal=True, width=700, height=400)
			# dpg.add_file_extension('.npz', color=(255, 255, 255, 255//2), parent='dirselector')
			dpg.add_button(label="Select Directory...", callback=lambda s, d: dpg.show_item('dirselector'))
			dpg.add_text(tag="Dirtext", default_value="Selected: (no directory selected)")
			with dpg.group(horizontal=True):
				dpg.add_text("Filter:")
				dpg.add_input_text(tag="Filter", label='', hint="eg. *.npy", callback=filter_cb, enabled=False)
				dpg.add_button(tag='clearfilter', label='Clear filter', callback=clearfilter_cb, enabled=False)

			dpg.add_text("Files found:")
			with dpg.group(horizontal=True):
				dpg.add_listbox(items=[], tag="burstselect", label='', num_items=10, width=520, callback=burstselect_cb)
				with dpg.tooltip(dpg.last_item()):
					dpg.add_text("Select to load burst...")

			dpg.add_text("Burst Metadata:")
			dpg.add_input_text(tag='burstname', label='Burst Name')
			with dpg.group(horizontal=True):
				dpg.add_input_float(tag='DM', label='DM (pc/cm3)', callback=dmchange_cb, step=0.100000)
				dpg.add_button(tag="SaveDMButton", label='Save', callback=savedm_cb, enabled=False)
			dpg.add_input_float(tag='dt', label='Time Resolution (ms)', enabled=False)
			dpg.add_input_float(tag='df', label='Freq Resolution (MHz)', enabled=False)
			dpg.add_input_float(tag='center_f', label='Center Frequency (MHz)', enabled=False)
			dpg.add_input_float(tag='bandwidth', label='Bandwidth (MHz)', enabled=False)
			dpg.add_input_float(tag='duration', label='Data Duration', enabled=False)
			dpg.add_input_float(tag='burstSN', label='Burst SNR', enabled=False)
			dpg.add_input_text(tag='telescope', label='Telescope', enabled=False)
			## TODO: Use these units to populate the resolution inputs
			dpg.add_input_text(tag='freq_unit', label='Frequency Unit', enabled=False)
			dpg.add_input_text(tag='time_unit', label='Time Unit', enabled=False)
			dpg.add_input_text(tag='int_unit', label='Intensity Unit', enabled=False)
			dpg.add_input_int(tag='twidth', label='Display width (# chans)',
				default_value=twidth_default,
				step=10,
				callback=twidth_cb
			)
			dpg.add_input_float(tag='twidth_ms', label='Display width (ms)', enabled=False)
			with dpg.group(horizontal=True):
				dpg.add_input_int(tag='pkidx', label='Peak Index', enabled=False, callback=pkidx_cb)
				dpg.add_checkbox(tag='pkidxbool', label='Use ?', enabled=False, default_value=False,
					callback=lambda s, d: dpg.configure_item('pkidx', enabled=dpg.get_value('pkidxbool')))

			dpg.add_text('Dedispersion Range for all Bursts: ')
			dpg.add_text(default_value='Warning: Range chosen does not include burst DM', tag='DMWarning',
				color=[255, 0, 0], show=False)
			with dpg.group(horizontal=True):
				dpg.add_drag_floatx(tag='dmrange', label='DM range (pc/cm^3)', callback=dmrange_cb,
					size=2, min_value=0, max_value=0)
				dpg.add_text(label="(?)", color=[150, 150, 150])
				with dpg.tooltip(dpg.last_item()):
					dpg.add_text("Double click to edit")
			dpg.add_input_int(tag='numtrials', label='# of Trial DMs', default_value=10, callback=dmrange_cb)
			dpg.add_text(' or ')
			dpg.add_input_float(tag='dmstep', label='DM Step (pc/cm^3)', default_value=0.1,
				min_value=0.001,
				min_clamped=True,
				callback=dmrange_cb)

		with dpg.collapsing_header(label="2. Waterfall Cleanup", default_open=True):
			with dpg.tree_node(label='Subtract Background', default_open=True):
				dpg.add_checkbox(tag='EnableSubBGBox', label='Subtract background sample', callback=subtractbg_cb)
				dpg.add_drag_floatx(tag='SubtractBGRegion', size=2,
					label='t_start (ms), t_end (ms)',
					width=280,
					enabled=False,
					max_value=100,
					speed=0.5,
					default_value=[0, 25],
					callback=subtractbg_cb
				)

			with dpg.tree_node(label='Masking', tag='Masking', default_open=True):
				dpg.add_text("Click on the waterfall plot to begin masking frequency channels.")
				dpg.add_text("NOTE: only mask on the original waterfall (todo: add a 'mask' button)")

				with dpg.group(horizontal=True):
					dpg.add_file_dialog(tag='maskdialog', show=False, callback=importmask_cb,
						width=700, height=400, modal=True)
					dpg.add_file_extension('.npy', parent='maskdialog')
					dpg.add_button(label='Export Masks', callback=exportmask_cb, enabled=True)
					dpg.add_button(label='Import Masks', callback=lambda s, d: dpg.show_item('maskdialog'))
					dpg.add_text(tag='maskstatus', label='', color=(255, 255, 0))

				dpg.add_table(tag='Masktable', height=100)

				with dpg.group(tag="MaskRangeGroup"): # added to by the user later
					pass

				dpg.add_button(tag="MaskRangeButton", label='Add Mask Range', callback=addMaskRange_cb,
								enabled=True)

			with dpg.tree_node(label='Auto Masking', default_open=True):
				dpg.add_checkbox(tag='EnableSKSGMaskBox', label='Enable SK-SG Filter', default_value=False,
					callback=sksgmask_cb)
				with dpg.group(horizontal=True):
					dpg.add_input_int(tag="SKSGSigmaInput", width=100, default_value=3, label="sigma",
						callback=sksgmask_cb, enabled=False, min_value=1)
					dpg.add_input_int(tag="SKSGWindowInput", width=100, default_value=15, label="window",
						callback=sksgmask_cb, enabled=False, min_value=1)
				# dpg.add_text(tag='SKSGStatus', default_value=' ')

			with dpg.tree_node(label='Downsampling', default_open=True):
				dpg.add_text(tag="Wfallshapelbl", default_value="Original Size: (no burst selected)")
				dpg.add_text(tag="Subfallshapelbl", default_value="Current Size: (no burst selected)")
				with dpg.group(horizontal=True):
					dpg.add_input_int(tag="numfreqinput", width=100, label="numf", callback=subsample_cb, enabled=False)
					dpg.add_input_int(tag="numtimeinput", width=100, label="numt", callback=subsample_cb, enabled=False)
				dpg.add_button(tag='ResetSamplingBtn', label='Reset', callback=subsample_cb, enabled=False)

		with dpg.collapsing_header(tag="SplittingSection", label="3. Burst Splitting", default_open=True):
			dpg.add_checkbox(tag='MultiBurstBox', label='Are there multiple bursts in this waterfall?',
				default_value=False, enabled=True,
				callback=enablesplitting_cb,
				user_data={'kwargs': ['enabled']}
			)
			regionSelector()
			dpg.add_button(tag='AddRegionBtn', label="Add Region", callback=addregion_cb, enabled=False)

		with dpg.collapsing_header(tag='SlopeSection', label="4. Spectro-Temporal Measurements", default_open=True):
			with dpg.group(horizontal=True):
				dpg.add_button(label="Measure over DM Range", callback=slope_cb)
				dpg.add_progress_bar(tag="SlopeStatus", default_value=0,
					overlay="Status: Click to calculate",
					width=300
				)
			dpg.add_checkbox(tag='RepeatBox', label='Repeat measurements',
				default_value=True, enabled=True
			)
			with dpg.group(horizontal=True):
				dpg.add_file_dialog(tag='resultsdialog', show=False, callback=loadresults_cb,
						width=700, height=400, modal=True)
				dpg.add_file_extension('.csv', parent='resultsdialog')
				dpg.add_button(label="Load Results", callback=lambda s, d:dpg.show_item('resultsdialog'))
				dpg.add_text(tag="LoadResultsWarning", color=(255, 0, 0),
					default_value="Warning: loading will overwrite unsaved results")
			with dpg.group(horizontal=True):
				dpg.add_button(tag='DeleteBurstBtn',
					label="Delete this burst's results",
					callback=deleteresults_cb,
					user_data='burst',
					enabled=False
				)
				dpg.add_button(tag='DeleteAllBtn',
					label="Delete ALL results",
					callback=deleteresults_cb,
					user_data='all',
					enabled=False
				)

			dpg.add_text(tag='NumMeasurementsText', default_value="# of Measurements for this burst: (none)")
			dpg.add_text(tag='TotalNumMeasurementsText', default_value="Total # of Measurements: (none)")
			with dpg.group(tag="DMselector", horizontal=True):
				dpg.add_button(tag="PrevDM", arrow=True, direction=dpg.mvDir_Left, enabled=False,
					callback=displayresult_cb)
				dpg.add_button(tag="NextDM", arrow=True, direction=dpg.mvDir_Right, enabled=False,
					callback=displayresult_cb)
				dpg.add_text("Burst Displayed: ")
				dpg.add_text(tag="burstdisplayed", default_value=str(dpg.get_value('burstname')))
				dpg.add_text("DM Displayed (pc/cm^3): ")
				dpg.add_text(tag="dmdisplayed", default_value=str(dpg.get_value('DM')))
				dpg.add_text(label="(?)", color=[150, 150, 150])
				with dpg.tooltip(dpg.last_item()):
					dpg.add_text('Selecting a row also displays that DM')

			with dpg.group(tag='P0Group'):
				with dpg.tree_node(label='Fit Initial Guess', default_open=False):
					items = ['P0AllDMsBox','AmplitudeDrag','AngleDrag','x0y0Drag','SigmaXYDrag', 'RedoBtn']
					with dpg.group(horizontal=True):
						dpg.add_checkbox(tag="EnableP0Box", label='Use Initial Guess', default_value=False,
							enabled=False,
							callback=enablep0_cb,
							user_data={'kwargs': ['enabled'], 'items': items}
						)
						dpg.add_checkbox(tag='P0AllDMsBox', label='Use for all DMs', default_value=False, enabled=False)
					dpg.add_drag_float(tag="AmplitudeDrag", label="Amplitude", enabled=False,
						min_value=0, max_value=0, callback=updatep0_cb)
					dpg.add_drag_float(tag="AngleDrag", label="Angle", enabled=False,
						min_value=0, max_value=0, speed=0.01, callback=updatep0_cb)
					dpg.add_drag_floatx(tag="x0y0Drag", size=2, label="x0, y0", enabled=False,
						min_value=0, max_value=0, speed=1, callback=updatep0_cb)
					dpg.add_drag_floatx(tag="SigmaXYDrag", size=2, label="sigmax, sigmay", enabled=False,
						min_value=0, max_value=0, callback=updatep0_cb)
					dpg.add_button(tag="RedoBtn", label="Redo this DM", callback=redodm_cb, enabled=False)

			with dpg.group(tag='ResultsGroup'):
				dpg.add_table(tag='Resulttable', header_row=True, height=10, callback=resulttable_cb)

			dpg.add_input_text(tag='ExportPrefix', label='Filename Prefix', default_value="Test")
			with dpg.group(horizontal=True):
				dpg.add_button(tag='ExportCSVBtn', label="Export Results CSV", callback=exportresults_cb, enabled=False)
				dpg.add_text(tag='ExportCSVText', default_value=' ', show=True, color=(0, 255, 0))
			with dpg.group(horizontal=True):
				dpg.add_button(tag='ExportPDFBtn', label="Export Results PDF", callback=exportresults_cb, enabled=False)
				dpg.add_text(tag='ExportPDFText', default_value=' ', show=True, color=(0, 255, 0))

	### Plotting window
	with dpg.window(label="FRB Plots", width=1035, height=745, pos=[600,30]):
		with dpg.group(horizontal=True):
			dpg.add_slider_floatx(tag="wfallscale", size=2, label='Wfall Min/Max', enabled=False,
								  width=400, callback=plotdata_cb)
			dpg.add_slider_floatx(tag="corrscale", size=2, label='Corr Min/Max', enabled=False,
								  width=400, callback=plotdata_cb)

		with dpg.group(horizontal=True):
			with dpg.plot(tag="WaterfallPlot", height=480, width=500):
				dpg.add_plot_axis(dpg.mvXAxis, label="Time (ms)", tag="time_axis")
				dpg.add_plot_axis(dpg.mvYAxis, label="Frequency (MHz)", tag="freq_axis")
				dpg.bind_colormap("WaterfallPlot", dpg.mvPlotColormap_Viridis)
			with dpg.handler_registry():
				dpg.add_mouse_click_handler(callback=mousemask_cb)

			with dpg.plot(tag="Corr2dPlot", height=480, width=500):
				dpg.add_plot_axis(dpg.mvXAxis, label="Time lag (ms)", tag="time_lag_axis")
				dpg.add_plot_axis(dpg.mvYAxis, label="Frequency lag (MHz)", tag="freq_lag_axis")
				dpg.bind_colormap("Corr2dPlot", 9) # Pink. "Hot" is good too

		with dpg.plot(tag="TimeSeriesPlot", height=200, width=500):
			dpg.add_plot_axis(dpg.mvXAxis, label="Time (ms)", tag="tseries_t_axis")
			dpg.add_plot_axis(dpg.mvYAxis, label="Intensity (arb.)", tag="tseries_int_axis")

	### Main Menu Bar
	with dpg.window(tag="main"):
		with dpg.menu_bar(tag="MenuBar"):
			with dpg.menu(label="Menu"):
				pass
			with dpg.menu(label="Themes"):
				dpg.add_menu_item(label="Default", callback = lambda sender, data: dpg.bind_theme(global_theme))
				dpg.add_menu_item(label="Dark", callback = lambda sender, data: dpg.bind_theme(darktheme))
				dpg.add_menu_item(label="Light", callback = lambda sender, data: dpg.bind_theme(lighttheme))
			with dpg.menu(label="Tools"):
				dpg.add_menu_item(label="Show Logger", callback=logger.mvLogger)
				dpg.add_menu_item(label="Show About", callback=dpg.show_about)
				dpg.add_menu_item(label="Show Metrics", callback=dpg.show_metrics)
				dpg.add_menu_item(label="Show Documentation", callback=dpg.show_documentation)
				dpg.add_menu_item(label="Show Debug", callback=dpg.show_debug)
				dpg.add_menu_item(label="Show Style Editor", callback=dpg.show_style_editor)
				dpg.add_menu_item(label="Show Texture Registry", callback=lambda:dpg.configure_item("TextureRegistry", show=True))
				dpg.add_menu_item(label="Show Item Registry", callback=dpg.show_item_registry)
				# dpg.add_menu_item("Show Demo", callback=dpg.show_demo)

	# dpg.show_documentation()

	# Load defaults
	dpg.set_value('Filter', filefilter)
	directory_cb('user', {'file_path_name': datadir})
	if len(gdata['files']) != 0:
		burstselect_cb('burstselect', None)
		importmask_cb('user', maskfile)
		burstselect_cb('burstselect', None) # silly trick to load regions on start

		## dm range defaults
		if dmrange: dpg.set_value('dmrange', dmrange)
		dpg.set_value('numtrials', numtrials)
		dpg.set_value('dmstep', dmstep)
		dmrange_cb('user', None)

	# loadresults_cb('user2', ['B:\\dev\\frbrepeaters\\results', 'FRB121102_oostrum_525rows_Aug28.csv'])
	dpg.set_global_font_scale(1)
	dpg.setup_dearpygui()
	dpg.show_viewport()
	dpg.set_primary_window('main', True)

	# Start manual render loop. Matplotlib doesn't like running in threads so we will
	# pause DPG and run the pdf export code on the main thread when needed.
	# for the corr2dtexture matplotlib runs fine so long as a non-interactive backend is used
	# like 'agg'. This workaround however doesn't seem to work for pdf export,
	# possibly because of the additional use of the PDF backend.
	while dpg.is_dearpygui_running():
		if exportPDF:
			dpg.configure_item('ExportPDFText', show=True)
			dpg.set_value('ExportPDFText', 'Saving PDF...')
			dpg.render_dearpygui_frame() # update text
			exportresults_cb("MainThread", None)
			exportPDF = False
		dpg.render_dearpygui_frame()
	# dpg.start_dearpygui() # is replaced by above block

	dpg.destroy_context()

	return gdata

def main(): # windows cli
	frbgui(datadir=os.getcwd())

if __name__ == '__main__':
	frbgui(
		datadir='/Users/mchamma/dev/SurveyFRB20121102A/data/scholz2016',
		# maskfile='B:\\dev\\frbgui\\SurveyFRB20121102A\\aggarwalmasks_Jun17_2022.npy',
		dmstep=0.5,
		dmrange=[555, 575],
	)
