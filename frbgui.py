import dpg
import frbrepeaters
import driftrate
import os, glob, itertools
import matplotlib.pyplot as plt
import numpy as np

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
	'globfilter'     : '*.npy',
	'masks'          : {}, # will store lists of masked channel numbers
	'datadir'        : 'B:\\dev\\frbrepeaters\\data\\luo2020\\180813_ar_file\\ar_file\\converted'
}

def getscale(m, M):
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

def log_cb(sender, data):
	dpg.log_debug(f"{sender}, {data}")
def error_log_cb(sender, data):
	dpg.log_error(f"{sender}, {data}")

def loaddata_cb(sender, data):
	wfall = None
	if 'fileidx' in data.keys():
		filename = gdata['files'][data['fileidx']]
		gdata['currfile'] = filename
		wfall = np.load(filename)
		gdata['wfall']          = wfall
		# cache the original waterfall, since subsample needs it
		gdata['wfall_original'] = np.copy(wfall)

		# update subsample controls
		dpg.set_value('Wfallshapelbl', 'Maximum Size: {}'.format(np.shape(wfall)))
		dpg.set_value('Subfallshapelbl', 'Current Size: {}'.format(np.shape(wfall)))
		dpg.configure_item('dfreqinput', enabled=True, min_value=0, max_value=wfall.shape[0])
		dpg.configure_item('dtimeinput', enabled=True, min_value=0, max_value=wfall.shape[1])
		dpg.set_value('dfreqinput', wfall.shape[0])
		dpg.set_value('dtimeinput', wfall.shape[1])
	elif 'subsample' in data.keys() and data['subsample']:
		wfall = gdata['wfall']
		dpg.set_value('Subfallshapelbl', 'Current Size: {}'.format(np.shape(wfall)))
	else:
		wfall = gdata['wfall']

	if gdata['currfile'] not in gdata['masks'].keys():
		gdata['masks'][gdata['currfile']] = [] # initialize list

	if wfall.shape == gdata['wfall_original'].shape:
		wfall = applyMasks(wfall)

	gdata['ts']      = np.nanmean(wfall, axis=0)
	pkidx = np.nanargmax(gdata['ts'])
	gdata['pkidx']   = pkidx

	plotdata_cb(sender, data)

def plotdata_cb(sender, data):
	if not data:
		data = {}

	wfall, pkidx = gdata['wfall'], gdata['pkidx']
	twidth = 150
	wfall = wfall[..., pkidx-twidth:pkidx+twidth]

	# print('zeroing channels for ', gdata['currfile'], gdata['masks'][gdata['currfile']])
	for mask in gdata['masks'][gdata['currfile']]:
		if mask < len(wfall):
			wfall[mask] = 0

	corr = driftrate.autocorr2d(wfall)

	## enable scale sliders
	mostmin, mostmax = np.min(wfall), np.max(wfall)
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
		values=list(np.flipud(wfall).flatten()),
		rows=wfall.shape[0], columns=wfall.shape[1],
		scale_min=smin, scale_max=smax,
		bounds_min=(0,0), bounds_max=(wfall.shape[1], wfall.shape[0]), format='')

	dpg.add_heat_series("Corr2dPlot", "Corr2d",
		values=list(np.flipud(corr).flatten()),
		rows=corr.shape[0], columns=corr.shape[1],
		scale_min=scmin, scale_max=scmax,
		bounds_min=(0,0), bounds_max=(corr.shape[1], corr.shape[0]), format='')

	tseries = gdata['ts'][pkidx-twidth:pkidx+twidth]
	dpg.add_line_series("TimeSeriesPlot", "TimeSeries",
						list(range(0, len(tseries))), tseries)

	if 'keepview' in data.keys() and data['keepview']:
		dpg.set_plot_xlimits('WaterfallPlot', wx[0], wx[1])
		dpg.set_plot_ylimits('WaterfallPlot', wy[0], wy[1])
		dpg.set_plot_xlimits_auto('WaterfallPlot')
		dpg.set_plot_ylimits_auto('WaterfallPlot')

def subsample_cb(sender, data):
	df, dt = dpg.get_value("dfreqinput"), dpg.get_value("dtimeinput")
	print(df, dt)

	try:
		# Make a copy of the original fall, apply the masks, then downsample
		wfall = np.copy(gdata['wfall_original'])
		wfall = applyMasks(wfall)
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

def masktable_cb(sender, data):
	# dpg makes working with tables impossible so we will delete the table and re-add it every time
	dpg.delete_item('Masktable')

	tableheight = round(min(25*len(gdata['masks'][list(gdata['masks'].keys())[0]]), 250))
	dpg.add_table('Masktable', [], height=tableheight, parent='Masking')

	columns = [s.split('.')[0][-8:] for s in gdata['masks'].keys()]
	for key, col in zip(gdata['masks'].keys(), columns):
		dpg.add_column('Masktable', col, gdata['masks'][key])

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

### Analysis window
with dpg.window('FRB Analysis', width=560, height=700, x_pos=10, y_pos=30):
	with dpg.collapsing_header("1. Data", default_open=True):
		dpg.add_button("Select Directory...", callback = lambda s, d: dpg.select_directory_dialog(directory_cb))
		dpg.add_text("Dirtext", default_value="Selected: (no directory selected)")
		dpg.add_text("Filter:"); dpg.add_same_line()
		dpg.add_input_text("Filter", label='', hint="eg. *.npy", callback=filter_cb, enabled=False)
		dpg.add_same_line();
		dpg.add_button('Clear filter', callback=clearfilter_cb, enabled=False)

		dpg.add_text("Files found:")
		dpg.add_listbox("burstselect", label='', items=[], num_items=10, width=520,
						callback=burstselect_cb, tip="Select to load burst...")
	with dpg.collapsing_header("2. Burst Cleanup", default_open=True):
		with dpg.tree_node('Masking', default_open=True):
			dpg.add_text("Click on the waterfall plot to begin masking frequency channels.")
			dpg.add_text("NOTE: only mask on the original waterfall (todo: add a 'mask' button)")
			dpg.add_table('Masktable', [], height=50)

		with dpg.tree_node('Downsampling', default_open=True):
			dpg.add_text("Wfallshapelbl", default_value="Maximum Size: (no burst selected)")
			dpg.add_text("Subfallshapelbl", default_value="Current Size: (no burst selected)")
			dpg.add_input_int("dfreqinput", width=100, label="df", callback=subsample_cb, enabled=False)
			dpg.add_same_line()
			dpg.add_input_int("dtimeinput", width=100, label="dt", callback=subsample_cb, enabled=False)

### Plotting window
with dpg.window("FRB Plots", width=1035, height=745, x_pos=560, y_pos=30):
	dpg.set_mouse_click_callback(mousepos_cb)
	dpg.add_slider_float2("wfallscale", label='Wfall Min/Max', enabled=False,
						  width=400, callback=plotdata_cb,
						  callback_data=lambda: {'scale': dpg.get_value('wfallscale')})
	dpg.add_same_line()
	dpg.add_slider_float2("corrscale", label='Corr Min/Max', enabled=False,
						  width=400, callback=plotdata_cb,
						  callback_data=lambda: {'cscale': dpg.get_value('corrscale')})

	# Colors: From implot.cpp: {"Default","Deep","Dark","Pastel","Paired","Viridis","Plasma","Hot","Cool","Pink","Jet"};
	dpg.add_plot("WaterfallPlot", no_legend=True, height=480, width=500)
	dpg.set_color_map("WaterfallPlot", 5) # Viridis
	dpg.add_same_line()
	dpg.add_plot("Corr2dPlot", no_legend=True, height=480, width=500)
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


# dpg.show_documentation()
dpg.set_main_window_size(1600, 800)
dpg.set_main_window_title("FRB Repeaters")
# dpg.start_dearpygui()

# Load defaults
dpg.set_value('Filter', gdata['globfilter'])
directory_cb('user', [gdata['datadir']]) # default directory

dpg.start_dearpygui(primary_window='main')
