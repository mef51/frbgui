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

def log_cb(sender, data):
	dpg.log_debug(f"{sender}, {data}")

def loaddata_cb(sender, fileidx):
	filename = gdata['files'][fileidx]
	gdata['currfile'] = filename
	if gdata['currfile'] not in gdata['masks'].keys():
		gdata['masks'][gdata['currfile']] = [] # initialize list

	# filename = 'data/oostrum2020/R1_frb121102/R1_B07.rf'
	# _, pkidx, wfall = frbrepeaters.loadpsrfits(filename)
	subfall = np.load(filename)
	pkidx = np.nanargmax(np.nanmean(subfall, axis=0))
	# gdata['subfall'] = subfall
	gdata['pkidx']   = pkidx
	gdata['wfall']   = subfall

	plotdata_cb(sender, None)

def plotdata_cb(sender, scale):
	print('in plotdata')
	wfall, pkidx = gdata['wfall'], gdata['pkidx']
	width = 150
	wfall = wfall[..., pkidx-150:pkidx+150]
	# subfall = driftrate.subsample(wfall, 128, wfall.shape[1]//4)
	print('zeroing channels for ', gdata['currfile'], gdata['masks'][gdata['currfile']])
	for mask in gdata['masks'][gdata['currfile']]:
		wfall[mask] = 0

	# enable scale slider
	mostmin, mostmax = np.min(wfall), np.max(wfall)
	dpg.configure_item('wfallscale', enabled=True, min_value=mostmin, max_value=mostmax,
						format='%.{}f'.format(getscale(mostmin, mostmax)+1))

	if scale != None:
		smin, smax = scale
	else:
		smin, smax = mostmin, mostmax
	dpg.set_value('wfallscale', [smin, smax])

	dpg.add_heat_series("WaterfallPlot", "Waterfall",
		values=list(np.flipud(wfall).flatten()),
		rows=wfall.shape[0], columns=wfall.shape[1],
		scale_min=smin, scale_max=smax,
		bounds_min=(0,0), bounds_max=(wfall.shape[1], wfall.shape[0]), format='')

def subsample_cb(sender, data):
	df, dt = dpg.get_value("num freq"), dpg.get_value("num time")
	log_cb('subsample_cb', (df, dt))

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
	loaddata_cb(sender, fileidx)
	log_cb(sender, 'Opening file {}'.format(gdata['files'][fileidx]))

def masktable_cb(sender, data):
	# dpg makes working with tables impossible so we will delete the table and re-add it every time
	dpg.delete_item('Masktable')
	dpg.add_table('Masktable', [], height=100, parent='Masking')
	columns = [s.split('.')[0][-8:] for s in gdata['masks'].keys()]
	for key, col in zip(gdata['masks'].keys(), columns):
		dpg.add_column('Masktable', col, gdata['masks'][key])

def mousepos_cb(sender, data):
	isOnWaterfall = dpg.is_item_hovered('WaterfallPlot')
	if isOnWaterfall:
		tchan, fchan = dpg.get_plot_mouse_pos()
		gdata['masks'][gdata['currfile']].append(round(fchan)-1)
		plotdata_cb(sender, None)
		masktable_cb(sender, None)
		log_cb('mousepos_cb ', [[tchan, fchan], isOnWaterfall])
	else:
		return

### Analysis window
with dpg.window('FRB Analysis', width=560, height=600, x_pos=100, y_pos=50):
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
			dpg.add_table('Masktable', [], height=100)

		with dpg.tree_node('Subsampling', default_open=True):
			dpg.add_input_int("num freq", width=100, label="df", callback=subsample_cb)
			dpg.add_same_line()
			dpg.add_input_int("num time", width=100, label="dt", callback=subsample_cb)

### Plotting window
with dpg.window("FRB Plots", width=570, height=480, x_pos=900, y_pos=50):
	dpg.set_mouse_click_callback(mousepos_cb)
	dpg.add_slider_float2("wfallscale", enabled=False, callback=plotdata_cb,
		callback_data=lambda: dpg.get_value('wfallscale'))
	dpg.add_plot("WaterfallPlot", no_legend=True)
	dpg.set_color_map("WaterfallPlot", 5) # Viridis

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
