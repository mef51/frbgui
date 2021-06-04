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
datastore = {}
# dpg.show_debug()
dpg.show_logger()

def log_cb(sender, data):
	dpg.log_debug(f"{sender}, {data}")

def loaddata_cb(sender, fileidx):
	filename = datastore['files'][fileidx]
	# filename = 'data/oostrum2020/R1_frb121102/R1_B07.rf'
	print('in load callback')
	# _, pkidx, wfall = frbrepeaters.loadpsrfits(filename)
	subfall = np.load(filename)
	pkidx = np.nanargmax(np.nanmean(subfall, axis=0))
	# datastore['subfall'] = subfall
	datastore['pkidx']   = pkidx
	datastore['wfall']   = subfall
	plotdata_cb(sender, fileidx)

def plotdata_cb(sender, data):
	wfall, pkidx = datastore['wfall'], datastore['pkidx']
	width = 150
	wfall = wfall[..., pkidx-150:pkidx+150]
	# subfall = driftrate.subsample(wfall, 128, wfall.shape[1]//4)
	dpg.add_heat_series("Plot", "Waterfall",
		values=list(np.flipud(wfall).flatten()),
		rows=wfall.shape[0], columns=wfall.shape[1],
		scale_min=np.min(wfall), scale_max=np.max(wfall),
		bounds_min=(0,0), bounds_max=(wfall.shape[1], wfall.shape[0]), format='')

def subsample_cb(sender, data):
	df, dt = dpg.get_value("num freq"), dpg.get_value("num time")
	log_cb('subsample_cb', (df, dt))

def directory_cb(sender, data):
	dpg.set_value('Dirtext', 'Selected: {}'.format(data[0]))
	files = glob.glob(data[0]+'/*')
	dpg.configure_item('burstselect', items=[os.path.basename(x) for x in files])
	datastore['datadir'] = data[0]
	datastore['files']   = files
	log_cb(sender, data[0])

def burstselect_cb(sender, data):
	fileidx = dpg.get_value('burstselect')
	loaddata_cb(sender, fileidx)
	log_cb(sender, 'Opening file {}'.format(datastore['files'][fileidx]))

with dpg.window("frb plots", width=570, height=480, x_pos=100, y_pos=50):
	dpg.add_plot("Plot", no_legend=True)
	dpg.set_color_map("Plot", 5) # Viridis

with dpg.window('frb analysis', width=560, height=600, x_pos=900, y_pos=50):
	with dpg.collapsing_header("1. Data", default_open=True):
		dpg.add_button("Select Directory...", callback = lambda s, d: dpg.select_directory_dialog(directory_cb))
		dpg.add_same_line()
		dpg.add_text("Dirtext", default_value="Selected: ")
		dpg.add_text("Files found:")
		dpg.add_listbox("burstselect", label='', items=[], num_items=10, width=520,
						callback=burstselect_cb, tip="Select to load burst...")
	with dpg.collapsing_header("2. Subsampling", default_open=True):
		dpg.add_input_int("num freq", width=100, label="df", callback=subsample_cb)
		dpg.add_same_line()
		dpg.add_input_int("num time", width=100, label="dt", callback=subsample_cb)
		# dpg.add_button("Load and Plot Data", callback=loaddata_cb, tip="This can block the interface for large data files")

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
dpg.start_dearpygui(primary_window='main')
