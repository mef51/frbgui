import dpg
import frbrepeaters
import matplotlib.pyplot as plt
import numpy as np

subfall, pkidx = frbrepeaters.loadpsrfits('data/oostrum2020/R1_frb121102/R1_B07.rf')
width = 150
subfall = subfall[:, pkidx-width:pkidx+width]
# plt.imshow(subfall[:, pkidx-width:pkidx+width], origin='lower',interpolation='none', aspect='auto')
# plt.show()
# dpg.show_demo()

def log_callback(sender, data):
	dpg.log_debug(f"{sender}, {data}")

with dpg.window("frb plots", width=750, height=280, x_pos=100, y_pos=50):
	dpg.add_plot("Plot", no_legend=True)
	dpg.set_color_map("Plot", 5) # Viridis
	dpg.add_heat_series("Plot", "DynamicSpectrum",
		values=list(np.flipud(subfall).flatten()),
		rows=subfall.shape[0], columns=subfall.shape[1],
		scale_min=np.min(subfall), scale_max=np.max(subfall),
		bounds_min=(0,0), bounds_max=(subfall.shape[1], subfall.shape[0]), format='')

with dpg.window('frb analysis', width=560, x_pos=900, y_pos=50):
	with dpg.collapsing_header("Subsampling", default_open=True):
		dpg.add_input_int("num freq", callback=log_callback, width=100, label="df")
		dpg.add_same_line()
		dpg.add_input_int("num time", callback=log_callback, width=100, label="dt")

### Main Menu Bar
with dpg.window("main"):
	with dpg.menu_bar("MenuBar##demo"):
		with dpg.menu("Menu##demo"):
			pass
		with dpg.menu("Themes##demo"):
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
		with dpg.menu("Tools##demo"):
			dpg.add_menu_item("Show Logger##demo", callback=dpg.show_logger)
			dpg.add_menu_item("Show About##demo", callback=dpg.show_about)
			dpg.add_menu_item("Show Metrics##demo", callback=dpg.show_metrics)
			dpg.add_menu_item("Show Documentation##demo", callback=dpg.show_documentation)
			dpg.add_menu_item("Show Debug##demo", callback=dpg.show_debug)
			dpg.add_menu_item("Show Style Editor##demo", callback=dpg.show_style_editor)


# dpg.show_documentation()
dpg.set_main_window_size(1600, 800)
dpg.set_main_window_title("FRB Repeaters")
# dpg.start_dearpygui()
dpg.start_dearpygui(primary_window='main')
