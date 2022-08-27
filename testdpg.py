import dpg
dpg.set_main_window_size(500, 500)
dpg.set_main_window_title("Group Test")

with dpg.window('FRB Analysis', width=200, height=200, x_pos=10, y_pos=30):
	with dpg.group("Hello"):
		pass

	with dpg.group("Bye", parent="Hello"):
		dpg.add_button("A button", parent="Hello")
		dpg.add_button("B button", parent="Hello")
		dpg.add_button("C button", parent="Hello")

dpg.delete_item("Hello", children_only=True)
dpg.start_dearpygui()
