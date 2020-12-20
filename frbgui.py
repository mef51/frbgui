from dearpygui import core, simple, demo

demo.show_demo()

def save_callback(sender, data):
    print("Save Clicked")

with simple.window("Example Window"):
    core.add_text("Hello world")
    core.add_button("Save", callback=save_callback)
    core.add_input_text("string")
    core.add_slider_float("float")

core.start_dearpygui()
