from dearpygui.core import *
from dearpygui.simple import *
import time
import numpy as np

def slow_callback(sender, data):
	print('in slow_cb')
	x, y = 5, 5
	s = np.random.random((x, y))
	for i in range(0, 4096):
		s = s + np.random.random((x,y))
	# print("sum:", s.sum())

def async_cb(s, d):
	print('in async_cb')
	set_threadpool_high_performance()
	t1 = time.time()
	run_async_function(slow_callback, None, return_handler=lambda s, d: print('dpg thread time:', time.time() - t1, flush=True))

t1 = time.time()
slow_callback(None, None)
t2 = time.time()
print('main thread time:', t2-t1)

with window("async test"):
	add_button("test async", callback=async_cb)
	add_button("test sync", callback=slow_callback)

show_debug()
# show_demo()
start_dearpygui()
