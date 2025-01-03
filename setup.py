import setuptools

def readme():
	with open("README.md") as f:
		return f.read()

setuptools.setup(
	name="frbgui",
	version="0.11.0",
	author="Mohammed Chamma",
	author_email="mchamma@uwo.ca",
	description="GUI and utilities for processing Fast Radio Burst waterfalls",
	long_description=readme(),
	long_description_content_type="text/markdown",
	url="https://github.com/mef51/frbgui",
	#packages=setuptools.find_packages(),
	py_modules=['frbgui', 'driftrate', 'driftlaw', 'arrivaltimes'],
	install_requires=[
		'matplotlib>=3.7.2',
		'numpy>=1.24.4',
		'dearpygui>=1.9.1',
		'dearpygui_ext>=0.9.5',
		'your>=0.6.7',
		'pandas>=2.1.0',
		'tqdm>=4.65.0'
	],
	license='MIT',
	classifiers=[
		"Programming Language :: Python :: 3",
		"License :: OSI Approved :: MIT License",
		"Operating System :: OS Independent",
		"Development Status :: 3 - Alpha",
		"Intended Audience :: Science/Research",
		"Topic :: Scientific/Engineering :: Astronomy"
	],
	python_requires='>=3.6',
	scripts=['bin/frbgui'],
	entry_points={
		'console_scripts': ['frbgui=frbgui:main']
	}
)
