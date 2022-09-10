import setuptools

setuptools.setup(
	name="frbgui",
	version="0.9",
	author="Mohammed Chamma",
	author_email="mchamma@uwo.ca",
	description="GUI and utilities for processing Fast Radio Burst waterfalls",
	# long_description=long_description,
	# long_description_content_type="text/markdown",
	url="https://github.com/mef51/frbgui",
	#packages=setuptools.find_packages(),
	py_modules=['frbgui', 'driftrate', 'driftlaw'],
	install_requires=[
		'dearpygui==0.6.415',
		'your>=0.6.5',
		'pandas'
	],
	classifiers=[
		"Programming Language :: Python :: 3",
		# "License :: OSI Approved :: MIT License",
		"Operating System :: OS Independent",
	],
	python_requires='>=3.6',
	scripts=['bin/frbgui']
)
