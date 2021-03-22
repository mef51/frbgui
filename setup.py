import setuptools

setuptools.setup(
	name="frbrepeaters", # Replace with your own username
	version="0.0.1",
	author="Mohammed Chamma",
	author_email="mchamma@uwo.ca",
	description="utilities for processing FRBs from repeating sources",
	# long_description=long_description,
	# long_description_content_type="text/markdown",
	# url="https://github.com/pypa/sampleproject",
	packages=setuptools.find_packages(),
	classifiers=[
		"Programming Language :: Python :: 3",
		# "License :: OSI Approved :: MIT License",
		"Operating System :: OS Independent",
	],
	python_requires='>=3.6',
)
