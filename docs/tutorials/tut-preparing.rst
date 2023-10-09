.. _tutpreparation:

Preparing Data
==============

Fast Radio Burst observations are taken with different telescopes from different countries around the world. Because of this, different studies will make their radio observations available in all kinds of formats and knowing how to adapt to and deal with different formats is a skill that one invariably develops when working with FRBs (and astronomical data in general).

This tutorial is intended as an example of how FRB data from different sources can be prepared for use in FRBGui and is intended for researchers new to FRBs with some technical experience programming in Python. Parts of the tutorial will require FRBGui to be installed (see :ref:`installation`).

In this tutorial we will load FRB data from two different studies and prepare them in FRBGui's :ref:`burstformat` and save them as numpy zipped ``.npz`` archive files.

These ``.npz`` files will be used in the tutorial ":ref:`tutmeasure`" to obtain spectro-temporal measurements of the FRBs. These tutorials can be done out of order if preferred, as the ``.npz`` files will be provided in the next tutorial for download.

The two studies we will process data from today will feature two different data formats and are

* `Gajjar et al. (2018) <https://iopscience.iop.org/article/10.3847/1538-4357/aad005>`_:
* `Aggarwal et al. (2021) <https://iopscience.iop.org/article/10.3847/1538-4357/ac2577>`_


Gajjar et al. (2018) features observations of the repeating source FRB 20121102A using the Green Bank Telescope. As stated in their paper, the data are available in the `PSRFITS <https://www.atnf.csiro.au/research/pulsar/psrfits_definition/Psrfits.html>`_ format at http://seti.berkeley.edu/frb121102/technical.html. The bursts in these data are already conveniently cut-out into individual files.

Aggarwal et al. (2021) also features observations of the repeating source FRB 20121102A, taken with the Arecibo telescope. These data are available in the `filterbank <https://sigproc.sourceforge.net/sigproc.pdf>`_ format, and their results are available at https://github.com/thepetabyteproject/FRB121102. For these data, we will download the entire observation and need to locate the bursts and cut them out ourselves.

Following each paper's respective documentation, we will use Python packages that can read the data and prepare scripts to read the PSRFITS and filterbank data and save the bursts into ``.npz`` files that we can use with FRBGui.

Our goal is to prepare the following two bursts so that we may obtain their spectro-temporal properties in the next tutorial.

.. _paperfigure:

.. figure:: /imgs/tut-prep-1.png

	On the left, burst 11A from Fig. 2 of Gajjar et al. (2018). On the right, Burst B6 from Fig. 1 of Aggarwal et al. (2021).

Reading PSRFITS (``.fits``)
---------------------------

In this section we will develop a script that reads a PSRFITS file and prepares the FRB contained within it for measurement.

**>>** Open the `documentation <http://seti.berkeley.edu/frb121102/technical.html>`_ for Gajjar et al. (2018) and download burst 11A which has the filename ``11A_16sec.calib.4p`` and is 308 MB. Save this file into your folder of choice.

As stated in their documentation, the package `PyPulse <https://github.com/mtlam/PyPulse>`_ can be used to read their PSRFITS format.

**>>** Install PyPulse from the command line

.. code-block:: bash

	pip install --user pypulse

.. note::

	Note that version 0.1.1 of PyPulse will not work with recent versions of numpy (>=1.20).

	If a newer release of PyPulse is not yet available on pip, you will need to install PyPulse from a local copy of its repository which has resolved its issue with newer versions of numpy. You can do so in the following way:

	.. code-block:: bash

		git clone https://github.com/mtlam/PyPulse.git
		cd PyPulse
		pip install --user --editable .

With PyPulse installed you can open the data and start preparing it according to FRBGui's :ref:`burstformat` with the following code:

.. code-block:: python
	:linenos:

	import pypulse
	import numpy as np
	import matplotlib.pyplot as plt

	ar = pypulse.Archive('11A_16sec.calib.4p', prepare=False) # load without dedispersing
	ar.pscrunch() # average polarizations if any
	ar.center() # center pulse in array

	wfall = ar.getData()
	burstmetadata = {
		'dt'        : ar.getTimes(),
		'dfs'       : ar.getFreqs(),
		'DM'        : ar.getDM(),
		'bandwidth' : ar.getBandwidth(),
		'duration'  : ar.getDuration(), # usually in seconds
		'center_f'  : ar.getCenterFrequency(),
		'freq_unit' : ar.getFrequencyUnit(),
		'time_unit' : ar.getTimeUnit(),
		'int_unit'  : ar.getIntensityUnit(),
		'telescope' : ar.getTelescope(),
		'burstSN'   : ar.getSN(),
		'tbin'      : ar.getTbin(),
	}

In the above we have used the methods in ``pypulse.Archive`` to extract the information we will need from the fits file. For more information on these methods visit `PyPulse's documentation <https://mtlam.github.io/PyPulse/archive.html>`_.

We can now perform a few checks to ensure the values from the fits file will work with FRBGui and adjust them if not.

.. code-block:: python
	:lineno-start: 25

	for item in burstmetadata.items():
		print(*item)

Which outputs each key and value in ``burstmetadata``:

.. code-block:: python

	dt [0.01048303]
	dfs [8188.87304688 8188.68994141 8188.50683594 ... 4626.92236328 4626.73925781
	 4626.55615234]
	DM 557.91
	bandwidth -3562.5
	duration 0.020966058666666648
	center_f 6407.7148
	freq_unit MHz
	time_unit SEC
	int_unit Jy
	telescope GBT
	burstSN 7.9028664
	tbin 1.0239731655700501e-05

We can also see that ``wfall.shape = (19456, 2048)``, indicating the data have 19456 frequency channels and 2048 time channels.

These values align well with FRBGui's :ref:`burstformat` since the ``freq_unit`` and ``time_unit`` are already in MHz and seconds, respectively. However the ``dfs`` array and the ``bandwidth`` field are in slightly different formats. The ``dfs`` array should be sorted from smallest to largest and the negative sign in ``bandwidth`` is ignored by FRBGui, so it is a good idea to remove it. We can correct these by doing

.. code-block:: python
	:lineno-start: 28

	burstmetadata['bandwidth'] = abs(burstmetadata['bandwidth'])
	burstmetadata['dfs'].sort()

At this stage we could save the wfall and burstmetadata into an .npz file and the file will load in FRBGui.

However, due to the large size of the data array (with 19456 frequency channels) we may want to reduce the amount of data that is loaded into FRBGui for the sake of signal to noise as well as performance. We can also plot the data and check that we have correctly read the file.

Let's start by plotting the data using matplotlib to see what we are working with.

.. code-block:: python
	:lineno-start: 31

	plt.imshow(wfall, aspect='auto', origin='lower')

We use ``aspect='auto'`` to get a more square figure and ``origin='lower'`` to specify we want the lowest frequency at the bottom. Feel free to experiment with these options.

The resulting plot is shown below.

.. figure:: /imgs/tut-prep-2.png
	:width: 450
	:align: center

	Burst waterfall (frequency on y-axis, time on x-axis)

Oops! We don't see anything. The large number of channels may be washing out our signal, so let's temporarily downsample the burst using the :py:mod:`driftrate` module (which ships with FRBGui) to see if the signal to noise will increase.

.. code-block:: python
	:lineno-start: 31
	:emphasize-lines: 3

	import driftrate

	wfall = driftrate.subsample(wfall, 304, 512)
	plt.imshow(wfall, aspect='auto', origin='lower')

Here we've downsampled the original waterfall size of 19456 freq. channels and 2048 time channels to 19456/64 = 304 frequency channels and 2048/4 = 512 time channels.

.. figure:: /imgs/tut-prep-3.png
	:width: 450
	:align: center

	Downsampled burst waterfall (frequency on y-axis, time on x-axis)

While still faint we can begin to see the burst (though it appears flipped) and more importantly what appears to be two bands of radio frequency interference (RFI) near the bottom and top of the waterfall at around y = 20 and 245.

These RFI bands may be responsible for washing out the burst. If we can locate them precisely we can remove them and then plot the waterfall without them.

To precisely locate them we will plot the waterfall's spectrum, averaged over all time samples, to obtain a plot of intensity vs. frequency.

.. code-block:: python

	plt.plot(np.nanmean(wfall, axis=1))

.. figure:: /imgs/tut-prep-4.png
	:width: 450
	:align: center

	Burst spectrum (intensity on y-axis, frequency on x-axis)

The two spikes at around x = 1000 and 15300 are the RFI we saw in the earlier waterfall. By inspecting the figure more closely using matplotlib's graphical interface (you can use ``plt.show()`` in a script or ``%matplotlib qt`` in a jupyter notebook) we can more precisely determine the offending channel numbers and remove them from the waterfall with the following code. We will also flip the waterfall so that the burst is the right way up.

.. code-block:: python
	:lineno-start: 32
	:emphasize-lines: 1,2,3

	wfall[1053:1110] = 0
	wfall[15360:15471] = 0
	wfall = np.flipud(wfall)
	wfall = driftrate.subsample(wfall, 304, 512)

	plt.imshow(wfall, aspect='auto', origin='lower', interpolation='none')

.. figure:: /imgs/tut-prep-5.png
	:width: 450
	:align: center

	Burst waterfall with high signal-to-noise (frequency on y-axis, time on x-axis)

Great! We can now clearly see the burst.

Computing Axes
^^^^^^^^^^^^^^

One last check we can perform is to add the units of the frequency (MHz) and time (ms) axes so that we can be sure we have loaded the data as it is presented in the paper (See the figure :ref:`above <paperfigure>`).

The PSRFITS file contains information about the bandwidth, duration, frequency axis, and duration in the ``burstmetadata['bandwidth']``, ``burstmetadata['duration']``, ``burstmetadata['dfs']``, ``burstmetadata['dt']``, and ``burstmetadata['tbin']`` fields (which we printed above). The ``dt`` and ``tbin`` fields seem to slightly differ, and if we compute the resolution with the duration and original waterfall shape (before subsampling)

.. code-block:: python

	>>> burstmetadata['duration']/wfall.shape[0]*1000
	0.010237333333333324

we get yet another slightly different value. In this case it may be best to just try all three and see which best matches up with the publication or contact the paper author for clarification.

For now we will simply compute the frequency and time resolutions using the ``'bandwidth'`` and ``'duration'`` fields from the file and the shape of the original waterfall. Downsampling the burst will change the resolutions since we are decreasing the number of channels, so we will store the original shape of the waterfall and update the resolutions based on the factor we downsampled by.

.. code-block:: python
	:lineno-start: 32
	:emphasize-lines: 5,6,7,9,10

	wfall[1053:1110] = 0
	wfall[15360:15471] = 0
	wfall = np.flipud(wfall)

	df = burstmetadata['bandwidth']/wfall.shape[0] # MHz
	dt = burstmetadata['duration']/wfall.shape[1]*1000 # ms
	origshape = wfall.shape
	wfall = driftrate.subsample(wfall, 2432, 2048)
	df *= origshape[0]/wfall.shape[0]
	dt *= origshape[1]/wfall.shape[1] # no change

	plt.imshow(wfall, aspect='auto', origin='lower', interpolation='none')

Note here that we have changed our downsampling to 2432 by 2048 channels (more than before but less than the original) to preserve some of the data resolution. This will decrease the signal-to-noise but we will be able to modify this more dynamically later from inside FRBGui. A higher data resolution also helps with measurements stability and accuracy.

Using the frequency and time resolutions, we can now add axis labels and display the waterfall axes using ``imshow``'s ``extent`` keyword, which takes a list of ``[left, right, bottom, top]`` limits for the axes.

.. code-block:: python
	:lineno-start: 42

	lowest_freq = min(burstmetadata['dfs'])
	extent = [
		0,
		dt*wfall.shape[1],
		lowest_freq,
		lowest_freq + df*wfall.shape[0]
	]
	# convenience function from driftrate module, same as above:
	# extent, _ = driftrate.getExtents(wfall, df=df, dt=dt, lowest_freq=lowest_freq)

	plt.imshow(wfall, aspect='auto', origin='lower', interpolation='none', extent=extent)
	plt.xlabel("Time (ms)")
	plt.ylabel("Frequency (MHz)")
	plt.title(f"Burst 11A ({wfall.shape = })")

.. figure:: /imgs/tut-prep-6.png
	:width: 800

	The prepared waterfall, ready for saving. On the left the burst is faintly seen compared to on the right due to the difference in waterfall downsampling.

From here we can see that the frequency axes match with the figure in the paper and that the burst lasts just over 2 ms, also consistent with what is shown in the paper figure.

There are other tasks we could perform if desired, we could crop the waterfall to include less of the data before and after the burst or we could decrease the bandwidth. This however is sufficient for our purposes and we can now save the burst for measurement in FRBGui with the following command.

.. code-block:: python
	:lineno-start: 56

	np.savez('burst11A.npz', wfall=wfall, **burstmetadata)

This creates a file ``burst11A.npz`` that can be loaded into FRBGui (see :ref:`tutmeasure`) in the directory of your script. Congratulations!

This script can now be used as the basis for preparing PSRFITS files, especially the other bursts made available from Gajjar et al. (2018). For example, it can be adapted to prepare all the bursts automatically in a loop over the downloaded data files.

Complete Code
^^^^^^^^^^^^^

For your reference, below is the complete script we developed for reading burst 11A in its PRSFITS format and preparing it as a numpy zipped ``.npz`` file, ready for further analysis in FRBGui or other Python scripts. Included are some optional lines for removing the bandwidth above and below the burst.

.. code-block:: python
	:linenos:

	import pypulse
	import numpy as np
	import matplotlib.pyplot as plt
	import driftrate

	ar = pypulse.Archive('11A_16sec.calib.4p', prepare=False) # load without dedispersing
	ar.pscrunch() # average polarizations if any
	ar.center() # center pulse in array

	wfall = ar.getData()
	burstmetadata = {
		'dt'        : ar.getTimes(),
		'dfs'       : ar.getFreqs(),
		'DM'        : ar.getDM(),
		'bandwidth' : ar.getBandwidth(),
		'duration'  : ar.getDuration(), # usually in seconds
		'center_f'  : ar.getCenterFrequency(),
		'freq_unit' : ar.getFrequencyUnit(),
		'time_unit' : ar.getTimeUnit(),
		'int_unit'  : ar.getIntensityUnit(),
		'telescope' : ar.getTelescope(),
		'burstSN'   : ar.getSN(),
		'tbin'      : ar.getTbin(),
	}

	burstmetadata['bandwidth'] = abs(burstmetadata['bandwidth'])
	burstmetadata['dfs'].sort()

	wfall[1053:1110] = 0
	wfall[15360:15471] = 0
	wfall = np.flipud(wfall)

	df = burstmetadata['bandwidth']/wfall.shape[0] # MHz
	dt = burstmetadata['duration']/wfall.shape[1]*1000 # ms
	origshape = wfall.shape
	wfall = driftrate.subsample(wfall, 2432, 2048)
	df *= origshape[0]/wfall.shape[0]
	dt *= origshape[1]/wfall.shape[1] # no change

	# Optional: crop frequency band to match Fig 1 in Gajjar+2018
	# burstmetadata['bandwidth'] = burstmetadata['bandwidth']*((wfall.shape[0] - ((wfall.shape[0]-2200)+510))/wfall.shape[0])
	# lowest_freq = burstmetadata['center_f'] - wfall.shape[0]/2*(df*(19456/2432)) + (df*(19456/2432))*510
	# wfall = wfall[510:2200, ...]
	# df = burstmetadata['bandwidth']/wfall.shape[0]
	# burstmetadata['dfs'] = np.linspace(lowest_freq, lowest_freq+burstmetadata['bandwidth'], num=wfall.shape[0])

	lowest_freq = min(burstmetadata['dfs'])
	extent = [
		0,
		dt*wfall.shape[1],
		lowest_freq,
		lowest_freq + df*wfall.shape[0]
	]
	# convenience function from driftrate module, same as above:
	# extent, _ = driftrate.getExtents(wfall, df=df, dt=dt, lowest_freq=lowest_freq)

	plt.imshow(wfall, aspect='auto', origin='lower', interpolation='none', extent=extent)
	plt.xlabel("Time (ms)")
	plt.ylabel("Frequency (MHz)")
	plt.title(f"Burst 11A ({wfall.shape = })")

	np.savez('burst11A.npz', wfall=wfall, **burstmetadata)


Reading Filterbank (``.fil``)
-----------------------------

Our goal in this section is to read the bursts from Aggarwal et al. (2021), specifically burst B6 and prepare if for measurement in FRBGui.

The data are available in a large filterbank file that includes the entire observational scan. It's size is multiple gigabytes. Accompanying the paper is a spreadsheet with the timestamps of when the bursts were observed in the filterbank file. The paper also explains that it used a custom package called `BurstFit <https://github.com/thepetabyteproject/burstfit/>`_ for its analysis. Our strategy then will be to download the files, extract the timestamps from the paper's provided spreadsheet, and use BurstFit to cutout the data from the filterbank file. In order to extract basic information about the data file such as the resolutions and bandwidth, we will use another package called `Your <https://github.com/thepetabyteproject/your>`_ for convenience.

**>>** Download the `filterbank data <https://zenodo.org/record/5029530>`_ from Aggarwal et al. (2021) and extract it (6.2 GB). In it will be two filterbank files in deeply nested folders. Copy them to the directory of your choice and start a script.

.. code-block:: python
	:linenos:

	files = [
		'puppi_57644_C0531+33_0021_subs_0001.fil',
		'puppi_57645_C0531+33_0029_subs_0001.fil'
	]

Note that the '57644' and '57645' in the filenames refer to the day of the observation (the Modified Julian Date).

**>>** Download `all_bursts_bary.csv <https://github.com/thepetabyteproject/FRB121102/blob/main/data/all_bursts_bary.csv>`_ which will contain metadata about the bursts as well as the burst timestamps.

**>>** Install the `BurstFit <https://github.com/thepetabyteproject/burstfit/>`_ and `Your <https://github.com/thepetabyteproject/your>`_ packages on the command-line. Read the `"BurstData" class section <https://thepetabyteproject.github.io/burstfit/BurstData/#burstdata-class>`_ of the BurstFit documentation.

.. code-block:: bash

 	pip install --user burstfit your

**>>** Load the spreadsheet using pandas (``pip install --user pandas`` if you do not already have it)

.. code-block:: python
	:linenos:

	import pandas as pd

	canddf = pd.read_csv('all_bursts_bary.csv')

The spreadsheet does not have an explicit column for the "candidate times" but it does contain all the information for locating each burst in the "cand_id" column. For example,

.. code-block:: python

	>>> canddf['cand_id'][0]
	'cand_tstart_57644.407719907409_tcand_61.2631000_dm_565.30000_snr_8.12529'

In this field we can find the day of the observation, the candidate time (the time from the start of the filterbank file where the candidate burst can be found), the DM, and the signal-to-noise ratio (SNR). Notice that the values are separated by a "_" character.

We will use the ``BurstData`` class from the BurstFit package to load the bursts. This class requires filename, DM, candidate time, width (number of samples to load), and the SNR of the burst.

**>>** Extract the properties needed for the ``BurstData`` class from the ``'cand_id'`` column and add these columns to the spreadsheet

.. code-block:: python
	:lineno-start: 4

	canddf['tstart'] = [float(candid.split('_')[2]) for candid in canddf['cand_id']]
	canddf['tcand']  = [float(candid.split('_')[4]) for candid in canddf['cand_id']]
	canddf['dm']     = [float(candid.split('_')[6]) for candid in canddf['cand_id']]
	canddf['snr']    = [float(candid.split('_')[8]) for candid in canddf['cand_id']]
	canddf['label']  = canddf['bidx']
	canddf['width']  = 256
	canddf['file']   = [files[0] if '57644' in str(tstart) else files[1] for tstart in canddf['tstart']]
	canddf = canddf.set_index('label')

Here we have used the ``.split()`` function and the fact that ``'cand_id'`` contains the "_" character between values to extract the information we want. We have also used the ``'bidx'`` column to serve as a label for each burst, and chosen to extract 256 time samples around each candidate burst.

At this point we could just load the single burst we want, but we have done all this work and it is perfectly setup to just extract every burst from the study in a nice for-loop. So let's do that.


.. code-block:: python
	:lineno-start: 13

	import numpy as np
	import matplotlib.pyplot as plt
	from burstfit.data import BurstData

	for bid, row in canddf.iterrows():
		bd = BurstData(
			fp=row['file'],
			dm=row['dm'],
			width=row['width'],
			snr=row['snr'],
			tcand=row['tcand'],
		)
		bd.prepare_data(time_window=0.1) # 0.1 seconds
		wfall = np.flipud(bd.sgram.astype(np.float64).copy())
		if bid == 'B6.1':
			plt.imshow(wfall, aspect='auto', interpolation='none')

Here we have added an if-condition to plot the waterfall if it is burst B6 so we can see some output.

.. figure:: /imgs/tut-prep-7.png
	:width: 450
	:align: center

	Burst B6 of Aggarwal et al. (2021) loaded.

To save these bursts we need the bandwidth and resolution information of the file, which is not readily available from the BurstFit package. The Your package is a general library for reading filterbank files (amongst others) and will provide this information.

Then, we will create the ``burstmetadata`` object and save the ``.npz`` for each burst, as in the previous section.

.. code-block:: python
	:lineno-start: 16

	import your
	yourdata = your.Your(files[1]) # assuming both files have same res and band
	datameta = yourdata.your_header

	for bid, row in canddf.iterrows():
		bd = BurstData(
			fp=row['file'],
			dm=row['dm'],
			width=row['width'],
			snr=row['snr'],
			tcand=row['tcand'],
		)
		bd.prepare_data(time_window=0.1) # 0.1 seconds
		wfall = np.flipud(bd.sgram.astype(np.float64).copy())
		if bid == 'B6.1':
			plt.imshow(wfall, aspect='auto', interpolation='none')

		burstmetadata = {
			'dt'        : datameta.tsamp,
			'dfs'       : np.linspace(datameta.fch1-abs(datameta.bw), datameta.fch1, num=datameta.nchans),
			'DM'        : row['dm'],
			'bandwidth' : abs(datameta.bw),
			'duration'  : 0.1,
			'center_f'  : datameta.center_freq,
			'freq_unit' : 'MHz',
			'time_unit' : 's',
			'int_unit'  : 'arb',
			'telescope' : 'Arecibo',
			'burstSN'   : row['snr'],
			'raw_shape' : wfall.shape
		}

		wfallout = f'{bid}.npz'
		np.savez(wfallout, wfall=wfall, **burstmetadata)


.. note::
	You may, as in the previous section, compute the axes and plot the figures to ensure that they match with the publication. Another option however is to load the burst into FRBGui and check the axes from there, as FRBGui will compute them for you from the information provided in the .npz file.

You have now completed the "Preparing Data" tutorials. These have hopefully given you an impression of the kind of tasks that may be involved in loading data from different research groups using different telescopes and different formats. The techniques shown here are certainly not exhaustive, but will hopefully have given you some experience in dealing with the data formats that FRBs can be found in.

Complete Code
^^^^^^^^^^^^^

For your reference, below is the complete script we developed for reading the bursts from Aggarwal et al. (2021) in filterbank format and preparing it as a numpy zipped ``.npz`` file, ready for further analysis in FRBGui or other Python scripts.

.. code-block:: python
	:linenos:

	import pandas as pd
	import numpy as np
	import matplotlib.pyplot as plt
	from burstfit.data import BurstData
	from burstfit.utils.plotter import plot_me
	import your

	canddf = pd.read_csv('all_bursts_bary.csv')
	files = [
		'puppi_57644_C0531+33_0021_subs_0001.fil',
		'puppi_57645_C0531+33_0029_subs_0001.fil'
	]

	canddf['tstart'] = [float(candid.split('_')[2]) for candid in canddf['cand_id']]
	canddf['tcand']  = [float(candid.split('_')[4]) for candid in canddf['cand_id']]
	canddf['dm']     = [float(candid.split('_')[6]) for candid in canddf['cand_id']]
	canddf['snr']    = [float(candid.split('_')[8]) for candid in canddf['cand_id']]
	canddf['label']  = canddf['bidx']
	canddf['width']  = 256
	canddf['file']   = [files[0] if '57644' in str(tstart) else files[1] for tstart in canddf['tstart']]
	canddf = canddf.set_index('label')

	yourdata = your.Your(files[1]) # assuming both files have same res and band
	datameta = yourdata.your_header
	for bid, row in canddf.iterrows():
		bd = BurstData(
			fp=row['file'],
			dm=row['dm'],
			width=row['width'],
			snr=row['snr'],
			tcand=row['tcand'],
		)
		bd.prepare_data(time_window=0.1)
		wfall = np.flipud(bd.sgram.astype(np.float64).copy())
		if bid == 'B6.1':
			plt.imshow(wfall, aspect='auto', interpolation='none')

		burstmetadata = {
			'dt'        : datameta.tsamp,
			'dfs'       : np.linspace(datameta.fch1-abs(datameta.bw), datameta.fch1, num=datameta.nchans),
			'DM'        : row['dm'],
			'bandwidth' : abs(datameta.bw),
			'duration'  : 0.1,
			'center_f'  : datameta.center_freq,
			'freq_unit' : 'MHz',
			'time_unit' : 's',
			'int_unit'  : 'arb',
			'telescope' : 'Arecibo',
			'burstSN'   : row['snr'],
			'raw_shape' : wfall.shape
		}

		wfallout = f'{bid}.npz'
		np.savez(wfallout, wfall=wfall, **burstmetadata)

