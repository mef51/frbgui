.. _tutarrtimes:


.. show
	* 1 component burst (just setup script and run) (B006)
	* multiple component bursts (xos, and clicking in the window) (11A)
		* read console
	* blended bursts (cuts) (right-click)


Arrival Times Pipeline
============================================

This tutorial introduces using the arrival times pipeline for performing spectro-temporal measurements
of FRBs from their dynamic spectra, or waterfalls. The outputs of the pipeline are essentially the same as with the GUI, however some definitions
are different and the methodology is also different. For those reasons, the arrival times pipeline is
an independent submodule of FRBGUI and can be used independently of the GUI. This tutorial can also be followed
independently of the ':ref:`tutmeasure`' tutorial.

.. For this tutorial we will again use the burst files from the previous tutorials, namely Burst B006 and Burst 11A.

Introduction
------------

The arrival times pipeline works differently than the autocorrelation analysis but the goals are the same: measure the sub-burst slope, multi-component drift rate, duration, bandwidth, and other spectro-temporal properties in a reproducible and convenient way.

The arrival times pipeline fits a Gaussian model in each frequency channel of a waterfall, and tags the beginning of the Gaussian profile (i.e., the 'arrival time'). The signal to noise in each channel is computed, and arrival times in channels with less than 3σ of signal are filtered out. To obtain a sub-burst slope measurement, a line is fit to the frequency vs. arrival time data. The duration and bandwidth are found from 1D Gaussian fits to the integrated time series and spectrum, respectively.
The center frequency is defined as the frequency at the peak of the spectrum fit. For bursts with multiple components, the multi-component frequency drift rate is measured by fitting a line to the center frequency vs component time data.

The pipeline aims to gracefully handle FRB waterfalls with multiple burst components. Such waterfalls can vary from having a few well-separated pulses to being a complicated forest of overlapping and blended components with varying durations. Therefore the pipeline has some features and tools for tagging components, cutting them out, increasing SNR, filtering RFI and/or bad channels, and more.

See `Chamma et al. (2025) <https://academic.oup.com/mnras/article/544/2/2537/8293225>`__ for a more detailed description of the pipeline and ":ref:`docarrtimes`" for the full reference documentation.

Setup
-----

If you don't have them already, download the burst files into your working directory:

    * `aggarwal2021_B006.npz <../aggarwal2021_B006.npz>`_
    * `gajjar2018_11A.npz <../gajjar2018_11A.npz>`_

Create a file named ``measurebursts.py`` to hold your measurement script. This script will serve as your 'interface'.

Type the following code into ``measurebursts.py``:

.. code-block:: python
    :linenos:

    import arrivaltimes

    burst_results = arrivaltimes.measureburst('aggarwal2021_B006.npz')
    print(burst_results)

This is the minimal code needed to perform a measurement. The return value of the :py:func:`arrivaltimes.measureburst` function, stored in ``burst_results``, is an array with an array of each burst component's measurements. For example, if a waterfall has only one component, the ``burst_results`` variable will look like ``[[val1, val2, val3,...]]``.

Run the script with ``python measurebursts.py`` and you should see an interactive matplotlib plot window appear:

.. figure:: /imgs/tut-arrtimes-1.png
    :width: 1000

    Interactive plot window showing the burst measurement.

In this plot you will see the burst waterfall in the center panel with the arrival times it found overlaid (the white points), along with the dispersion measure (DM) the measurement was performed at. The top panel shows the integrated time series, the right panel shows the integrated spectrum in a window around the burst, and the bottom panel shows the sub-burst slope (or single-component drift rate) measurement (i.e. :math:`dt/d\nu`) with the line of best fit.

The interactive plot allows you to zoom and pan the burst to see it more clearly, and also lets you right and left-click and double click parts of the plot for more information about the waterfall. We will use these interactive features later in the tutorial when measuring a more complicated waterfall.

If you close the window by pressing Q or clicking the x you will see some additional output in the terminal::

    user@local:~ $ python measurebursts.py
    Info: Measuring burst aggarwal2021_B006 with 1 components...
    Info: res_freq = 12.500 MHz res_time_ms = 0.08197 ms
    Info: min(data['dfs']) = 974.53 -- max(data['dfs']) = 1774.53 MHz
    Info: aggarwal2021_B006: data['wfall'].shape = (64, 1220), wfall.shape = (64, 1220).
    Info: Using bw_filter = 'data_cutoff' and snr_cutoff = 3
    dtdnu = -4.66369e-04 +/- 3.46984e-04 ms/MHz dtdnu_redchisq = 2.2
    bwidth = 139.013 +/- 3.693 MHz
    Info: Saved aggarwal2021_B006_DM562.056.png.
    burst_results = [['aggarwal2021_B006', array(562.056), 50.306985895993556, 0.01235944557431291, 1393.8700197297696, 3.6101046254094764, 1.0159791308329635, 0.012359445738438335, 139.01317181320326, 3.692720433443701, -0.0004663692105864059, 0.0003469839079227972, 1.1353710822516152, 1.9324598743707861, 2.2040519056015, -0.5951950200809106, 0.4806610536835352, 22]]

.. note::
    Some warnings have been removed from the output for brevity. These primarily arise due to fitting issues in channels with low/no signal and can often be ignored or turned off.

The output includes some information about the measurement, including the time and frequency resolutions of the waterfall, the spectral extent, the dimensions of the waterfall, and some measurement options. It also saved the figure in an image named ``aggarwal2021_B006_DM562.056.png`` and printed the ``burst_results`` array, which contains all the measurements in a row-like format.

To see the names of each value we can print ``arrivaltimes.results_columns``, which we will use later to save a spreadsheet of our measurements::

    >>> print(arrivaltimes.results_columns)
    ['name', 'DM', 't0 (ms)', 't0_err', 'center_f (MHz)', 'center_f_err', 'duration (ms)', 'duration_err', 'bandwidth (MHz)', 'bandwidth_err', 'dtdnu (ms/MHz)', 'dtdnu_err', 'redchisq_t', 'redchisq_nu', 'redchisq_dtdnu', 'tb (ms)', 'tb_err', 'num_arrtimes']

Based on these outputs, we can see the measurement was performed at a dispersion measure (DM) of 562.056 pc/cm\ :sup:`3`, found a center frequency of 1393 MHz, and so on.

Saving measurements
-------------------

It is often important so save measurements in a convenient format for figure making, future analysis, or data sharing.

To demonstrate how to do this, we will imagine that we need to measure this burst at different DMs and then save the set of measurements in a ``.csv`` file. Since the DM is hard to know for certain given that it depends on the electron and ion density along the line of sight, repeating measurements at different DMs is a good way to investigate the covariance of spectro-temporal properties with the DM.

In your ``measurebursts.py`` file, type/update your code to the following:

.. code-block:: python
    :linenos:

    import arrivaltimes
    import pandas as pd

    DMs = [561.0, 561.5, 562.056, 562.5, 563.0]
    results = [] # for storing measurement results
    for dm in DMs:
        burst_results = arrivaltimes.measureburst(
            'aggarwal2021_B006.npz',
            show=True, # Enable the interactive plot window
            save=True, # Save an image file of the plot
            targetDM=dm
        )
        for row in burst_results:
            results.append(row)

    resultsdf = pd.DataFrame(
        data=results,
        columns=arrivaltimes.results_columns
    ).set_index('name')

    fileout = f"results.csv"
    resultsdf.to_csv(fileout)
    print(f"Saved {fileout}")

Here we have added a list of DMs we want to perform our measurements over and placed the :py:func:`arrivaltimes.measureburst` call inside a for loop to repeat the measurements.
We have also explicitly added the ``show=True`` and ``save=True`` keyword arguments (often called ``kwargs``, or options) to the function call in order to control whether or not we want the interactive plot to appear and the ``.png`` file of the plot to be saved.

To save the measurement results in a spreadsheet, we add each array in ``burst_results`` to a list of results (in the ``results`` variable) and then use the Pandas library to save a ``.csv`` file.
Note that we tell pandas to name the columns of the spreadsheet using the ``arrivaltimes.results_columns`` variable.

Run the script again (``python measurebursts.py``) and verify as the plot windows popup that the measurements look reasonable, pressing 'Q' to quickly move from one plot to the next. At the end we should see the line ``Saved results.csv``. Indeed, we should now see a ``results.csv`` file that we can open in any spreadsheet program.

.. figure:: /imgs/tut-arrtimes-2.png
    :width: 1000

    Example output spreadsheet.

We can also see there is an image saved for each measurement at each DM::

    user@local:~ $ ls -1 *.png
    aggarwal2021_B006_DM561.000.png
    aggarwal2021_B006_DM561.500.png
    aggarwal2021_B006_DM562.056.png
    aggarwal2021_B006_DM562.500.png
    aggarwal2021_B006_DM563.000.png


.. figure:: /imgs/tut-arrtimes-3-grid.png
    :width: 1000

    The outputted measurement plots.

Notice that the size of the uncertainty on :math:`dt/d\nu` increases as the burst approaches vertical (and therefore as a function of DM).

Multi-Component Bursts
----------------------

So far we have looked at a very bright, single component burst with simple burst structure.
FRB repeaters can have much more complicated bursts, and a measurement pipeline needs to be able to handle this if we want to efficiently take measurements.

We will adapt our script to add the second burst for this tutorial without deleting the first burst from it.
To do this we will replace the for-loop over the DMs to be a for-loop over filenames (and put aside the analysis over different DMs for the rest of this tutorial).

Type the following into your ``measurebursts.py`` file:

.. code-block:: python
    :linenos:

    import arrivaltimes
    import pandas as pd
    # other imports ...

    files = {
        'aggarwal2021_B006.npz': {
            'show': False, # We are done with the measurement so we don't need to look at it anymore
            'save': True   # Still, save the measurement plot as an image
        },
        'gajjar2018_11A.npz': {
            'show': True,  # We will use the interactive window to improve the measurement
            'save': False  # We are still working on this so no need to save the measurement plot
        },
        # Other bursts can be added here...
    }

    results = []
    for filename, kwargs in files.items():
        burst_results = arrivaltimes.measureburst(
            filename,
            **kwargs
        )
        for row in burst_results:
            results.append(row)

    resultsdf = pd.DataFrame(
        data=results,
        columns=arrivaltimes.results_columns
    ).set_index('name')

    fileout = f"results.csv"
    resultsdf.to_csv(fileout)
    print(f"Saved {fileout}")

In addition to the restructuring, notice that we are now listing the burst files using the ``files`` dictionary, where the key is the filename and the value is another dictionary holding the keyword arguments that are passed into the :py:func:`arrivaltimes.measureburst` function. We have also moved the ``show`` and ``save`` options into each file's value dictionary. In this way we can choose how each waterfall's measurements happens and whether or not we want to see or save the outputs on a burst-by-burst basis.

Running this, we see:

.. figure:: /imgs/tut-arrtimes-4.png
    :width: 1000

    Burst 11A with no good measurements.

There are several issues with the measurement shown.

#. There is a lot of data with no signal, and the signal is hard to see in the middle and narrow
#. The signal-to-noise is low, and the burst appears faint
#. There are no tagged arrival times, and therefore no drift rate measurement
#. There are four resolved components, but the integrated time series only fits the group of them (and fits poorly). Based on what we see in the waterfall, we want four sub-burst slope measurements and one drift rate measurement spanning the frequency drift across the four components.

Let us first make the signal more visible. We can use the zoom tool in the interactive plot window but we would have to do this every time we run the script, and it would not help the measurement. Instead, we'll use the ``crop`` option of :py:func:`arrivaltimes.measureburst`, which allows us to consistently resize the waterfall before measurement.

.. note::
    Click :py:func:`arrivaltimes.measureburst` to see detailed documentation and usage information on each of the options available for the ``measureburst`` function.

The ``crop`` option takes the time indices of the underlying 2D array that represents the waterfall and crops the waterfall to them. We could manually try numbers until the waterfall looks right. Instead, the pipeline allows you to right-click on the waterfall at the left and right points of the window we want to crop around, and the pipeline will print the indices we need to the terminal, as shown in the following animation:

.. figure:: /imgs/tut-arrtimes-5.gif
    :width: 1000

    Right-click on the waterfall to print the frequncy and time index to the terminal, at the position clicked.

The result of right-clicking on the waterfall at good spots before and after the signal prints the following information to the terminal::

    freq chan: 1055 time chan: 803 amp: 0.0258
    freq chan: 1010 time chan: 1305 amp: -0.0141

From this we see that ``(800, 1300)`` would be good values to give to the crop option.
To use the crop option, we just need to add ``'crop': (800, 1300)`` to the options block for the 11A burst.

First however, we also want to increase the signal to noise. To increase the signal to noise, we can use the ``downfactors`` option, which accepts a tuple of integers ``(N, M)`` and will downsample the frequency axis by a factor of ``N`` and the time axis by a factor of ``M``, where ``N`` and ``M`` evenly divide the number of frequency channels and time channels available. From the output in the terminal, we know the watefall has dimensions ``(1690, 2048)``, so ``(5,2)`` is a decent option. We'll add ``'downfactors': (5, 2)`` to the options block and now we should have the following in our script:

.. code-block:: python
    :lineno-start: 5

    files = {
        # ...
        'gajjar2018_11A.npz': {
            'show': True,  # We will use the interactive window to improve the measurement
            'save': False,  # We are still working on this so no need to save the measurement plot
            'crop': (800, 1300),
            'downfactors': (5, 2),
        },
    }
    # ...

.. note::
    In general, we don't want to downsample by too large of a factor, as we will start to lose information from the data and over-influence the measurement. You will need to find a balance between increasing signal-to-noise and preserving the information in the data.

Once we re-run the script we should see a much better view of the burst:

.. figure:: /imgs/tut-arrtimes-6.png
    :width: 1000

    After cropping and downsampling.

Things have improved, and we can see the pipeline is starting to tag some arrival times, but it still doesn't know about the individual components. We can see by the blue vertical dashed line and the gaussian profile in the time-series plot that it has tried to fit the first two components and completely missed the third and fourth components.

To fix this, we will again use the interactive features of the plot window to tell the pipeline where the components are.

We do this by **double left-clicking** on (or close to) the peaks of each component. We will also **double-right click** in between each component to 'cut' each component out from its neighboring component to mitigate the effects of the signal of one component bleeding into it's neighbouring component. Each double-left click adds a dashed blue line while each double-right click adds a black dash-dot line. The information we then need to pass to the pipeline will be printed to the terminal.

.. figure:: /imgs/tut-arrtimes-7.gif
    :width: 1000

    Setting component positions and cuts.

.. note::
    You can also double left and right click on the waterfall plot if the peak of the burst is clearer there.

In the terminal, we will see the following information has been printed after each double-click::

    'xos' : [2.16, 2.03]
    'xos' : [2.16, 2.03, 2.3]
    'xos' : [2.16, 2.03, 2.3, 3.16]
    'xos' : [2.16, 2.03, 2.3, 3.16, 3.8]
    'cuts' : [2.16]
    'cuts' : [2.16, 2.73]
    'cuts' : [2.16, 2.73, 3.41]

Each value is the time in ms where the click happened. Note that there are three values in the ``cuts`` array because there are three black dash-dotted lines, but *five* values in the ``xos`` array because the first value corresponds to the erroneous blue line that was initially there. So, we remove the first value of ``xos`` and paste the last ``xos`` and ``cuts`` lines as measurement options for this burst:

.. code-block:: python
    :lineno-start: 5

    files = {
        # ...
        'gajjar2018_11A.npz': {
            'show': True,  # We will use the interactive window to improve the measurement
            'save': False,  # We are still working on this so no need to save the measurement plot
            'crop': (800, 1300),
            'downfactors': (5, 2),
            'xos' : [2.03, 2.29, 3.13, 3.76], # Don't forget to take out the first value when pasting from the terminal!
            'cuts' : [2.17, 2.78, 3.43],
        },
    }
    # ...

When we run the pipeline again, ``xos`` will be used as initial guesses for the component positions when fitting the time series, which will now be a sum of 4 Gaussian profiles, corresponding to the number of values in the ``xos`` array. The measurement now has 4 reasonable looking sub-burst slopes, shown with the white lines, and a drift rate, shown with the red line.

.. figure:: /imgs/tut-arrtimes-8.png
    :width: 1000

    Setting component positions and cuts.

The measument looks quite good now and it might be tempting to call it complete.
However, there is still a big problem.

We can see that the burst bandwidths, shown with the gaussian profiles and the blue shaded regions on the spectrum (right-most) plot, are much narrower than the bursts. The first two components are the main problem; the profiles have fit to the bright core of the burst at around 7000-7200 MHz but missed the parts of the signal that extends up to 7500 MHz.

This is a common issue with FRBs, as their measured intensities can vary significantly as a function of frequency. These intensity variations can be caused by multiple sources, including scintillation, the (unknown) emission mechanism, and the frequency response of the observing telescope.

In this case, we can see that if the brightest scintles, or intensity peaks, for the first two components are masked out, the fit may be more realistic. Using the intensity axis on the spectrum, the spikes seem to begin at an intensity value of around ``0.11`` (using the interactive plot to get precise values).

If we now add the ``'bandmask_thres': 0.11`` option to our measurement script, we will get the following measurement:

.. figure:: /imgs/tut-arrtimes-9.png
    :width: 1000

    Using the ``bandmask_thres`` option to mask scintles improves bandwidth measurements.

A hatched region has appeared on the plot indicating where scintles are masked out and little red 'x's denote the masked peaks. We can also see the bandwidth measurement is more accurate, and the blue shaded region representing the bandwidth is now much broader than before.

At this point we can call the measurement complete, and should have all the information we need to review the measurement in the output plot.

Summary and Template Script
---------------------------

If you've made it this far, congratulations! You now have a good idea of how to use the pipeline and how to take, review, and modify measurements of FRB waterfalls with single or multiple components.

At this stage its a good idea to read through all the measurement options that :py:func:`arrivaltimes.measureburst` allows, as they cover many more features that are not used in this tutorial.

Some of these additional features include:

* Masking
* Masking components individually
* Changing the signal-to-noise cutoff
* Tweaking what is shown in the output plot, the colormap, the figure size, location of the legend, etc.
* Setting initial guesses for burst widths
* Forcing the pipeline to fix your initial guesses and not move them
* Different types of filters to apply to the arrival times
* Manually setting the noise sampling window
* Saving and loading time series solutions (very useful for waterfalls that take a long time)
* Change where the drift is measured from (beginning, middle, or end of the profile)
* Display diagnostic plots by frequency channel
* And others...

For your reference, below is the code we developed during this tutorial, modified slightly so that it can be used as a template script for larger more complicated measurement scripts with hundreds of bursts.

.. code-block:: python
    :linenos:

    import arrivaltimes
    from collections import namedtuple
    import pandas as pd
    # other imports ...

    Dataset = namedtuple('Dataset', 'run save path')
    datasets = {
        'tutorial' : Dataset(True, False, ''),
        # Other datasets
    }

    if datasets['tutorial'].run:
        prefix = datasets['tutorial'].path
        files = {
            'aggarwal2021_B006.npz': {},
            'gajjar2018_11A.npz': { # scintillation/scalloping, note use of `bandmask_thres` option
                'show': True,  # We will use the interactive window to improve the measurement
                'save': False,  # We are still working on this so no need to save the measurement plot
                'crop': (800, 1300),
                'downfactors': (5,4),
                'xos' : [2.03, 2.29, 3.13, 3.76],
                'cuts' : [2.17, 2.78, 3.43],
                'bandmask_thres': 0.11,
            },
            # Other bursts...
        }

        save = datasets['tutorial'].save
        results = []
        for filename, kwargs in files.items():
            filename = f'{prefix}{filename}'
            burst_results = arrivaltimes.measureburst(
                filename,
                outdir='',
                show=not save,
                save=save,
                **kwargs
            )
            for row in burst_results:
                results.append(row)

        resultsdf = pd.DataFrame(data=results, columns=arrivaltimes.results_columns).set_index('name')
        fileout = f"results.csv"
        if save:
            resultsdf.to_csv(fileout)
            print(f"Saved {fileout}")

Here, we introduce a named tuple called ``Dataset`` that gives an additional organizational tool to the measurement script. Using it, we can control whether a dataset is run, saved, and specify the path to the data. For example, you could organize the script around different FRB repeating sources, and use the ``Dataset`` tuple and ``datasets`` dictionary to manage your data.
