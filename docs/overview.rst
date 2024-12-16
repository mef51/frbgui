.. _overview:

FRBGui
======

.. preview with make html && open _build/html/index.html or use livereloadx -s docs/_build/html
.. rm -rf _build/* to regenerate the index

FRBGui is a graphical user interface for measuring spectro-temporal properties of `Fast Radio Bursts`_ (FRBs) from their waterfalls using 2D autocorrelation functions (ACF). It can be used to split bursts with multiple components, change the dispersion measure (DM), add noise filters, and other preparation tasks before measurement. FRBGui is built with the `Dear PyGui <https://github.com/hoffstadt/DearPyGui>`_ interface library.

After measurement, FRBGui can be used to review and fix incorrect measurements by supplying initial guessues, producing an output PDF of all measurements and spreadsheets of measured values. Spreadsheets produced in this way can be loaded back into the FRBGui for further work or different sessions of measurements.

Measurements can be performed over a range of DMs and include the burst frequency, the sub-burst slope and/or drift rate, burst duration and burst bandwidth.

In this documentation you will find instructions on installing and getting started with FRBGui, an overview of its features, a guide on using FRBGui's scripting capabilities along with an API reference, tutorials for preparing data and performing burst measurements, and a reference table of FRBGui's measurement outputs. For advanced users, a section will introduce how custom tools, windows, and other interface elements can be added using Dear PyGui.

These pages are still being added, and if there is a specific topic you would like to see that is not yet available, please open an issue on `GitHub <https://github.com/mef51/frbgui>`_ or email me directly at chammam at mcmaster . ca

.. _Fast Radio Bursts: https://en.wikipedia.org/wiki/Fast_radio_burst

Features
--------

* Change the DM of a burst
* Crop waterfalls
* Mitigate noise and RFI via background subtraction, SK-SG filter, manual channel zapping, and mask ranges
* Import and export of noise masks
* Measure over user-defined ranges of DMs
* Downsample waterfalls in frequency and time
* Split bursts with arbitrary numbers of sub-burst components
* Define inital fit guesses
* Review measurements via the output table
* Correct individual fits by DM
* Output measurements as a CSV spreadsheet and/or PDF with plots of each waterfall with its measurements.
* Provides :py:mod:`driftrate` and :py:mod:`driftlaw` python modules for scripting and automation.
* Automatic backups of measurements as they are made

.. image:: ../imgs/screen1.JPG

.. _installation:

Installation
------------

Install FRBGui with::

	pip install --user frbgui

For a local, editable installation with bleeding edge updates you may clone the repo and install locally::

	git clone https://github.com/mef51/frbgui.git
	cd frbgui
	pip install --user --editable .


Usage
------

.. image:: ../imgs/demo.gif

.. Getting Started
.. ---------------

Run from the command-line with the following command to start in your current working directory::

	frbgui

In a python script, you can invoke the gui in the following way:

.. code-block:: python

	from frbgui import frbgui
	frbgui() # starts the GUI

.. _burstformat:

Burst Format
------------

FRBGui works with burst waterfalls that are prepared as python ``.npz`` archives.

The following snippet shows the format of the archive and an example of how a burst can be saved in the right format:

.. code-block:: python
	:linenos:

	import numpy

	wfall = # 2D numpy array with shape (num freq channels, num time channels)
	burstmetadata = {
	   ### required fields:
	   'dfs'       : # 1D array of frequency channels in MHz, lowest first
	   'DM'        : # dispersion measure (DM) in pc/cm^3, float
	   'bandwidth' : # bandwidth of `wfall` in MHz, float (negative is ignored)
	   'duration'  : # duration of `wfall` in seconds, float
	   ### optional fields:
	   'center_f'  : # burst frequency in MHz, optional,
	   'freq_unit' : # string of freqeuncy unit, e.g. 'MHz', optional,
	   'time_unit' : # string of time unit, e.g. 'ms', optional,
	   'int_unit'  : # string of intensity unit, e.g. 'Jy', optional,
	   'telescope' : # string of observing telescope, e.g. 'Arecibo', optional,
	   'burstSN'   : # float of signal to noise ratio, optional,
	   'tbin'      : # float of time resolution, optional
	}

	np.savez('burst.npz', wfall=wfall, **burstmetadata)

Optional fields are used for display purposes and do not otherwise affect measurements from within Frbgui.

Acknowledgements
----------------

If used in an academic study please cite:

* *A broad survey of spectro-temporal properties from FRB 20121102A*, Chamma, Mohammed A.; Rajabi, Fereshteh; Kumar, Aishwarya; Houde, Martin. `MNRAS`_, 522, 2, 3036-3048, June 2023. `arXiv:2210.00106 <https://arxiv.org/pdf/2210.00106.pdf>`_

.. _MNRAS: https://academic.oup.com/mnras/article/522/2/3036/7120059

Publications
------------

In addition to the above paper FRBGui has been used in the following studies:

* *Validating the Sub-Burst Slope Law: A Comprehensive Multi-Source Spectro-Temporal Analysis of Repeating Fast Radio Bursts*, Brown, Katie; Chamma, Mohammed A.; Rajabi, Fereshteh; Kumar, Aishwarya; Rajabi, Hosein; Houde, Martin. (Jan 2024) MNRAS:Letters, Vol. 529 No. 1 p.L152-158. `doi:10.1093/mnrasl/slae012 <https://doi.org/10.1093/mnrasl/slae012>`_

* *Characterization of the repeating FRB 20220912A with the Allen Telescope Array*, Sheikh, Sofia Z. et al. (Jan 2024). MNRAS , Vol. 527 No. 4 p. 10425-10439. `doi:10.1093/mnras/stad3630 <https://doi.org/10.1093/mnras/stad3630>`_

Attributions
------------

`Meteor icon created by Freepik - Flaticon <https://www.flaticon.com/free-icons/meteor>`_

.. Overview

.. Getting Started

.. Tutorial: Preparing data
.. Tutorial: Measurements of an FRB waterfall

.. Scripting and Modules
   driftrate.py
   driftlaw.py

.. Output Table Reference

.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
