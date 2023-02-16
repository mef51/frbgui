# Frbgui

Frbgui is a graphical user interface for measuring spectro-temporal properties of Fast Radio Bursts from their waterfalls using 2D autocorrelations. It can be used to split bursts with multiple components, change the dispersion measure (DM), add noise filters, and other preparation tasks before measurement. 

After measurement, Frbgui can be used to review and fix incorrect measurements by supplying initial guessues, producing an output pdf of all measurements and spreadsheets of measured values. Spreadsheets produced in this way can be loaded back into the Frbgui for further work or different sessions of measurements.

Measurements can be performed over a range of DMs and include:

* Burst frequency
* Sub-burst slope and/or drift rate
* Duration
* Bandwidth

![Screenshot of Frbgui](imgs/screen1.JPG)

The GUI is extensible and pull requests are welcome.

## Status

Currently Frbgui works best with waterfalls saved as 2d numpy arrays (See [Usage](#usage) for more details). Frbgui is functional and can produce thousands of measurements but is quirky, buggy, and not tested on different platforms. Frbgui will run on any platform but with varying performance.

## Installation

Install Frbgui with

```pip install --user frbgui```

## Usage

Run from the command-line with the following command to start in your current working directory:

```bash
frbgui
```

In a python script, you can invoke the gui in the following way:

```python
from frbgui import frbgui

frbgui() # starts the GUI
```

At the moment frbgui works best with burst waterfalls that are prepared as python `.npz` archives. The following snippet shows the format of the archive and an example of how a burst can be saved in the right format:
```python
wfall = # 2D numpy array with shape (num freq channels, num time channels)
burstmetadata = {
    ### required fields:
    'dfs'       : # array of frequency channels in MHz,
    'DM'        : # dispersion measure (DM) in pc/cm^3, float
    'bandwidth' : # bandwidth of `wfall` in MHz, float
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
```

Optional fields are used for display purposes and do not otherwise affect measurements from within `frbgui`.

## Acknowledgements

If used in an academic study please cite

"A broad survey of spectro-temporal properties from FRB20121102A", Chamma, Mohammed A. ; Rajabi, Fereshteh ; Kumar, Aishwarya ; Houde, Martin. Oct. 4 2022. Submitted to MNRAS.
[arxiv:2210.00106](https://arxiv.org/abs/2210.00106),
[ADS:2022arXiv221000106C](https://ui.adsabs.harvard.edu/abs/2022arXiv221000106C/abstract)
