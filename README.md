# Frbgui

Frbgui is a graphical user interface for measuring spectro-temporal properties of Fast Radio Bursts from their waterfalls using 2D autocorrelations. It can be used to split bursts with multiple components, change the dispersion measure (DM), add noise filters, and other preparation tasks before measurement.

Measurements can be performed over a range of DMs and include:

* Burst frequency
* Sub-burst slope and/or drift rate
* Duration
* Bandwidth

![Screenshot of Frbgui](imgs/screen1.JPG)

The GUI is extensible and pull requests are welcome.

## Status

Currently Frbgui works best with waterfalls saved as 2d numpy arrays (See [Usage](#usage) for more details). Frbgui works but is quirky, buggy, and not tested on different platforms. Frbgui will run on any platform but with varying performance.

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
    'dt'        : # array of time axis, unused, optional,
    'dfs'       : # array of frequency channels in MHz,
    'DM'        : # float of dispersion measure (DM),
    'bandwidth' : # float of bandwidth in MHz,
    'duration'  : # burst duration in seconds
    'center_f'  : # burst frequency in MHz, unused, optional,
    'freq_unit' : # string of freqeuncy unit, e.g. 'MHz', optional,
    'time_unit' : # string of time unit, e.g. 'ms', optional,
    'int_unit'  : # string of intensity unit, e.g. 'Jy', optional,
    'telescope' : # string of observing telescope, e.g. 'Arecibo', optional,
    'burstSN'   : # float of signal to noise ratio, optional,
    'tbin'      : # float of time resolution, unused, optional
}
np.savez('burst.npz', wfall=wfall, **burstmetadata)
```
