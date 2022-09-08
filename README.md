# Frbgui

Frbgui is a graphical user interface for measuring spectro-temporal properties of Fast Radio Bursts from their waterfalls using 2D autocorrelations. It can be used to split bursts with multiple components, change the dispersion measure (DM), add noise filters, and other preparation tasks before measurement. 

Measurements can be performed over a range of DMs and include

* Burst frequency
* Sub-burst slope and/or drift rate
* Duration
* Bandwidth

![Screenshot of Frbgui](imgs/screen1.JPG)

The GUI is extensible and pull requests are welcome.

## Status

Currently Frbgui works best with waterfalls saved as 2d numpy arrays. Frbgui works but is quirky, buggy, and not tested on different platforms. Frbgui will run on any platform but with varying performance. 

## Installation

Install Frbgui with 

```pip install --user frbgui```

## Usage

<!-- You can run the gui from the command line with the command `frbgui` and it will start in the current directory -->

```python
from frbgui import frbgui

frbgui() # starts the GUI
```

## Dependencies

* Light dependency on `your`
