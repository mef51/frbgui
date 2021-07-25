## nov 4

papers:
	* chime https://www.nature.com/articles/s41586-020-2863-y
	* bochenek et al. https://www.nature.com/articles/s41586-020-2872-x
		* https://github.com/cbochenek/STARE2-analysis
	* franz et al. https://arxiv.org/abs/2007.05101
	* Lin et al. https://www.nature.com/articles/s41586-020-2839-y

## nov 8

* i can read the chime bursts, snr seems a little low but i think its fine
* stare2 data:
	* sigproc throws an error
	* sigpyproc doesnt run anywhere

## nov 9

* i can read the kirsten bursts
* the STARE2 burst is the third burst in the pulse train where the 2 CHIME bursts happen earlier. Very happy trombone
* not sure what dfdt is giving me

## nov 10

* chime detections yield stripey autocorrelations due to the noise so I'm going to use the ARO detection instead
	* ARO snr has less noise but snr is much lower so autocorrelation is bad.. maybe I can remove some noise? The plot in the paper seems clearer
	* component 1 seems to have a thin and fat component
* for meeting:
	* dfdt
	* bursts

## nov 18

* DMs used
	* CHIME Detection  : 332.7206 +/- 0.0009
	* STARE2 Detection : 332.702(8)
	* VLBI Detection   : 332.7206

## nov 19

* loaded stare2 bursts.

## nov 23

* VLBI bursts still have positive drift at STARE2 DM
* Additional detections:
	* http://www.astronomerstelegram.org/?read=14074
	* http://www.astronomerstelegram.org/?read=14080
	* http://www.astronomerstelegram.org/?read=14084
* at DM = 332.5ish (less than 0.1% difference from the quoted DMs) the fit is pretty good.

## nov 24

* rewriting to make fig 1 in terms of 1/drift since we are near an asymptote

## nov 25

* kink in spectra is most likely caused by integer cast in dedispersion
* Savitski-Golay filter for width preserving filtering

## nov 26

* factor of sqrt(2) difference between the tau_w error equations

## dec 3

* redshift_host and frb_name (internal_name on TNS) should be ok "is this a repeater?" filters
* https://www.chime-frb.ca/repeaters

## dec 8

* luo et al. 2020 (FAST, FRB180301) only provides pol'n data. sent an email

## dec 9

* got data from Oostrum2020 et al. for frb121102 in psrchive format
* requested data from Marthi2020 et al
* requested data from Cruces2020 et al

## dec 15

* du et al. 2021
	* NS/companion models don't fit with frb121102 and frb180916.j015(..)
	* suggests that other triggers should exist.
		* What could those be? maser flares?
			* Likely unobservable, or difficult to observe, given that the frb outshines everything. need to observe a close by FRB
			* if the rate for frbs is 1/50yrs/galaxy then trying to observe a trigger in our galaxy or a nearby galaxy will be a difficult effort.

## dec 16

* pipeline:
	* [.fil, .fits, .ar] --> .npy --> burst finding/windowing, masking, splitting, simple noise removal --> dm variations, processburst --> parameter normalizing --> save, model values, figures
	* pysigproc3 and `your` can handle .fil, pypulse can handle .ar if you first convert it to fits with psrconv

## dec 18

* omg I just discovered DearPyGui and it is glorious. better than the other python wrappers for Dear Imgui

## dec 21

* dearpygui's heat series can display dynamic spectra from the data directly.

## dec 25

* two columns in one window is unwieldy. Do one window for the plots and another window for analysis
* dpg's async utility is too slow, opened (https://github.com/hoffstadt/DearPyGui/issues/407)

## jan 4 2021

* hypothesis: the closer to 1.4ghz (or whatever nu_0 is) an FRB is the smaller the bandwidth of the burst since the source is less relativistic
* subsample fails on dimensions that are a third of the original? (if not evenly divisible)
	* `len(a.flatten()) == dim1*dim2*dim3`
	* `dim3 == len(a)/(dim1*dim2)` and dim3 must be an int
	* i dunno how to handle this, I'm just gonna handle the error by forcibly incrementing the dimension by one which will handle all the cases that matter since the point is just to increase the snr in the dynamic spectrum

## jan 20

* what happened to the background in 4 and 11?

## may 27
* ligo's allsky plotting thing should be merged into astropy https://lscsoft.docs.ligo.org/ligo.skymap/plot/allsky.html

## may 31
* FRB20200120E (Bhardwaj et al. 2021, Majid et al. 2021)
	* 5 CHIME bursts (400-800 MHz) and few ms long
	* Majid et al. has a pulse train on the microsecond scale, with very steep slopes
	* Nimmo et al. 2021 has 4 nice bursts on the microsecond scale
		* rough calculation of the slop in burst B4 appears to be close to our trend (a little below)
* psrconv makes .jar into .rf?? >> .rf is .fits. pypulse can read it
* what if i made a website and asked people to upload their downsampled and dedispersed waterfalls to it?
	* what if TNS required that?

## june 1
* FRB20201124A >> atel
	* https://www.astronomerstelegram.org/?read=14538
	* http://www.ncra.tifr.res.in/~viswesh/Atel_burstpanorama_timeordered.pdf

## june 3
* waiting for data to load while guessing what downsample factors i want is such a headache
	* also a problem if trying to remove bg or mask bad channels before downsampling
	* just keep a burst in memory while figuring it out
	* subtractbg streaks the luo et al data
	* add burst selection and data loading into gui

## june 6
* add very basic clickmasking via a click callback
	* handle crashes
	* handle view changing
	* list masked channels
	* clicking same channel should remove mask
	* masks should be index by burst, so reset the mask if the burst changes
* this is going nowhere i still can't see some of these bursts
* dpg table data shape is transpose of dictionary shape
	* set_table_data is useless. massaging sparse 2d arrays for it is painful. add columns and their data one by one

## june 9
* burst 3 was just hidden due to low snr in the time axis
* double click to reset view isn't working anymore
	>> happens because subsampling changes bounds but plot still has same bounds
* interface assumes masking happens before subsampling, which is best for retaining as much data as possible, but not necessary.
TODO:
	* find a better subsampling algorithm
	* add a way to remove mask
	* add a way to mask after subsampling
	* add time series
	* add autocorrelation
	* just analyse the luo data manually and then go back to gui
* keeping the view when masking seems to break because of set_plot_xlimits_auto

## june 12
* based on the breadcrumbs found in the paper and the data files they provide, the time resolution is either 49.152us or twice that
	* the number returned from the fits file using .getTimes() is the full duration in the units returned by .getTimeUnit()
		* ie. 0.70464307 seconds over 28672 samples = 24.576 us, which is half the time resolution from the paper. Not sure how it could be half but its one of those two numbers. Pick one and see if it lines up with their figures
* added burst metadata to load, will use savez to store subsampled waterfalls with metadata

## june 14
* lol i think my corr extents are wrong and have been wrong for two papers
	* why did i have a factor of 2??
* not sure i can do contour plots easily in DPG, will plot lines that represent the semi-major and -minor axes
* adding bounds breaks masking
	>> its because i use physical units on the axes instead of channel nums

## june 15
meeting notes:
	* consider removing the one frb121102 fit that only uses the one burst at the high DM, this will fix the range
	* fix the correlation extents in the figures that are affected. Old figures are unnaffected

## july 12
* DM range for FRB180301: 516.3 - 525.5 (there's a lower DM but its on a tiny scrap of FRB so whatever)

## july 14
* played with aggarwal et al.'s burstfit package. might use it for noise primarly and burst modelling. autocorr of model will always be cleaner

## july 15
* 'keepview' feature is very bugged
* can produce results csv over DM range now. need to be able to specify a p0 tho, and still no visual way to evaluate fits over the range (just the fit at the burstdm atm)

## july 18
* dpg.get_value produces rounding errors. for e.g. burst dm is 518.3, but becomes 518.2999877 after dpg.get_value
* implemented buttons for switching DM of displayed burst
* results seem repeated? also left dm select button
	* yikes was passing the wrong waterfall to processDMRange. need to pass the wfall at the burstdm, as processDMRange will dedisperse to each trial DM
* scipy.curve_fit is getting completely stuck?? >> was passing huge waterfall instead of windowed waterfall

## july 19
* ask about redshift

## july 21
* wip on multi-burst results

## july 24
* make image from matplotlib corr2d plot and use that in the gui. this is a workaround that can show the fit contours since dpg cant do contour plots
