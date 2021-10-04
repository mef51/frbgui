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
* revising DM range for frb180301 to 515.2 - 525.5. Low end based on burst 5 DM error
* add ability to specify initial fit guess and redo DMs

## july 25-26
* add burst splitting for frbs with multiple sub-bursts
* redo dm and p0 for sub-bursts is buggy

## aug 14
* first pass at getting measurements, fixing some bugs

## aug 16
* add burst region importing and exporting

## aug 17
* revising dm range for frb180301 to exlude extremes when burst is very low snr: 516.3 - 522.7 (521.1)

## aug 18
* frb180813:
	* snr too low: bursts 2, 3, 4, 6, 14, 15
	* signal likely cutoff: 8
	* worth retrying: 9c, 10b
	* burst 12 is possibly fit to the wrong feature
* new bugs:
	* result table for multiburst doesnt show subburst results after navigating to different burst
	* cant disable p0 interface group on single burst after performing multiburst (sometimes.. do burst 9 and then burst 1 to reproduce)

## aug 21
* results look pretty scattered for this source but its points vaguely follow the other sources... i want to try getting better statistics on frb121102 instead of adding more sources i think..
* for this source drive the higher end of the dm range even lower and maybe do more trials. exclusion logger shows most bursts are positive at 518pc/cm3
* did measurements from DM 516.3-520
* let's try 516-518.5

## aug 23
* consider telescope selection effects. maybe all these bursts agree cuz telescopes with similar capabilities are performing these observations. This is why it might be more useful to focus on FRB121102, since it is observed in many time resolutions and at different frequencies

## aug 24
* I use a subset of the bursts from Luo et al based on SNR. the dm range implied by that subset is 515.9-518.3. With the error bars, thats 513.6-519.5
* TODO:
	* fix axes
	* fix masking when displaying true axes
	* fix resolution update after subsampling

## aug 26
* might need fancier noise removal than point-and-click. transient short spikes are hard to spot, even when tweaking the colorbar. Can't seem to reproduce some figure plots:
	* R1_B04
	* R1_B05
	* R1_B11
	* R1_B12
* HMMM looks like some of them are flipped in the time-axis??
	* R1_B13

## aug 27
* sksg filter failed to reveal the missing bursts, even when applied on the whole waterfall

561-569
570
572
* I will start with a dm range of 560-570 with a step of 0.5
	* redo B14.. there's a bug when redoing an entire burst
* fits were found really fast for this data, maybe somethign to do with the intensity axis?
* either way I have data for about 25 new points. I feel kind of lost so I want to write some diagnostic plots that take in the giant csvs and make a document that can be reviewed. I also want to automate the slope vs duration figure a little better.
	>> plotResults from results csv
	>> automate best DM selection
	>> altair version of figure 1? would be nice to click a point and see its underlying fits. but can do this with annotations and plotResults to start
* instead of "determining a good DM range" it feels more like taking a huge DM range and then "finding" the DM range where most of the bursts have valid measurements. Maybe the range found this way means something?

## aug 29
* added plotResults to make pdf of every fit

## sept 01
* drift range becoming negative is common in oostrum dataset. seems to be related to slope errors > 40%
	* could add filter based on proportion of error or could just explicitly exclude any measurement that flips the dm range sign
* there are several horizontal tracking points as well, these are problematic because they indicate that the DM choice can cause the trend

## sept 6
* working on reading aggarwal et al. 2021 data
	* parse cand_id from the csv and use BurstData to load the waterfall from the filterbank

## sept 8
* subsample_cb and applyMasks and subtractbg can probably be refactored into a single 'cleanWaterfall' function
* check resolution of aggarwal.. subsamplecb affects the resolutions differently for this dataset

## sept 11
* batch measure feature: add bursts to a list and measure all with one button
* its very easy to discard low snr bursts. this is a bias
* splitting is broken for burst 121 aggarwal and splitting is not implemented for pdf export
	* add more metadata to results csv (ie subtract bg, region split, subsampling, etc. Everything that is needed to reproduce)
* num results is not updated after loading results
* export result is broken for cropped waterfalls

## sept 13
* regions, subtract bg, and twidth should be stored by burst but are currently global. storing by burst will let me export the measurement parameters by burst

## sept 16
* done measurement param export/reload, fix multiburst region splitting next

## sept 18
* fix multiburst. was reading the wrong time extent when splitting up subbursts

## sept 19
* check time downsampling
* B025 at DM 563.6315, 564.2, 564.7 has the wrong shape
* trial DM range is wrong with dm step = 0.5

## sept 20
* B025 problem seems to be that pkidx is finding a different peak as the dm changes. and it does kinda seem like there's another burst in there. but the real issue is that's its low snr, and i should just downsample in time to ensure i get the right burst
* bug: downsampling in time changes length of burst
	* redo twidth handling to respect displayed duration instead of displayed # number of channels

## sept 21
* added a backup feature
* check twidth loading
* check nans in output csv
* paginate pdf files so that you dont have a humongo pdf
* oops i never implemented regions in plotResults

## sept 27
* turns out GBT has REALLY good spectral resolution. When I saw that the data shape was (2048, 19456), I assumed it was 2048 freq channels and 19456 time channels. But checking the length of the dfs array in the fits file matches 19456 frequency channels. This corresponds to about 183khz resolution, or half the 366Khz resolution mentioned in the paper

## sept 28
* why does `driftrate._dedisperse` do nothing? why do the gajjar data look overdedispersed/unanchored?
* turns out the gajjar data was dedispersed already, so passing `prepare=True` to pypulse was messing it up. Also data is upside down, just flip it.

## sept 30
* todo:add masking by range for 11B, 11G, 11M, 12A
* load 11J with 256x256 subsampling to trigger autocorr2d stack trace

## oct 1
* fixing bugs revealed by first measurement pass
	* remove raw_shape attribute in npz file (ie. the stored shape), just infer it from the stored waterfall
* there's a problem with region saving
	> load a burst, set regions, change burst, go back, export regions, saved key will be wrong
* with gajjar set fits in pdf dont line up with fit found in gui

## oct 2
* regions dont get loaded at gui start if first burst in list has regions
* regions dont make sense without the twidth they were made with
* rounding from time to channel number is needed to display the waterfall, but using this extent to update the time resolution was leaving a tiny error in the resolution that was causing a PDF display bug. Using the full shape and duration to update resolution always now
	> just as a general rule respect the value of duration and time in the fits file. I was inferring them from dfs but since that can be a list of the center frequencies of the channel you could have slight errors in the bandwidth/duration that affects the downsampling
