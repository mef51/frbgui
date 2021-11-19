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

## oct 4
* I measured 11A at two different wfall resolutions: 256x512 and 512x1024. Most of the slope measurements agreed, but the measurements of slopes that were close to vertical differed by an order of magnitude. Looking at the figures, you can see that dedispersion for the lower res waterfall was steeper than for the higher res waterfall, because the channels are more coarse.
	* Measuring at a time resolution that is too low will lead to a larger slope measurement, and make the range of slope measurements larger over the DM range. Dedispersion becomes more inaccurate when the time resolution is low (fewer channels) which means larger measurements as well as larger errors at the vertical end of the slope range.
	* Use as high a time resolution as possible for higher quality measurements
* need to fix:
	* loading a csv of multiburst results displays only the full burst results and not the subburst results >> fixed
	* pdf of multiburst at high res (maybe low res too?) doesn't split burst properly
		* see 11A_b @ DM=556. low res is fine, higher res doesnt split properly
	* clicking between subburst results in table after measurement leads to:
		```
		Traceback (most recent call last):
			File "B:\dev\frbrepeaters\frbgui.py", line 494, in resulttable_cb
			    displayresult_cb('User', {'resultidx': newsel[0]})
			  File "B:\dev\frbrepeaters\frbgui.py", line 762, in displayresult_cb
			    plotdata_cb(sender, data)
			  File "B:\dev\frbrepeaters\frbgui.py", line 299, in plotdata_cb
			    subburst = subbursts[subname]
			TypeError: list indices must be integers or slices, not str
		```

## oct 5
* slope and duration measurements are resolution dependent. So how should you compare burst measurements across different resolutions?
	* Should burst measurements across sources be compared at the same resolution?
* gajjar resolution in first paper (ascii data from authors) was at 1.46484375 MHz/chan and 4.0958926622802e-5 s/chan. That corresponds to a shape of 2048x512. The data they provide online was 19456x2048 so I can prepare wfalls of shape 2432x2048 and use that to reproduce the old measurements at the same resolution

## oct 6
* try ftol and xtol in `curve_fit` to speed up fit finding https://stackoverflow.com/questions/31070618/how-to-speed-up-python-curve-fit-over-a-2d-array
* differences:
	* different shapes
	* wrong res used in old 11A measurements? measurement seems right by manual calculation
	* dedisperse before subsample (old) vs subsample before dedisperse (new)
* redo with shape  to get resolution 7.32421875mhz/chan and 4.09e-5 s/chan (5x base freq resolution)
* write an automated test that given a results csv will, for each row, repeat the measurement and compare the results with those in the csv.
	* If they differ, raise an error.

## oct 7
* everything looks the same in the plot, but the solutions are different. try plotting each solution on the other, are they equivalent?
	* they are not
	* oh snap if you divide the gui's sigmax, sigmay and angle by 4 you retrieve paper 1's solution. extents bug? time res is off by 4?
* why are the gaussian amplitudes so different?

## oct 8
* tres goes through fine. Turned out I copied 11A instead of 11A_a and it turned out the solution for 11A/4 is identical to 11A_a. Why does that work?
* minor difference seen between the results when plotting with gui's 11A_a and paper 1's 11A_a

## oct 10
* to maximize precision the order of steps for treating a waterfall before measuring should be:
	1. Apply masks --> Do at high res to preserve as much signal as possible
	2. Dedisperse --> do this at high res because the channels can be rolled more finely. Low res dedispersions can result in huge slope ranges
	3. Downsample --> This step is necessary to increase the SNR enough for a measurement (ie. a fit) to be performed.
	4. Center/crop/pad burst --> for speed and display and measurement checking
	5. Measure
* If there are regions, then they should be time stamped (not channel numbers) and ideally interpolated into channel numbers and split before downsampling
* In the jupyter notebook if I subsample before dedispersing in an attempt to mimic the gui, the results still dont match. There's another difference.
	* In fact subsampling (in just frequency at least) before dedispersing only has a small effect on the measurement haha.
* Width differences?
	> nope
* maybe its a regions thing? Compare 11D between the two since its a single burst.

## oct 11
* 11D is also different. reproduce that since its a simpler case and then do 11A. dont gotta deal with the regions stuff that way
* 11D: gui solution and notebook solution both look good on either correlation. So why the difference in slope?
* notebook with gui wfall finds same slope as notebook wfall but different duration measurement
* notebook with notebook wfall with NO p0 finds same slope and duration as notebook with gui wfall with NO p0
* fuck it ill just make a table

11D @ DM = 555
_________________________________ | slope        | t_w_ms
notebook + gui wfall      + no p0 | -1162.344651 | 1.590012
notebook + notebook wfall + no p0 | -1162.344651 | 1.590012
notebook + gui wfall      + p0    | -1162.347946 | 0.354748
notebook + notebook wfall + p0    | -1162.347946 | 0.354748

* the angle is not "cleaned" before calculating t_w_ms. I want to recheck that calculation
* Specifically these two lines:
	```
	theta = popt[5] if abs(popt[3]) > abs(popt[4]) else popt[5] - np.pi/2
	popt[5] = theta
	```

* replacing angle in popt with theta reproduces the duration measurements for 11D. This makes sense as `driftlaw.computeModelDetails` assumes that 'angle' is the correctly defined theta
* that line (popt[5] = theta) was in the michilli script but not the chime script, so those bursts probably have bad time ranges lmao i need to re-run them
	* actually looks like the scripts for both 180916 and 180814 stores the correct angle and not the raw angle, which means their durations should be fine
	* so the gajjar script was the only one with that bug and luckily it just worked out to always use the correct angle.
	* When I wrote the gui and the libs i didn't copy from michilli so the key line was missing, and I handled the storage differently than from the chime scripts, which left the bug exposed
	* btw for the record the sigma_t_ms and t_w_ms columns should basically be the same. they are calculated differently but sigma_t turns out to be a pretty good estimate of t_w_ms
* notebook produces the same slope and duration for the ascii waterfall and the gui waterfall (ie fits)
* maybe dump the ascii waterfall and measure it in the gui?

## oct 12
* 11D: notebook reproduces gui slope and duration result when given the fits data. This points to the differences in the ascii and fits data as the source of the discrepancy. This could also be due to victor's treatment of the ascii data
* there is a definite rotation between the two versions of the data. One of them is maybe at the wrong reported DM?

## oct 13
* The gajjar data in ascii form that we obtained through email is in fact not at 565 pc/cm^3. It is at 557.91 pc/cm^3, matching the DM the fits data is at.
* This makes the measurements for 11D agree. There are still differences with 11A

## oct 14
* width of region affects slope and center_f measurement for 11A_a in notebook
* region also affects slope measurement in gui!
* this can be easily reproduced by the choice of background

## oct 15
* zero padding instead of background sampling padding removes width dependence
* burst region width only seems to matter if region is not generously wide. For 11Aa there isn't much room, which seems to be the reason for the measurement instability
* 11A Regions:
regionname = 'Region1' regiontype = 1 region = [53, 70]
regionname = 'Region2' regiontype = 1 region = [70, 83]
regionname = 'Region3' regiontype = 1 region = [83, 101]
regionname = 'Region4' regiontype = 1 region = [101, 144]
regionname = 'Region5' regiontype = 0 region = [0, 41]

* measurements are much closer together now when regions are duplicated. Differences are due to differences in signal in the window, so this means we should be careful when making regions, and make them as wide as possible.
* to summarize: take waterfalls as wide as possible, if you need more room, pad with zeroes. There are things in the "noise" that affect the measurement.

## oct 18
* no more regions file! for FRBs that are split for measurement write the region info into the right row in the csv

## oct 19
* 11G almost looks like a happy trombone. Hard to tell, the second burst is extremely low snr.. not even sure its there
	* needs noise removal
* when data at the edge the fit is misaligned (11O)

## oct 23
* FRB121102
	* 161 day periodicity (cruces 2020)
	* localized to milliarcsecond precision in 2018, in a dwarf galaxy at redshift 0.193 (see first para of cruces 2020)
	* z = 0.193

## oct 25
* dont put float logic in if statements use np.isclose

## oct 29
* redoing some aggarwal bursts to check effect of gui changes on measurement:
	* B016: has two bursts but both are hard to measure
	* B023: looks good but for the life of me I cant get a fit. Noise at the top of the burst
	* B037 shows time resolution affects the measurement. Presumably the higher the resolution the better. Will need to look at this effect in more detail.
	* slope measurements seem to agree if the resolutions are the same
	* B06 slope measurements agree well.
	* Duration measurements are different, for B006 its because of the theta fix. So I reran the model details for those results but it doesn't really matter since in the plot I do `computeModelDetails` every time, so as soon as I fixed it there it fixed it on the fly for the plot
* B025 at 561.5pc/cm3 Has a weird fit. Looks fine in GUI, looks clearly wrong in PDF, and doesn't fit the rest of the fits

## oct 31
* fix B025 fit, and gui was displaying the wrong fit because it was using the wrong dm list

## nov 1
* B028 is cutoff, maybe that's why its off trend?

## nov 2
* I did a test where I made a burst and then cutoff the top and measured the slope. The slope changed up to a quarter of its initial measurement (e.g. a slope measurement of -200 mhz/ms became -30 mhz/ms and a slope measurement of -8000 mhz/ms became -1000 mhz/ms). On that basis I think cutoff bursts should be excluded from the sample.

## nov 3
* B109 is a happy trombone with two bursts

## nov 4
* more negative slope is smaller in value? Noticed for B109, same deal in B006. Problem isn't in the Gajjar results

## nov 5
* angle seems to mean something different for the aggarwal dataset, but behaves normally for all the other datasets.

## nov 7
* forcing angle to be positive (the other datasets all have positive angles) doesn't fix the problem. Slope measurements are the same

## nov 8
* it might be more robust to fit a tan function to the slope measurements over the DM range.. but meh I dont know if that's worth the effort yet
* ~aggarwal set needs tan(-theta) not tan(theta)~
* investigate angle trend vs. aspect ratio. Based on the results I have i think its linear with equal aspect ratio but with a weird aspect ratio (like aggarwal's) its that weird kinked sigmoid thing
* add slope check to results pdf
* check the handedness of the angle between datasets (CCW or CW)

## nov 11
* handedness is CCW for both datasets. What is the aspect ratio between gajjar and aggarwal? gajjar is wider horizontally, aggarwal is wider vertically
* check slope with pdf across gajjar and aggarwal --> both look fine
* check full pdf

## nov 12
* different aspect ratios produce expected angle trend (ie some kind of tan(x))
* squaring the aspect ratio restores the trend
* in the aggarwal data the raw shape is 64x1220, but we clip to like 64x280
* i get good measurements if i decrease the shape from 64x1220 to 64x305. Why?
* angle offset controls slope trend but doesn't explain why when i went to 64x305 i get good measurements (ie. went very negative at vertical then very positive, through the tan asymptote)
* good measurements at 64x610, limited extent and full extent
* good measurements at 64x305, limited extent and full extent
* bad measurements at 64x1220, limited extent and full extent
* B011 has bad measurements and is downsampled by downf, downt = 2, 4, fchans, tchans = 32, 76, tsamp = 38
	* has sigmax < sigmay and sigmax > 0  and sigmay > 0 for all measurements
	* angle > 0
* B128 has good measurements has downf, downt = 2, 4, fchans, tchans = 32, 74, tsamp = 38
	* has sigmax > 0, sigma < 0 for all msrments and sigmax < abs(sigmay) for the first half, but then flips
	* angle < 0

## nov 15
* B006: slopes follow -tan when measured at 64x610 but not when at 64x1220
* B046: slopes follow -tan at 64x115
* 11D at any downt (full res to low res) produces a -tan trend in the slopes
* take new measurements of all bursts and reduce the time factor until you get the -tan trend. I suspect the very unequal aspect ratio (ie. way more time channels than frequency channels) is messing with the nonlinear lsq, but not sure.

## nov 16
* martin pointed out that the fits are always good, so there should be a corresponding correct slope calculation. Look at just one measurement and see what the correct slope is
* numerical precision limits the verticality of the slope, which is obviously apparent at big differences in aspect ratio: https://www.desmos.com/calculator/9nkuvrfzkb
	* the problem angles are at around pi/2, which if you look at the derivative of tan(x), blows up at pi/2
