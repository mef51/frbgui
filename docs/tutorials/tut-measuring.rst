.. _tutmeasure:

Measuring FRBs
==============

This tutorial will describe how to perform measurements of FRB waterfalls using two examples. At the end of the tutorial we will have produced a CSV spreadsheet with all of our measurements that can be used for further analysis and a PDF that displays plots of each burst with each of its measurements overlaid on top for review. This tutorial is intended for those with some knowledge on the observational characteristics of FRBs and with basic experience programming and using a terminal. For an introduction to FRBs, see (e.g.) `Petroff et al. (2022) <https://link.springer.com/article/10.1007/s00159-022-00139-w>`_

The first FRB is a simple bright pulse with some noise in its waterfall, and the second FRB consists of multiple sub-bursts at a high data resolution that we would like to separate and measure independently.

These examples will broadly demonstrate the capabilities of FRBGui and what a typical workflow looks like.

To follow along, you may download the two burst files which have been prepared in ``.npz`` format for FRBGui here (see :ref:`tutpreparation`):

	* `aggarwal2021_B006.npz <../aggarwal2021_B006.npz>`_
	* `gajjar2018_11A.npz <../gajjar2018_11A.npz>`_

The first burst is burst B006 published in Aggarwal et al. (2021) and discovered in data taken by the Arecibo telescope while observing the repeating source FRB 20121102A. With FRBGui installed (see :ref:`installation`), navigate to the folder where your bursts have been downloaded and saved and run the following command in the terminal::

	frbgui

This will start the graphical user interface and you should see your first burst already loaded and displayed:

.. setting width will tell sphinx doc to generate a link for the image. So I pick a large number
.. figure:: /imgs/tut-meas-1.png
	:width: 1000

	The main FRBGui window. Click to enlarge.

The main interface consists of the "FRB Analysis" window which is where we will manipulate the waterfall and the "FRB Plots" window, which displays the FRB waterfall (left), the two-dimensional autocorrelation function (ACF) on the right, and the frequency integrated time-series of the waterfall at the bottom.

At this stage it is important to verify that the frequency and time axes are correctly displayed, usually by cross-referencing with the corresponding figure in the publication. If the axes of the waterfall are incorrect, this can indicate an issue with the way the ``.npz`` file was prepared, either with the frequencies, bandwidth, or duration supplied and should be double-checked against FRBGui's :ref:`burstformat`.

If the axes appear correct we can proceed with measurement.

For this burst we want to measure its sub-burst slope, or the rate of change of frequency with time, a quantity related to the drift rate and a characteristic feature of FRBs. Because the choice of Dispersion Measure (DM) affects the value of the sub-burst slope by "rotating" the burst as it appears in its waterfall we will repeat our measurements over a range of DMs. These measurements can then be used at a later stage of analysis to estimate and characterize the uncertainty on our measurements.

Measurement Properties
----------------------

Under the fold-out section "1. Data" in the FRB Analysis window you will find the Burst Metadata subsection. In here we see that the DM of the burst is 562.056 pc/cm :math:`^3`. If needed we could change this DM and save it back into the burst's ``.npz`` file using the +/- keys next to the DM and the grayed out "Save" button. We will leave this value as is for this burst.

.. figure:: /imgs/tut-meas-2.png
	:width: 500

Turning our attention to the "Display Width" field, we see the value 150, which is the number of time channels that are currently being displayed in the waterfall plot. This default value is not the full size of the waterfall saved in the ``.npz``, and we can increase this value to display more channels.

Since this burst will rotate as it is measured at different DMs, we want to ensure that there is enough room for it to do so:

**>>** Increase the display width for this burst to 200 channels by typing in the value manually or by clicking the +/- keys, which increments the value by 10 channels each time.

Below the Display Width we can set the Dedispersion Range and the # of Trial DMs which will determine the range of DMs that measurements are repeated over and the number of measurements obtained. The pair of values are displaying a default range of 556.435 to 567.677 pc/cm :math:`^3`. In the interface, the value on the left is the start of the DM range and the value on the right is the end.

**>>** Change the DM range to 555 to 575 pc/cm by **double-clicking** the values in the range input and set the DM step to 0.5 pc/cm :math:`^3`. This will tell FRBGui to perform measurements over the DMs 555.0, 555.5, 556.0, ... 574.5, 575.0 pc/cm :math:`^3`, as well as the burst DM from the loaded file, for a total of 41+1 measurements.

.. figure:: /imgs/tut-meas-3.png
	:width: 500

.. note::
	FRBGui will not allow you to specify the start of the range to be greater than the end of the range. It will overwrite the value for the start of the range if it detects this error. Therefore it is best to input the end of the DM range before inputting the beginning.

	In addition, if the DM range inputted does not include the burst's DM as read from the loaded file, a warning will be displayed. Ensure your chosen DM range includes the DM of the burst.

Cleanup (Background Subtraction)
--------------------------------

At this stage we will now look at the waterfall and its 2D autocorrelation function (ACF) in the plot windows in order to assess what kind of cleanup will be necessary before measurement. The burst now appears as in the following:

.. figure:: /imgs/tut-meas-4.png
	:width: 650

Notice in the waterfall plot, some channels have already been masked, and there is little radio frequency interference (RFI) remaining. RFI typically appear as bright horizontal bars across the waterfall. The data resolution is good and the SNR of the burst is high, as can be seen by the contrast between the colors of the bright burst signal and the background noise.

While the burst appears well defined, its ACF however appears blocky and filled with artifacts. Our goal is for the ACF to resemble an ellipse as much as possible. The artifacts in the ACF will be caused by irregularities in the burst waterfall. In this case, the masked out noise channels may not be properly zeroed out, resulting in the blocky structures we see.

To fix this we will use the "Subtract background" feature under the second fold-out section "2. Waterfall Cleanup". Subtract background will take a time averaged sample of the noise as a function of frequency and subtract this sample from each column of data, and is a good and simple way of removing frequency dependent noise that does not change in time from a waterfall. FRBGui allows you to pick the start time (t_start (ms)) and end time (t_end (ms)) of this background sample. Typically you will choose the first few milliseconds of a burst waterall before the burst starts as the background.

**>>** Check the "Subtract background sample" checkbox and set the second value of the number pair (which corresponds to t_end (ms)) to 5.0 ms. Feel free to scrub this value with your mouse to dynamically see the effect this has on the waterfall.

.. figure:: /imgs/tut-meas-5.png
	:width: 650

Immediately we see the colors of the burst waterfall update and the artifacts disappear from the ACF leaving behind a well-defined ellipse with a high contrast with its background, indicating a high signal-to-noise. Since there are no components to seperate, we can skip the third fold-out section "3. Burst Splitting". This burst is now ready for measurements.

.. figure:: /imgs/tut-meas-6.png
	:width: 650

Performing Measurements
-----------------------

**>>** In the last fold-out section "4. Spectro-Temporal Measurements", click "Measure over DM Range".

.. figure:: /imgs/tut-meas-7.png
	:width: 650

FRBGui will now begin the process of incoherently dedispersing the waterfall to each DM we specified in the DM range and fitting a 2D gaussian the ACF of the burst at that DM. From the parameters of the 2D gaussian, spectro-temporal properties such as the sub-burst slope, duration, bandwidth, and center frequency are computed (See :ref:`acfmethod`). Leaving the "Repeat Measurements" box checked (recommended) will tell FRBGui to perform two passes: the first finds a rough fit and then the rough fit is used to refine the fit further in the second pass to improve accuracy.

As FRBGui performs measurements the Status bar next to the "Measure..." button will fill up until measurements are complete. When done, a table with the results (the "results table") will be displayed in the fold-out section along with information on

* the total number of measurements
* the burst displayed and the DM it is displayed at

.. figure:: /imgs/tut-meas-8.png
	:width: 650

.. note::

	If FRBGui crashes for whatever reason, any measurements you have made will have automatically been saved in the ``backups`` folder FRBGui creates in the same directory you started FRBGui in. You can load the spreadsheet back into FRBGui and continue where you left off from using the "Load Results" button in the "4. Spectro-Temporal Measurements" section or continue afresh knowing that your measurements have already been saved.

Reviewing Measurements
----------------------

With our measurements complete, we will now review our measurements to ensure that the fits found by FRBGui are accurate to the ACF of the burst. With the measurements complete, the plot figures now show a blue elliptical contours overlaid onto the ACF and a dashed line that corresponds to the sub-burst slope found, as shown in the image below. The blue contours correspond to 90% and 25% of the peak of the 2D Gaussian model fit to the ACF.

.. figure:: /imgs/tut-meas-9.png
	:width: 1000

By clicking on rows of the results table you can choose which measurement is displayed on the plots. In addition, just above the results table are the left ◀ and right ▶ triangle scroll buttons, which can be clicked to go through each measurement one by one.

**>>** Click on the first row of the results table to display the first measurement. Check that the blue contours and dashed line line up well with the plot of the burst ACF as shown in the image above.

**>>** Click the ▶ button to load the each measurement one-by-one and note whether each contour aligns with its autocorrelation. The currently displayed DM will be listed in the "DM Displayed:" status text just above the results table.

Once you have reached the measurements for the DM of 575 pc/cm :math:`^3` you may notice the burst tilted so much that it has begun to roll around to the other side of its plot and blobby artifacts appearing on its ACF:

.. figure:: /imgs/tut-meas-10.png
	:width: 1000

In typical analyses this burst measurement would be excluded on the basis of it being greatly over dedispersed (as evidenced by the characteristic :math:`1/\nu^2` curve) and even excluded on the basis of its positively sloping sub-burst slope, which is usually assumed to be unphysical.

However, we may still want to ensure our set of measurements is accurate so we will adjust our measurement properties and re-measure to correct this issue.

**>>** Scrolling back up to the first fold-out section, increase the "Display Width" further under "Burst Metadata" to 360 channels. Notice the burst appears narrower as FRBGui loads more data into the displayed plot

**>>** Return to "4. Spectro-Temporal Measurements" and click "Measure over DM range" to repeat the measurements with the new measurement properties.

When reviewing the measurements now at a DM of 575 pc/cm :math:`^3` we should see a clean ACF without artefacts and a more accurate measurement:

.. figure:: /imgs/tut-meas-11.png
	:width: 1000

Manual Initial Fit Guess
------------------------

In this case the fitting algorithm found accurate measurements automatically without the need for us to input an initial approximate guess manually.

However, if needed, we could input an initial guess by expanding the "Fit initial Guess" section just above the results table, and checking "Use Initial Guess". This adds a green contour set based on the best fit so far to the ACF plot that represents our initial guess. We can modify its parameters using the Amplitude, Angle, x0, y0, sigmax, and sigmay inputs in the interface.

.. figure:: /imgs/tut-meas-12.png
	:width: 1000

To perform a measurement with our initial guess, we click "Redo this DM", which will redo the fit for the currently displayed measurement. Checking "Use for all DMs" will automatically repeat the measurement using the initial guess one-by-one and display the result as it works. This allows you to monitor the measurements as they are found.

For this burst, the fits found automatically are accurate without the need for a manual initial guess. The measurement of this burst is now complete and we can proceed to our second burst.

Downsampling an FRB
-------------------

The next FRB we will measure is burst 11A published in Gajjar et al. (2018) which was observed by the Green Bank Telescope from the same source as before, FRB 20121102A, using a 4--8 GHz receiver. This burst is at a much higher data resolution (more frequency channels and more time channels) and is made up of four sub-bursts. That is, four distinct pulses are seen spanning a duration of just 3 milliseconds.

**>>** Scroll to the top of the FRB analysis window and select the burst "gajjar2018_11A.npz" to load and display the burst.

.. figure:: /imgs/tut-meas-13.png
	:width: 1000

	FRBGui with burst 11A loaded

The burst is loaded and we see that the contrast in the figures is a little lower. This may indicate a low signal-to-noise ratio, but we see the burst pulses prominently in the time series plot below the waterfall so this is more likely due to the choice of color scale. Note that the color scale of the waterfall and ACF plots can be adjusted with the sliders just above each respective plot.

Nonetheless we can increase signal to noise and decrease computation time by downsampling the waterfall, i.e. by decreasing the number of frequency and time channels by averaging them together. While it is best to perform your measurements at full resolution when accuracy is a concern, many bursts may be adequately measured at a lower resolution and much more quickly by downsampling. For bursts with very low signal to noise, fits may be difficult to find, and downsampling is an essential strategy for boosting the signal to noise and securing a good fit.

The DM range settings persist from the first burst so we will move immediately to the waterfall cleanup stage, and downsample this waterfall. This will increase the signal-to-noise but in this case also make it a little easier to manipulate settings in the interface, as this large burst can slow down the response time of the interface.

**>>** Under "2. Waterfall Cleanup > Downsampling", set ``numt`` to 1024 and ``numf`` to 845. These correspond to the number of time and frequency channels in the downsampled waterfall.

"Original Size" shows a tuple of the waterfall's shape (frequency channels, time channels). The number you set must evenly divide the original shape. For example, since the original waterfall has 2048 channels, only 1024, 512, 256, 128, 64, etc. are valid inputs. Likewise with 1690, valid inputs are 845 (1690/2), or 338 (1690/5) (for example).

When downsampling the width of the waterfall may change as FRBGui replots the waterfall, so typically there is no need to change the "Display Width" until you have downsampled. In this case, the last pulse is cutoff, so we need to adjust the width.

**>>** Adjust the "Display Width" to show 150 channels after downsampling to 845 by 1024 channels.

.. figure:: /imgs/tut-meas-14.png
	:width: 600

Your waterfall now appears as in the following and you can notice a higher contrast in the colors of the waterfall and ACF plot after downsampling.

.. figure:: /imgs/tut-meas-15.png
	:width: 1000

	Click to enlarge

Mask Ranges
-----------

It will be difficult to see, but in the enlarged view of the above image you may notice a faint but distinct thin line streaking the waterfall at around 6600 MHz. This noise feature is minor and will likely not affect measurement, however we can zero it out using the "Mask Range" feature under "Masking" in the "2. Waterfall Cleanup" fold-out.

**>>** Under "Masking", click "Add Mask Range" to add a pair of number fields for specifying a mask range. The left and right values are the starting and ending frequency channels of the mask range. Drag the fields to change the values and notice the blank horizontal gap that appears in the waterfall. Drag the values until the mask range overlaps with the thin narrow feature, or input 880 and 905 into the fields.

.. figure:: /imgs/tut-meas-16.png
	:width: 600

After masking your waterfall will appear with a zeroed out band.

.. figure:: /imgs/tut-meas-17.png
	:width: 600

The mask range can be removed by selecting the "X" button next to it.

Burst Splitting
---------------

We may now choose to proceed with measurement, in which case the Gaussian fit will fit the ACF that is formed by the entire FRB, that is, the 4 sub-bursts that comprise it. This measurement is important and will include the so-called drift rate, or the change in frequency between consecutive resolved bursts. FRB drift rates obey characteristic relationships that are potentially a clue to the underlying emission mechanism that remains a mystery. However, we may also want to measure the change in frequency with time of the individual sub-bursts to obtain a sub-burst slope (also called intra-burst drift rate) as well as the individuated spectro-temporal properties (such as the durations of each pulse).

We can use FRBGui's burst splitting feature to encode the start and end times of each burst and FRBGui will automatically measure each individual sub-burst we specify after measuring the entire FRB. This is done by sectioning the waterfall into regions.

**>>** Under "3. Burst Splitting" check the "Are there multiple bursts in this waterfall" box. The grayed out fields will be enabled.

These fields consist of the beginning and end start times (the left and right values, respectively) as well the type of region you are specifying, either "Background" or "Burst". The first region will default to the "Background" type. Since the spaces between sub-bursts can be quite short, when FRBGui measures sub-bursts it will add a zero-padded region to the left and right of the sub-burst region before measuring. This "Background" region is the amount of zero-padding to add. Typically you can simply specify the region before the burst starts.

.. figure:: /imgs/tut-meas-18.png
	:width: 600

**>>** Slide the end of the "Background" region until the line that has appeared on the waterfall representing it is near the beginning of the burst, around 2.6 ms.

**>>** Click "Add Region" to add a new region. This second region will default to "Burst". Slide the end of this region to the very end of the first pulse. Use the time series plot beneath the waterfall to line up the region line that appears to the valley between the two pulses.

.. figure:: /imgs/tut-meas-19.png
	:width: 550

**>>** Continue adding regions until you have specified one background region and the four sub-bursts. Notice that as you add "Burst" regions, the end of the previous region is automatically set as the start of the new region.

.. note::

	When splitting sub-bursts it is important to specify the region around the sub-burst to be as large as possible while minimizing overlap. This improves measurement stability and reproducibility.

When finished your regions and waterfall may appear as in the following image:

.. figure:: /imgs/tut-meas-20.png
	:width: 600

.. figure:: /imgs/tut-meas-21.png
	:width: 550

Saving and Exporting Results
----------------------------

**>>** With the regions specified you click "Measure over DM Range" in the "4. Spectro-Temporal Measurements" section and sit back while FRBGui first obtains measurements for the entire waterfall over the DM range and then each of the four individual sub-burst you have specified. The terminal you started may offer additional feedback on the progress as this process can take a few minutes.

When complete all the measurements will again be displayed in a results table. This time, the table includes the sub-bursts suffixed by "_a", "_b", "_c", "_d" corresponding to each of the sub-burst regions specified. Note the number of measurements is 210 which is equal to the 1+4 regions specified times 42 measurements per region.

.. figure:: /imgs/tut-meas-22.png
	:width: 600

When reviewing measurements, notice that clicking on the row of a sub-burst measurement, such as "gajjar2018_11A_a", will display the extracted waterfall of the first isolated sub-burst region along with its own ACF. In this way you may verify the fit accuracy of the sub-bursts in addition to the fits of the entire waterfall.

**>>** Adjust the column sizes of the tables until you can clearly see the burst name and DM. Scroll down and select one of the measurements for "gajjar2018_11A_a". Notice the zero padded vertical regions in the waterfall when viewing sub-bursts. Click rows or use the ◀ and  ▶ buttons to review all the measurements in the table.


.. figure:: /imgs/tut-meas-23.png
	:width: 1000

Since we have finished the list of bursts we wanted to measure, we may now export the measurements to a CSV spreadsheet as well as a PDF of the measurements.

**>>** At the very bottom of the "FRB Analysis" window below the results table, enter a filename prefix (such as "B006_11A_results") and click "Export Results PDF". This automatically exports a CSV and then plots each measurement and saves it in a corresponding PDF. Note that saving a PDF with many measurements can take a lot of time. If you don't want the PDF you can just click "Export Results CSV".

.. figure:: /imgs/tut-meas-24.png
	:width: 1000

The results files will be saved in the same directory you started FRBGui in. The CSV file also doubles as a save file of your FRBGui session since all the information needed to reproduce your measurements will be saved in the CSV. You can reload your FRBGui session by starting FRBGui in the same directory, clicking "Load Results", and selecting the CSV you previously exported.

Example of the CSV Output. See :ref:`outputcsv` for the meanings of the columns in detail.

.. figure:: /imgs/tut-meas-25.png
	:width: 600

Example of the PDF Output. Note that the waterfall plots here show the sub-burst slope or drift rate (dashed line), the duration (the correlation length) with the horizontal bar, and the frequency bandwidth with the vertical bar.

.. figure:: /imgs/tut-meas-26.png
	:width: 1000

In this tutorial we have taken two FRBs and used FRBGui to prepare the bursts for measurement and measure their spectro-temporal properties, producing CSV and PDF outputs that are useful for later analysis outside of FRBGui and for reviewing measurements visually.

Thank you for following along this tutorial and best of luck on your FRB measuring days.

