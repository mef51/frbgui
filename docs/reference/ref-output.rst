.. _outputcsv:

Output CSVs
===========

FRBGui outputs measurements as a comma-separated-value (CSV) spreadsheet along with a corresponding PDF that contains plots of each burst waterfall with its measurements overlaid.

The output spreadsheet is arranged as one measurement per row under the following columns:

.. note::
    * Plaintext columns like "slope (mhz/ms)" indicate a measurement.
    * *Italicized* columns indicate information about how the measurement was taken, such as "time_res (s)".
    * :strike:`Stricken` column names indicate a deprecated column and normally should not be used.

.. list-table::
    :header-rows: 1

    *   - Column Name
        - Description
    *   - name
        - The name of the burst file measured, potentially suffixed with a letter (e.g., “a”, “b”, “c”) denoting the sub-pulse
    *   - DM
        - The particular trial DM used to perform measurements
    *   - center_f
        - The center frequency of the burst as measured by frbgui
    *   - center_f_err
        - The propagated uncertainty on center_f. For sub-pulses, this value will be the same across DMs, as the waterfall noise needed to compute this value will be sampled once from the full waterfall.
    *   - slope (mhz/ms)
        - The sub-burst slope, obtained from the cotangent of the fit angle
    *   - slope error (mhz/ms)
        - The propagated uncertainty on the slope
    *   - theta
        - The angle of the semi major axis of the fit ellipse measured counter-clockwise from the positive y-axis
    *   - red_chisq
        - Reduced chi-squared indicating the goodness-of-fit of the 2D Gaussian to the 2D ACF
    *   - amplitude
        - Amplitude of the 2D Gaussian fit. Not normalized.
    *   - xo
        - Central x-position (time) of the 2D gaussian fit to the 2D-ACF
    *   - yo
        - Central y-position (frequency) of the 2D gaussian fit to the 2D-ACF
    *   - sigmax
        - Standard deviation in x (time if vertical) of the 2D gaussian fit to the 2D-ACF
    *   - sigmay
        - Standard deviation in y (frequency if vertical) of the 2D gaussian fit to the 2D-ACF
    *   - angle
        - The fit angle. Equivalent to theta up to a difference of pi/2 due to the ambiguity between axes when fitting
    *   - amp_error
        - The fit uncertainty on the amplitude
    *   - xo_error
        - The fit uncertainty on xo
    *   - yo_error
        - The fit uncertainty on yo
    *   - sigmax_error
        - The fit uncertainty on sigmax
    *   - sigmay_error
        - The fit uncertainty on sigmay
    *   - angle_error
        - The fit uncertainty on angle
    *   - slope_abs
        - The absolute value of slope (mhz/ms). If negative, this indicates a potential measurement issue
    *   - slope_over_nuobs
        - 'slope_abs’ / ‘center_f’. In the TRDM, this is proportional to 1/’tau_w_ms’
    *   - slope_over_nuobs_err
        - The propagated 2D gaussian fit uncertainty for ‘slope_over_nuobs’.
    *   - recip_slope_over_nuobs
        - The reciprocal of slope_over_nuobs
    *   - slope_abs_nuobssq
        - 'slope_abs’ / ‘center_f**2’. In the TRDM, this is proportional to a constant.
    *   - min_sigma
        - Same as either sigmax or sigmay, whichever is smaller
    *   - max_sigma
        - Same as either sigmax or sigmay, whichever is larger
    *   - min_sigma_error
        - Same as either sigmax_error or sigmay_error, associated with whichever sigma is smaller
    *   - max_sigma_error
        - Same as either sigmax_error or sigmay_error, associated with whichever sigma is larger
    *   - tau_w
        - The sub-burst duration in milliseconds as defined in Chamma et al. (2022).
    *   - tau_w_error
        - The propagated uncertainty on tau_w (ms)
    *   - tau_w_ms
        - The sub-burst duration in milliseconds as defined in Chamma et al. (2022). An alias of tau_w. Units are implied by the choice of coordinates when fitting
    *   - bandwidth (mhz)
        - Best-fit bandwidth for the sub-pulse in MHz
    *   - bandwidth error (mhz)
        - Propagated uncertainty for the bandwidth in MHz
    *   - *f_res (mhz)*
        - The frequency resolution (MHz) of the final waterfall used for the FRBGUI measurements
    *   - *time_res (s)*
        - The time resolution (seconds) of the final data array used for the FRBGUI measurements
    *   - *downf*
        - The factor by which the file was downsampled in frequency from FRBGUI from the original waterfall
    *   - *downt*
        - The factor by which the file was downsampled in time from FRBGUI from the original waterfall
    *   - *fchans*
        - The number of frequency channels remaining after the downsample in FRBGUI
    *   - *tchans*
        - The number of time channels remaining after the downsample in FRBGUI
    *   - *tsamp_width*
        - The number of time channels set by the user in FRBGUI used for measurement (i.e., the width of the waterfall)
    *   - *subbg_start (ms)*
        - The time, in ms, from the start of the file that the user-defined background sample begins
    *   - *subbg_end (ms)*
        - The time, in ms, from the start of the file that the user-defined background sample ends
    *   - *sksigma*
        - The sigma chosen for the SK-SG filter, which performs RFI removal
    *   - *skwindow*
        - The window chosen for the SK-SG filter, which performs RFI removal
    *   - *regstart_a*
        - When using regions, the start of the “a” subpulse region in ms
    *   - *regend_a*
        - When using regions, the end of the “a” subpulse region in ms
    *   - *regstart_b*
        - When using regions, the start of the “b” subpulse region in ms
    *   - *regend_b*
        - When using regions, the end of the “b” subpulse region in ms
    *   - *background*
        - The size of the background subpulse region in ms. Used to pad sub-pulses with zeroes, which improves measurement stability. Background region is assumed to start from 0 ms
    *   - :strike:`sigma_t`
        - The product of min_sigma and time_res (s). This is a (poor) measure of burst duration. Use tau_w or tau_w_ms for duration instead
    *   - :strike:`tau_w_old`
        - The sub-burst duration as defined and used in Chamma et al. (2021). Due to the choice of coordinates when performing fits, this form can be subject to larger uncertainties when the burst is near vertical
    *   - :strike:`t_from_sigma`
        - The product of min_sigma and sin(theta). This is a (poor) measure of burst duration when finding fits with physical coordinates
    *   - :strike:`sigma_t_ms`
        - 'sigma_t’ in ms. See sigma_t
    *   - :strike:`tau_w_ms_old`
        - tau_w_old in milliseconds

Special thanks to Dr. Sofia Sheikh for contributions to this table.
