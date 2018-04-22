# moment-methods-for-foregrounds

We use the moment method to describe the SED 

# Data compression

Sum of SED of modified blackbodies derived from sum of modified black body SED's with three different distribution of parameters T and the slope. For the slope "alpha" we draw from a Gaussian with a mean at 1. For temperature we draw from a two Gaussians centered at 9.75 K and 15.7 K. For "NARROW" and "WIDE" SED models we basically draw from distributions with narrow and larger widths respectively. The "DELTA" SED model corresponds to a modified black body with slope=1 and T=15.7 K. The specific parameter distributions and the resultant SED's are depicted in the figures below:

<img src="https://github.com/adityarotti/moment-methods-for-foregrounds/blob/master/figures/data_compression/model_sed.jpeg" alt="Alpha distribution" width="33%" ></img>
<img src="https://github.com/adityarotti/moment-methods-for-foregrounds/blob/master/figures/data_compression/T_distribution.jpeg" alt="T distribution" width="33%"></img>
<img src="https://github.com/adityarotti/moment-methods-for-foregrounds/blob/master/figures/data_compression/alpha_distribution.jpeg" alt="Alpha distribution" width="33%"></img>


The figures below indicate how well differnet order moment fits perform on the different SED models:

<img src="https://github.com/adityarotti/moment-methods-for-foregrounds/blob/master/figures/data_compression/improvement_in_fits_delta_sed.jpeg" alt="Alpha distribution" width="33%"></img>
<img src="https://github.com/adityarotti/moment-methods-for-foregrounds/blob/master/figures/data_compression/improvement_in_fits_narrow_sed.jpeg" alt="Alpha distribution" width="33%"></img>
<img src="https://github.com/adityarotti/moment-methods-for-foregrounds/blob/master/figures/data_compression/improvement_in_fits_wide_sed.jpeg" alt="Alpha distribution" width="33%"></img>

The figures below show the relative errors

<img src="https://github.com/adityarotti/moment-methods-for-foregrounds/blob/master/figures/data_compression/relative_error_with_taylor_order_delta_sed.jpeg" alt="Alpha distribution" width="33%"></img>
<img src="https://github.com/adityarotti/moment-methods-for-foregrounds/blob/master/figures/data_compression/relative_error_with_taylor_order_narrow_sed.jpeg" alt="Alpha distribution" width="33%"></img>
<img src="https://github.com/adityarotti/moment-methods-for-foregrounds/blob/master/figures/data_compression/relative_error_with_taylor_order_wide_sed.jpeg" alt="Alpha distribution" width="33%"></img>


# Moment analyis on SED's provided by Colin Hill and Marcelo Alvarez

<img src="https://github.com/adityarotti/moment-methods-for-foregrounds/blob/master/figures/colin_marcelo_sed/intensity_spectrum_data.jpeg" alt="Alpha distribution" width="33%"></img>

Results from analysis on Colin's SED performed by interpolating on the SED:

<img src="https://github.com/adityarotti/moment-methods-for-foregrounds/blob/master/figures/colin_marcelo_sed/improvement_in_fits_colin_sed.jpeg" alt="Alpha distribution" width="33%"></img>
<img src="https://github.com/adityarotti/moment-methods-for-foregrounds/blob/master/figures/colin_marcelo_sed/relative_error_with_taylor_order_colin_sed.jpeg" alt="Alpha distribution" width="33%"></img>
<img src="https://github.com/adityarotti/moment-methods-for-foregrounds/blob/master/figures/colin_marcelo_sed/inu_fit_with_taylor_order_colin_sed.jpeg" alt="Alpha distribution" width="33%"></img>

Results from analysis on Marcelo's SED performed by interpolating on the SED:

<img src="https://github.com/adityarotti/moment-methods-for-foregrounds/blob/master/figures/colin_marcelo_sed/improvement_in_fits_marcelo_sed.jpeg" alt="Alpha distribution" width="33%"></img>
<img src="https://github.com/adityarotti/moment-methods-for-foregrounds/blob/master/figures/colin_marcelo_sed/relative_error_with_taylor_order_marcelo_sed.jpeg" alt="Alpha distribution" width="33%"></img>
<img src="https://github.com/adityarotti/moment-methods-for-foregrounds/blob/master/figures/colin_marcelo_sed/inu_fit_with_taylor_order_marcelo_sed.jpeg" alt="Alpha distribution" width="33%"></img>
