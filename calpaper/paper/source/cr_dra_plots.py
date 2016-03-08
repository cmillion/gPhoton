"""
.. module:: make_crdra_plots.py

   :synopsis: Create plots of flare events on CR Draconis from gAperture output.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

import matplotlib.pyplot as pyp
from gPhoton.gphoton_utils import read_lc
from gPhoton.analysis.gtool_utils import calculate_caldat
import numpy

def make_plots():
    """
    Creates a plot of the flares in CR Draconis from gAperture.
    """

    ifile_nuv = "cr_dra_lc_nuv.csv"
    ifile_fuv = "cr_dra_lc_fuv.csv"
    # which time in the output file to use?
    time_to_use = "t_mean"

    # Read in data.  First is NUV, second is FUV.
    data_nuv = read_lc(ifile_nuv)
    data_fuv = read_lc(ifile_fuv)

    # Define the time ranges for the flares in GALEX time.
    flare_1 = [741111251., 741111981.]
    flare_2 = [799279801., 799280849.]
    flare_3 = [806670389., 806670917.]
    flare_4 = [837496202., 837497520.]
    flare_5 = [929935832., 929936740.]
    flare_6 = [956922970., 956924331.]
    flare_7 = [961164443., 961165519.]
    flare_8 = [991378742., 991379988.]

    all_flares = [flare_1, flare_2, flare_3, flare_4, flare_5, flare_6,
                  flare_7, flare_8]
    this_figure, these_subplots = pyp.subplots(nrows=8, ncols=1,
                                               figsize=(1024./96., 2048./96.),
                                               dpi=96.)

    # If set, will include FUV and Welsh et al. labels.
    show_blue = True

    # Name of output file.
    file_name = "flares.eps"

    # Define x-axis plot limits.
    xlims = [(0., 12.), (0., 17.), (0., 8.5), (0., 22.), (0., 15.), (0., 22.5),
             (0., 18.), (0., 21.)]

    # This was used when doing 2-column subplots, not really needed anymore,
    # but kept without having to re-write entire loop.
    row = 0

    for i in xrange(len(all_flares)):
        if i != 0:
            row = row + 1

        # Label the y-axis if this is zero'th subplot.
        if i == 0:
            this_figure.text(0.06, 0.5, 'Flux (ergs/cm$^{2}\!$/sec/$\AA$)',
                             ha='center', va='center', rotation='vertical',
                             size="xx-large")

        # Get the data for this flare based on time range.
        these_indexes = [ii for ii, x in enumerate(data_nuv[time_to_use]) if
                         x >= all_flares[i][0] and x <= all_flares[i][1]]

        # NUV plot times and fluxes, always present.
        these_times = data_nuv[time_to_use][these_indexes][1:-2]
        x_offset = these_times[these_times.keys()[0]]
        plot_times = [(x - x_offset)/60. for x in these_times]
        these_fluxes = data_nuv['flux_bgsub'][these_indexes][1:-2]

        # FUV plot times and fluxes, not always present.
        these_indexes_fuv = [ii for ii, x in enumerate(data_fuv[time_to_use]) if
                             x >= all_flares[i][0] and x <= all_flares[i][1]]
        if len(these_indexes_fuv) > 0:
            these_times_fuv = data_fuv[time_to_use][these_indexes_fuv][1:-2]
            plot_times_fuv = [(x - x_offset)/60. for x in these_times_fuv]
            these_fluxes_fuv = data_fuv['flux_bgsub'][these_indexes_fuv][1:-2]
            n_fluxes_fuv = len(these_indexes_fuv)
        else:
            these_fluxes_fuv = [numpy.nan] * len(these_fluxes)
            n_fluxes_fuv = 0

        # Plot NUV curve.
        these_subplots[row].plot(plot_times, these_fluxes, '-ko')

        # Plot FUV curve if there is one.
        if n_fluxes_fuv > 0:
            these_subplots[row].plot(plot_times_fuv, these_fluxes_fuv, '-bo')

        if i == 0 and show_blue:
            these_subplots[row].text(2., 2.15E-14, "FUV = blue", color='b')
        if i == 1 and show_blue:
            these_subplots[row].text(10., 2.0E-13, "Welsh et al. (2006) Flare",
                                     color='r')
        if i == 3 and show_blue:
            these_subplots[row].text(4., 1.4E-14, "FUV = blue", color='b')

        # Get the approximate calendar date of the event.
        mean_date = calculate_caldat(numpy.mean(these_times))
        these_subplots[row].set_title(mean_date, fontsize=12)

        # Label the x-axis, if this is the last subplot.
        if row == len(all_flares)-1:
            these_subplots[row].set_xlabel('Time Elapsed (Minutes)')

        # Make sure x-axis plot limits are enforced.
        these_subplots[row].set_xlim(xlims[i])

    pyp.subplots_adjust(hspace=0.30)
    pyp.savefig(file_name)

if __name__ == "__main__":
    make_plots()
