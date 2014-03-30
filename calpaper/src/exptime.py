# Computing deadtime:
#  Two main ways to estimate the deadtime, each of which is implemented in the
#  compute_deadtime() funciton in CalibrationTools.py. The first "correct" way
#  compares the measured stim count rate against the nominal commanded rate.
#  This doesn't work very well for short time scales, though, because it's
#  subject to Poisson statistics. The second "less correct" way uses a conversion
#  fact between global detector count rate and deadtime. The second method is
#  far more consistent across exposure times (because the total number of
#  detector counts, even per second, is quite high.
#
#  The implementation in compute_deadtime() works straight from the photon CSV
#  files so we will have to adapt them to work with the database.
#
#  It would also be useful to recompute the scale fact for Method 2. Patrick
#  Morissey, who did the original work with Tim Conrow, said that the scale
#  factor is slightly different in each band, but they forced them to be
#  the same for consistency. They also only used a subset of the total mission
#  data, and he doesn't have any of their original code left.
#
#  We should test the possibility that the scale factor for the global dead time
#  may have changed over the mission.

refrate     = 79. # counts per second, nominal stim rate
feeclkratio = 0.966 # not entirely sure what detector property this adjusts for
tec2fdead   = 5.52e-6 # Conversion from TEC to deadtime correction (Method 2)

# We need an SQL query that will pull the stim count out of the database
stim_events = 0.

# Method 0
dead0 = tec2fdead*(len(t)/exptime)/feeclkratio

# Method 1 - Use a histogram to ignore crazy values.


# Method 2
# gQuery.deadtime1() + gQuery.deadtime2() = total global counts
# gQuery.deadtime() = deadtime correction using Method 2

