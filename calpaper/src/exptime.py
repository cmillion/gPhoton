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
