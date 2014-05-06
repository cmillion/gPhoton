#### Computing deadtime:
Two main ways to estimate the deadtime, each of which is implemented in the compute_deadtime() funciton in CalibrationTools.py. The first "correct" way compares the measured stim count rate against the nominal commanded rate. This doesn't work very well for short time scales, though, because it's subject to Poisson statistics. The second "less correct" way uses a conversion fact between global detector count rate and deadtime. The second method is far more consistent across exposure times (because the total number of detector counts, even per second, is quite high.
The implementation in compute_deadtime() works straight from the photon CSV files so we will have to adapt them to work with the database.
It would also be useful to recompute the scale fact for Method 2. Patrick Morissey, who did the original work with Tim Conrow, said that the scale factor is slightly different in each band, but they forced them to be the same for consistency. They also only used a subset of the total mission data, and he doesn't have any of their original code left.
We should test the possibility that the scale factor for the global dead time may have changed over the mission.

*Note: The CSP occurred on May 4, 2010. The post-CSP TAC switch occurred on eclipse 38150 (June 23, 2010), which corresponds to a time of t0=881881215.995.*

---
`from exptime import netdead
from dbasetools import fGetTimeRanges

band = 'NUV'

skypos = [294.86493, 47.30736]
tranges = fGetTimeRanges(band,skypos)
deadtime = [netdead(band,trange[0],trange[1]) for trange in tranges]

skypos = [323.06766667,0.25400000]
tranges = fGetTimeRanges(band,skypos)
deadtime0 = [netdead(band,trange[0],trange[1]) for trange in tranges[np.random.random_integers(0,tranges.shape[0]-1,100)]]

%pylab
# Two methods vs. global count rate
plt.plot(np.array(deadtime)[:,-1]/(np.array(deadtime)[:,1]-np.array(deadtime)[:,0]),np.array(deadtime)[:,2],'.')
plt.plot(np.array(deadtime)[:,-1]/(np.array(deadtime)[:,1]-np.array(deadtime)[:,0]),np.array(deadtime)[:,3],'.')

# Two methods vs. each other
plt.plot(np.array(deadtime)[:,2],np.array(deadtime)[:,3],'.')`

