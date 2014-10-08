from utils import *

skypos = {'LDS749B':[323.06766667,0.25400000],
          'PS_CDFS_MOS00':[53.1032558472, -27.7963826072],
          'CDFS_00':[53.1273118244, -27.8744513656]}
bands = ['NUV','FUV']

ra0,dec0=skypos['LDS749B']
band = bands[1]
outfile = 'test.csv'
datamaker(bands[1],skypos['LDS749B'],outfile)


#######################################


# Set up the plot
plt.figure()
plt.gca().invert_yaxis()
plt.title(str(src)+' - '+str(band))
plt.xlabel('Exposure time (s)')
plt.ylabel('AB Magnitude')
# Make the y dimensions a little bigger than the data
plt.axis([0,1600.,data['mag'].min()-0.01,data['mag'].max()+0.01])
# Plot the reference magnitude along with upper and lower bounds vs. expt
top,bot=gu.model_errors(refmag['FUV'],band,sigma=5)
plt.plot(top)
plt.plot(bot)
plt.plot(np.ones(1600)*refmag['FUV'])
# Plot the MCAT data
ix = np.where(mcat['fmag']>0)


#####################################
# Set up the plot
plt.figure()
plt.gca().invert_yaxis()
plt.title(str(src)+' - '+str(band))
plt.xlabel('Exposure time (s)')
plt.ylabel('AB Magnitude')
# Make the y dimensions a little bigger than the data
plt.axis([0,1600.,data['mag'].min()-0.01,data['mag'].max()+0.01])
