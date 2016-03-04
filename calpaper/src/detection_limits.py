import gPhoton.galextools as gt

magrange = np.arange(14,24,1)

expt_ratio = 0.8 # Estimate of ration of effective to raw exposure time
t_raw = np.arange(1,1600,1)
t_eff = expt_ratio*t_raw

for sigma in [1,3]:
    fig = plt.figure(figsize=(8*2,6*2))
    for i,band in enumerate(['FUV','NUV']):
        plt.subplot(2,1,i)
        plt.grid(b=True)
        plt.semilogx()
        plt.xlim(1,1600)
        plt.title('{b} Detection Threshholds'.format(b=band))
        plt.xlabel('Exposure Bin Depth (s)')
        plt.ylabel('{n} Sigma Error (AB Mag)'.format(n=sigma))
        for mag in magrange:
            cps = gt.mag2counts(mag,band)
            cps_err = sigma*np.sqrt(cps*t_eff)/t_eff
            mag_err = mag-gt.counts2mag(cps+cps_err,band)
            plt.plot(t_raw,mag_err,label=mag)
        plt.legend()
    plt.savefig('GALEX {n} Sigma Detection Limits.png'.format(n=sigma))
