This is mostly a list of known GALEX flares from the literature that we're able to reproduce with gAperture, along with the commands to do so.

One of the things that seems to happen in the literature is that the source positions get rounded off enough that they aren't really useful for positioning the aperture anymore, and we needed to go back to the MCAT to find the most likely nearby source. RA and Dec should really be reported with 5 or 6 significant figures (in decimal degrees).

####AF Psc
#Reported in [1].

python ./gAperture.py -b 'FUV' --skypos '[352.937473,-2.744136]' -s 10 -f 'AFPsc_FUV.csv' -a 0.002 -i 0.004 -o 0.006 -v 0 --maxgap 1000 --overwrite

python ./gAperture.py -b 'NUV' --skypos '[352.937473,-2.744136]' -s 10 -f 'AFPsc_NUV.csv' -a 0.002 -i 0.004 -o 0.006 -v 0 --maxgap 1000 --overwrite

####CR Dra
#Reported in [1]. Note: Only one flare was reported in the literature. The data appear to contain at least three and maybe four distinct flares.

python ./gAperture.py -b 'NUV' --skypos '[244.272750,55.268446]' -s 10 -f 'CRDra_NUV.csv' -a 0.002 -i 0.004 -o 0.006 -v 0 --maxgap 1000. --overwrite

python ./gAperture.py -b 'FUV' --skypos '[244.272750,55.268446]' -s 10 -f 'CRDra_FUV.csv' -a 0.002 -i 0.004 -o 0.006 -v 0 --maxgap 1000 --overwrite

####SDSS J084425
#Reported in [1].

python ./gAperture.py -b 'NUV' --skypos '[131.108178,51.642040]' -s 10 -f 'SDSSJ084425_NUV' -a 0.002 -i 0.004 -o 0.006 -v 0 --maxgap 1000 --overwrite

python ./gAperture.py -b 'FUV' --skypos '[131.108178,51.642040]' -s 10 -f 'SDSSJ084425_FUV' -a 0.002 -i 0.004 -o 0.006 -v 0 --maxgap 1000 --overwrite

####CDFS\_MOS00\_41
#Reported in [2] as an M dwarf with "stochastic variability." The classification was based upon variability per _epoch_. This light curve is a new result.

python ./gAperture.py -b 'NUV' --skypos '[53.3453,-27.3361]' -s 10 -v 2 -f CDFS_MOS00-41.csv --maxgap 1500 --trange '[975218628,975219882]' --overwrite

####Works Cited
[1] Welsh, Barry Y., et al. "GALEX high time-resolution ultraviolet observations of dMe flare events." arXiv preprint astro-ph/0608254 (2006).

[2] Gezari, S., et al. "The GALEX Time Domain Survey. I. Selection and Classification of Over a Thousand Ultraviolet Variable Sources." The Astrophysical Journal 766.1 (2013): 60.
