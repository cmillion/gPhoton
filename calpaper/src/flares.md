###Known Flares
This is mostly a list of known GALEX flares from the literature that we're able to reproduce with gAperture, along with the commands to do so.

One of the things that seems to happen in the literature is that the source positions get rounded off enough that they aren't really useful for positioning the aperture anymore, and we needed to go back to the MCAT to find the most likely nearby source. RA and Dec should really be reported with 5 or 6 significant figures (in decimal degrees).

####AF Psc
Reported in [1].

`./gAperture.py -b 'FUV' --skypos '[352.937473,-2.744136]' -s 10 -f 'AFPsc_FUV.csv' -a 0.002 -i 0.004 -o 0.006 -v 2 --maxgap 1000`

`./gAperture.py -b 'NUV' --skypos '[352.937473,-2.744136]' -s 10 -f 'AFPsc_NUV.csv' -a 0.002 -i 0.004 -o 0.006 -v 2 --maxgap 1000`

####CR Dra
Reported in [1]. Note: Only one flare was reported in the literature. The data appear to contain at least three and maybe four distinct flares.

`./gAperture.py -b 'NUV' --skypos '[244.272750,55.268446]' -s 10 -f 'CRDra_NUV.csv' -a 0.002 -i 0.004 -o 0.006 -v 2 --maxgap 1000.`

`./gAperture.py -b 'FUV' --skypos '[244.272750,55.268446]' -s 10 -f 'CRDra_FUV.csv' -a 0.002 -i 0.004 -o 0.006 -v 2 --maxgap 1000`

####SDSS J084425
Reported in [1].

`./gAperture.py -b 'NUV' --skypos '[131.108178,51.642040]' -s 10 -f 'SDSSJ084425_NUV' -a 0.002 -i 0.004 -o 0.006 -v 2 --maxgap 1000`

`./gAperture.py -b 'FUV' --skypos '[131.108178,51.642040]' -s 10 -f 'SDSSJ084425_FUV' -a 0.002 -i 0.004 -o 0.006 -v 2 --maxgap 1000`

####Works Cited
[1] Welsh, Barry Y., et al. "GALEX high time-resolution ultraviolet observations of dMe flare events." arXiv preprint astro-ph/0608254 (2006).

