import matplotlib.patches as mpatches
import matplotlib.pyplot as pyplot
import numpy
import pylab
import colorsys

"""Read in data."""
dbnum, mindec, maxdec, decrange, minzone, maxzone, numzones, minpartition, maxpartition, numpartitions = numpy.genfromtxt("db_info.txt", unpack=True, autostrip=True)

colors = [x/10. for x in range(10)]

this_figure = pyplot.figure(1)
pyplot.subplot(111, projection="aitoff")
pyplot.grid(True)
all_patches = []
for i in xrange(len(mindec)):
    pyplot.axhspan((numpy.pi/180.)*mindec[i], (numpy.pi/180.)*maxdec[i], color=colorsys.hsv_to_rgb(colors[i], 1.0, 1.0))
    all_patches.append(mpatches.Patch(color=colorsys.hsv_to_rgb(colors[i], 1.0, 1.0), label="DB {0:02.0F}: {1:+6.2F} < $\delta$ < {2:+6.2F}".format(i+1,mindec[i],maxdec[i])))

all_patches = all_patches[::-1]
pyplot.legend(handles=all_patches,bbox_to_anchor=(0.05,1.25,0.99,0.2), ncol=3, title="gPhoton Databases And Their DEC Ranges",fontsize="small", labelspacing=0.2, handletextpad=0.4,columnspacing=0.2)
this_figure.savefig("dbs_forpaper.eps", format="eps")
