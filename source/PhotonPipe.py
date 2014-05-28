import os
import csv
import time
from astropy.io import fits as pyfits
import numpy as np
from CalUtils import *
from FileUtils import *
from gnomonic import *
from MCUtils import *

def PhotonPipe(raw6file,scstfile,calpath,band,outbase,aspfile=0,ssdfile=0,nullfile=0,verbose=0,retries=20):

	startt = time.time()

	# Scale factor for the time column in the output csv so that it 
	#  can be recorded as an int in the database
	dbscale = 1000

	# This number determines the size of the chunks that gPhoton reads
	#  in from the raw6 for processing. Even if your machine has a lot
	#  of memory, making this number bigger is unlikely to improve the
	#  processing time much because so much is eaten up by the .csv write.
	chunksz = 1000000

	# These are just constants for the mission.
	detsize = 1.25 # Detector size in degrees
	pltscl = 68.754932 # Plate scale
	aspum = pltscl/1000.0
	ArcSecPerPixel = 1.5
	xi_xsc,xi_ysc,eta_xsc,eta_ysc = 0.,1.,1.,0.

	# Determine the eclipse number from the raw6 header.
	hdulist = pyfits.open(raw6file)
	hdr = hdulist[0].header
	hdulist.close()
	eclipse = hdr['eclipse']
	print "Eclipse is "+str(eclipse)+"."

	# Returns detector constants.
	print "Band is "+band+"."
	xclk, yclk, xcen, ycen, xscl, yscl, xslp, yslp = clk_cen_scl_slp(band,eclipse)

	# This determines the values for the post-CSP detector stim scaling
	#  and detector constant corrections.
	Mx,Bx,My,By,stimsep=1,0,1,0,0
	if (eclipse>37460):
		Mx,Bx,My,By,stimsep,yactbl=compute_stimstats(raw6file,band,eclipse)
		wig2, wig2data, wlk2, wlk2data, clk2, clk2data = postCSP_caldata(calpath)

	print "Loading wiggle files..."
	wiggle_x = get_fits_data(wiggle_filenames(band,calpath)['x'])
	wiggle_y = get_fits_data(wiggle_filenames(band,calpath)['y'])

	print "Loading walk files..."
	walk_x = get_fits_data(walk_filenames(band,calpath)['x'])
	walk_y = get_fits_data(walk_filenames(band,calpath)['y'])

	print "Loading linearity files..."
	linearity_x = get_fits_data(linearity_filenames(band,calpath)['x'])
	linearity_y = get_fits_data(linearity_filenames(band,calpath)['y'])

	# This is for the post-CSP stim distortion corrections.
	print "Loading distortion files..."
	if (eclipse>37460):
		print " Using stim separation of :"+str(stimsep)
	distortion_x = get_fits_data(distortion_filenames(band,calpath,eclipse,stimsep)['x'])
	distortion_y = get_fits_data(distortion_filenames(band,calpath,eclipse,stimsep)['y'])
	disthead = get_fits_header(distortion_filenames(band,calpath,eclipse,stimsep)['x'])
	cube_x0,cube_dx,cube_y0,cube_dy,cube_d0,cube_dd,cube_nd,cube_nc,cube_nr = disthead['DC_X0'],disthead['DC_DX'],disthead['DC_Y0'],disthead['DC_DY'],disthead['DC_D0'],disthead['DC_DD'],disthead['NAXIS3'],disthead['NAXIS1'],disthead['NAXIS2']

	if band == 'FUV':
		xoffset,yoffset = find_FUV_offset(scstfile,calpath)
	else:
		xoffset,yoffset = 0.,0.

	if os.path.isfile(str(ssdfile)):
		print "SSD file provided: "+str(ssdfile)
		stim_coef0,stim_coef1 = get_stim_coefs(ssdfile)
	elif ssdfile:
		print "SSD file requested: "+str(ssdfile)
		stim_coef0,stim_coef1 = create_SSD(raw6file,band,eclipse,ssdfile)
	else:
		print "No SSD file provided or requested."
		stim_coef0,stim_coef1 = create_SSD(raw6file,band,eclipse)
	print "		stim_coef0, stim_coef1 = "+str(stim_coef0)+", "+str(stim_coef1)

	print "Loading mask file..."
	mask = get_fits_data(mask_filename(band,calpath))
	maskinfo = get_fits_header(mask_filename(band,calpath))
	npixx = mask.shape[0]
	npixy = mask.shape[1]
	pixsz = maskinfo['CDELT2']
	maskfill = detsize/(npixx*pixsz)

	print "Loading aspect data..."
	# If not aspect file is provided, attempt to query the aspect database at MAST/STScI.
	if aspfile:
		aspra, aspdec, asptwist, asptime, aspheader, aspflags = load_aspect(aspfile)
	else:
		aspra, aspdec, asptwist, asptime, aspheader, aspflags = web_query_aspect(eclipse,retries=retries)

	minasp, maxasp = min(asptime), max(asptime)
	trange = [minasp,maxasp]
	print "			trange= ( "+str(trange[0])+" , "+str(trange[1])+" )"
	ra0, dec0, roll0 = aspheader['RA'], aspheader['DEC'], aspheader['ROLL']
	print "			[avgRA, avgDEC, avgROLL] = ["+str(aspra.mean())+", "+str(aspdec.mean())+", "+str(asptwist.mean())+"]"

	# This projects the aspect solutions onto the MPS field centers.
	print "Computing aspect vectors..."
	xi_vec, eta_vec = gnomfwd_simple(aspra, aspdec, ra0, dec0, -asptwist, 1.0/36000.0, 0.)

	print "Loading raw6 file..."
	raw6hdulist = pyfits.open(raw6file,memmap=1)
	raw6htab = raw6hdulist[1].header
	nphots = raw6htab['NAXIS2']
	print "		"+str(nphots)+" events"
	cnt = 0

	outfile = outbase+'.csv'
	print "Preparing output file "+outfile
	spreadsheet = csv.writer(open(outfile, 'wb'), delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

	# If specified, dump lines with NULLS into a separate csv file
	if nullfile:
		nullfile = outbase+'_NULL.csv'
		print "Preparing output file "+nullfile
		NULLspreadsheet = csv.writer(open(nullfile, 'wb'), delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

	print ""

	for i in xrange(int(nphots/chunksz)+1):
		a = time.time()

		csvrows = []
		chunkbeg, chunkend = i*chunksz, (i+1)*chunksz
		if chunkend > nphots:
			chunkend = nphots

		chunkid = " "+str(i+1)+" of "+str(int(nphots/chunksz)+1)+": "
		print_inline(chunkid+"Unpacking raw6 data...")
		t = np.array(raw6hdulist[1].data.field('t')[chunkbeg:chunkend])

                phb1 = np.array(raw6hdulist[1].data.field('phb1')[chunkbeg:chunkend],dtype='int64')
                phb2 = np.array(raw6hdulist[1].data.field('phb2')[chunkbeg:chunkend],dtype='int64')
                phb3 = np.array(raw6hdulist[1].data.field('phb3')[chunkbeg:chunkend],dtype='int64')
                phb4 = np.array(raw6hdulist[1].data.field('phb4')[chunkbeg:chunkend],dtype='int64')
                phb5 = np.array(raw6hdulist[1].data.field('phb5')[chunkbeg:chunkend],dtype='int64')

		q =    ((phb4 & 3) << 3) + ((phb5 & 224) >> 5)
		xb =   phb1 >> 5
		xamc = np.array( ((phb1 & 31) << 7), dtype='int16' ) + np.array( ((phb2 & 254) >> 1), dtype='int16') - np.array( ((phb1 & 16) << 8), dtype='int16')
		yb =   ((phb2 & 1) << 2) + ((phb3 & 192) >> 6)
		yamc = np.array( ((phb3 & 63) << 6), dtype='int16') + np.array( ((phb4 & 252) >> 2), dtype='int16') - np.array( ((phb3 & 32) << 7), dtype='int16')
		xa = ((phb5 & 16) >> 4) + ((phb5 & 3) << 3) + ((phb5 & 12) >> 1)
		xraw0 = xb*xclk + xamc
		yraw0 = yb*yclk + yamc
		ya = np.array( ((((yraw0/(2*yclk) - xraw0/(2*xclk)) + 10)*32) + xa), dtype='int64') % 32
		xraw = xraw0 + np.array((((xa+7) % 32) - 16), dtype='int64') * xslp
		yraw = yraw0 + np.array((((ya+7) % 32) - 16), dtype='int64') * yslp

		x = (xraw - xcen)*xscl
		y = (yraw - ycen)*yscl

		if (eclipse>37460):
			x = Mx*x+Bx
			y = My*y+By
			yac = rtaph_yac(yactbl,ya,yb,yamc,eclipse)
			y = y-yac
			yac = rtaph_yac2(q,xb,yb,ya,y,calpath,aspum,wig2,wig2data,wlk2,wlk2data,clk2,clk2data)
			y = y + yac

		# This and other lines like it below are for the purpose
		#  of memory management.
		phb1,phb2,phb3,phb4,phb5,xb,xamc,yb,yamc,xraw0,yraw0,xraw,yraw=[],[],[],[],[],[],[],[],[],[],[],[],[]

		flags = np.zeros(len(t))

		print_inline(chunkid+"Applying wiggle correction...")
		x_as = x*aspum
		y_as = y*aspum
		fptrx = x_as/10. + 240.
		fptry = y_as/10. + 240.

		x_as,y_as=[],[]

		# This and other lines like it below are to verify that the
		#  event is still on the detector.
		cut = ((fptrx>0.) & (fptrx<479.) & (fptry>0.) & (fptry<479.) & (flags==0))
		flags[np.where(cut==False)[0]] = 8
		ix = np.where(cut==True)[0]

		blt = fptrx-np.array(fptrx,dtype='int64')
		blu = fptry-np.array(fptry,dtype='int64')
		wigx,wigy = np.zeros(len(t)),np.zeros(len(t))
		wigx[ix] = (1-blt[ix])*(wiggle_x[xa[ix],np.array(fptrx[ix],dtype='int64')]) + (blt[ix])*(wiggle_x[xa[ix],np.array(fptrx[ix],dtype='int64')+1])
		wigy[ix] = (1-blu[ix])*(wiggle_y[ya[ix],np.array(fptry[ix],dtype='int64')]) + (blu[ix])*(wiggle_y[ya[ix],np.array(fptry[ix],dtype='int64')+1])

		xdig = x + wigx/(10.*aspum)
		ydig = y + wigy/(10.*aspum)

		wigx,wigy=[],[]

		print_inline(chunkid+"Applying walk correction...")
		xdig_as = xdig*aspum
		ydig_as = ydig*aspum
		fptrx = xdig_as/10. + 240.
		fptry = ydig_as/10. + 240.

		xdig_as,ydig_as=[],[]

		cut = ((fptrx>0.) & (fptrx<479.) & (fptry>0.) & (fptry<479.) & (flags==0))
		flags[np.where(cut==False)[0]] = 9
		ix = np.where(cut==True)[0]

		cut[ix] = ( (walk_x[q[ix],np.array(fptry[ix],dtype='int64'),np.array(fptrx[ix],dtype='int64')] != -999) |
					(walk_x[q[ix],np.array(fptry[ix],dtype='int64'),np.array(fptrx[ix],dtype='int64')+1] != -999) |
					(walk_x[q[ix],np.array(fptry[ix],dtype='int64')+1,np.array(fptrx[ix],dtype='int64')] != -999) |
					(walk_x[q[ix],np.array(fptry[ix],dtype='int64')+1,np.array(fptrx[ix],dtype='int64')+1] != -999) |
					(walk_y[q[ix],np.array(fptry[ix],dtype='int64'),np.array(fptrx[ix],dtype='int64')] != -999) |
					(walk_y[q[ix],np.array(fptry[ix],dtype='int64'),np.array(fptrx[ix],dtype='int64')+1] != -999) |
					(walk_y[q[ix],np.array(fptry[ix],dtype='int64')+1,np.array(fptrx[ix],dtype='int64')] != -999) |
					(walk_y[q[ix],np.array(fptry[ix],dtype='int64')+1,np.array(fptrx[ix],dtype='int64')+1] != -999) )
		flags[np.where(cut==False)[0]] = 9
		ix = np.where(cut==True)[0]

		blt = fptrx-np.array(fptrx,dtype='int64')
		blu = fptry-np.array(fptry,dtype='int64')
		walkx,walky = np.zeros(len(t)),np.zeros(len(t))
		walkx[ix] = (1-blt[ix])*(1-blu[ix])*(walk_x[q[ix],np.array(fptry[ix],dtype='int64'),np.array(fptrx[ix],dtype='int64')]) + (blt[ix])*(1-blu[ix])*(walk_x[q[ix],np.array(fptry[ix],dtype='int64'),np.array(fptrx[ix],dtype='int64')+1]) + (1-blt[ix])*(blu[ix])*(walk_x[q[ix],np.array(fptry[ix],dtype='int64')+1,np.array(fptrx[ix],dtype='int64')]) + (blt[ix])*(blu[ix])*(walk_x[q[ix],np.array(fptry[ix],dtype='int64')+1,np.array(fptrx[ix],dtype='int64')+1])
		walky[ix] = (1-blt[ix])*(1-blu[ix])*(walk_y[q[ix],np.array(fptry[ix],dtype='int64'),np.array(fptrx[ix],dtype='int64')]) + (blt[ix])*(1-blu[ix])*(walk_y[q[ix],np.array(fptry[ix],dtype='int64'),np.array(fptrx[ix],dtype='int64')+1]) + (1-blt[ix])*(blu[ix])*(walk_y[q[ix],np.array(fptry[ix],dtype='int64')+1,np.array(fptrx[ix],dtype='int64')]) + (blt[ix])*(blu[ix])*(walk_y[q[ix],np.array(fptry[ix],dtype='int64')+1,np.array(fptrx[ix],dtype='int64')+1])

		print_inline(chunkid+"Applying linearity correction...")
		xp = xdig - walkx
		yp = ydig - walky
		xp_as = xp*aspum
		yp_as = yp*aspum
		fptrx = xp_as/10. + 240.
		fptry = yp_as/10. + 240.

		xp,yp=[],[]
		walkx,walky=[],[]

		cut = ((fptrx>0.) & (fptrx<479.) & (fptry>0.) & (fptry<479.) & (flags==0))
		flags[np.where(cut==False)[0]] = 10
		ix = np.where(cut==True)[0]

		blt = fptrx-np.array(fptrx,dtype='int64')
		blu = fptry-np.array(fptry,dtype='int64')

		dx,dy = np.zeros(len(t)),np.zeros(len(t))
		dx[ix] = (1-blt[ix])*(1-blu[ix])*linearity_x[np.array(fptry[ix],dtype='int64'),np.array(fptrx[ix],dtype='int64')] + (blt[ix])*(1-blu[ix])*linearity_x[np.array(fptry[ix],dtype='int64'),np.array(fptrx[ix],dtype='int64')+1] + (1-blt[ix])*(blu[ix])*linearity_x[np.array(fptry[ix],dtype='int64')+1,np.array(fptrx[ix],dtype='int64')] + (blt[ix])*(blu[ix])*linearity_x[np.array(fptry[ix],dtype='int64')+1,np.array(fptrx[ix],dtype='int64')+1]
		dy[ix] = (1-blt[ix])*(1-blu[ix])*linearity_y[np.array(fptry[ix],dtype='int64'),np.array(fptrx[ix],dtype='int64')] + (blt[ix])*(1-blu[ix])*linearity_y[np.array(fptry[ix],dtype='int64'),np.array(fptrx[ix],dtype='int64')+1] + (1-blt[ix])*(blu[ix])*linearity_y[np.array(fptry[ix],dtype='int64')+1,np.array(fptrx[ix],dtype='int64')] + (blt[ix])*(blu[ix])*linearity_y[np.array(fptry[ix],dtype='int64')+1,np.array(fptrx[ix],dtype='int64')+1]

		print_inline(chunkid+"Applying stim distortion correction...")

		ss = stim_coef0 + (t * stim_coef1) # stim separation
		col,row,depth = np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t))
		col[ix] = ( xp_as[ix] - cube_x0) / cube_dx
		row[ix]   = ( yp_as[ix] - cube_y0) / cube_dy
		depth[ix] = ( ss[ix]    - cube_d0) / cube_dd

		# FIXME: This throws an error sometimes like the following...
		# PhotonPipe.py:262: RuntimeWarning: invalid value encountered in less
		# depth[((depth < 0)).nonzero()[0]] = 0.
		# PhotonPipe.py:263: RuntimeWarning: invalid value encountered in greater_equal
		#  depth[((depth >= cube_nd)).nonzero()[0]] = -1.
		# ERROR: IndexError: index -9223372036854775808 is out of bounds for axis 0 with size 17 [PhotonPipe]
		depth[((depth < 0)).nonzero()[0]] = 0.
		depth[((depth >= cube_nd)).nonzero()[0]] = -1.

		cut = ((col>-1) & (col<cube_nc) & (row>-1) & (row<cube_nr) & (flags==0))
		flags[np.where(cut==False)[0]] = 11
		ix = np.where(cut==True)[0]

		xshift,yshift = np.zeros(len(t)),np.zeros(len(t))
		xshift[ix] = distortion_x[np.array(depth[ix],dtype='int64'),np.array(row[ix],dtype='int64'),np.array(col[ix],dtype='int64')]
		yshift[ix] = distortion_y[np.array(depth[ix],dtype='int64'),np.array(row[ix],dtype='int64'),np.array(col[ix],dtype='int64')]

		xshift = (xshift*ArcSecPerPixel)+xoffset
		yshift = (yshift*ArcSecPerPixel)+yoffset

		print_inline(chunkid+"Applying hotspot mask...")

		# The detectors aren't oriented the same way.
		flip = 1.
		if band=='FUV':
			flip = -1.

		xi  =  xi_xsc*(flip*(yp_as + dy + yshift)*10.) +  xi_ysc*(flip*(xp_as + dx + xshift)*10.)
		eta = eta_xsc*(flip*(yp_as + dy + yshift)*10.) + eta_ysc*(flip*(xp_as + dx + xshift)*10.)

		xp_as,yp_as,depth,ss=[],[],[],[]
		dx,dy,xshift,yshift=[],[],[],[]

		col = ((( xi/36000.)/(detsize/2.)*maskfill + 1.)/2. * npixx)
		row = (((eta/36000.)/(detsize/2.)*maskfill + 1.)/2. * npixy)

		cut = ((col>0.) & (col<799.) & (row>0.) & (row<799.) & (flags==0))
		flags[np.where(cut==False)[0]] = 6
		ix = np.where(cut==True)[0]

		cut[ix] = ((mask[np.array(col[ix],dtype='int64'),np.array(row[ix],dtype='int64')] == 1.))
		flags[np.where(cut==False)[0]] = 6
		ix = np.where(cut==True)[0]

		col,row=[],[]

		# This gives the index of the aspect time that comes _before_
		#  each photon time. Without the '-1' it will give the index
		#  of the aspect time _after_ the photon time.
		print_inline(chunkid+"Mapping photon times to aspect times...")
		aspix = np.digitize(t,asptime)-1

		print_inline(chunkid+"Applying dither correction...")
		# Use only photons that are bracketed by valid aspect solutions
		#  and have been not themselves been flagged as invalid.
                cut = ((aspix>0) & (aspix<(len(asptime)-1)) & ( (flags==0) | (flags==6) ))
		flags[np.where(cut==False)[0]] = 7
		ix = np.where(cut==True)[0]

		print_inline(chunkid+"Interpolating aspect solutions...")
		dxi,deta = np.zeros(len(t)),np.zeros(len(t))
		dxi[ix] = (xi_vec[aspix[ix]+1]-xi_vec[aspix[ix]])*(t[ix]-asptime[aspix[ix]])/(asptime[aspix[ix]+1]-asptime[aspix[ix]])
		deta[ix]= (eta_vec[aspix[ix]+1]-eta_vec[aspix[ix]])*(t[ix]-asptime[aspix[ix]])/(asptime[aspix[ix]+1]-asptime[aspix[ix]])

		print_inline(chunkid+"Mapping to sky...")
		ra,dec = np.zeros(len(t)),np.zeros(len(t))
		ra[ix], dec[ix] = gnomrev_simple(xi[ix]+dxi[ix],eta[ix]+deta[ix],aspra[aspix[ix]],aspdec[aspix[ix]],-asptwist[aspix[ix]],1/36000.,0.)

		cut =  ( ((asptime[aspix[ix]+1]-asptime[aspix[ix]])==1) & (aspflags[aspix[ix]]%2==0) & (aspflags[aspix[ix]+1]%2==0) & (aspflags[aspix[ix]-1]%2==0) & (flags[ix]==0) & (flags[ix]!=7))
		flags[np.where(cut==False)[0]] = 12

		# NOTE: If you wish to add a hook that filters the gPhoton
		#	output (like perhaps by sky position or time range)
		#	then add it here. I reccomend that you use the
		#	"ix = np.where" technique used above.
		# TODO: Preprogram a (commented out) filter on RA/Dec.

		print_inline(chunkid+"Writing to spreadsheet...")

		# The issue is that we need to recombine the data into rows to
		#  feed to csv.writerow, hence the loop.
		# It might be possible to do without a loop. One way would be
		#  with numpy.column_stack except that this requires stupid
		#  amounts of memory and therefore takes longer to run than the
		#  loop iffen it manages to complete at all without a memory
		#  error or segmentation fault.
		for i in xrange(len(t)):
			cnt+=1
			# To avoid repeat indexing of flags...
			thisflag = flags[i]
			# This substitutes empty strings for RA and Dec
			#  values so that when they're dumped into the database
			#  they are correctly recorded as NULL
			if (thisflag == 2) or (thisflag == 5) or (thisflag == 7) or (thisflag == 8) or (thisflag == 9) or (thisflag == 10) or (thisflag == 11) or (thisflag == 12):
				if nullfile:
					NULLspreadsheet.writerow([int(t[i]*dbscale),x[i],y[i],xa[i],ya[i],q[i],xi[i],eta[i],"","",flags[i]])
				else:
					spreadsheet.writerow([int(t[i]*dbscale),x[i],y[i],xa[i],ya[i],q[i],xi[i],eta[i],"","",flags[i]])
			else:
				spreadsheet.writerow([int(t[i]*dbscale),x[i],y[i],xa[i],ya[i],q[i],xi[i],eta[i],ra[i],dec[i],flags[i]])

	raw6hdulist.close()
	stopt = time.time()

	print_inline("")
	print ""
	print "Runtime statistics:"
	print "	runtime		=	"+str(stopt-startt)+" sec. = ("+str((stopt-startt)/60.)+" min.)"
	print "	processed	=	"+str(cnt)+" of "+str(nphots)+" events."
	if cnt<nphots:
		print "		WARNING: MISSING EVENTS!"
	print "	rate		=	"+str(nphots/(stopt-startt))+" photons/sec."
	print ""

	return

