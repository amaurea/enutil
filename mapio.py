import numpy as np, pyfits, sys, h5py as h5, warnings
warnings.filterwarnings('ignore')

# Write an (...,ndec,nra) map in the flat-sky equatorial
# CEA projection with limits given by box [[ra1,ra2],[dec1,dec2]]
# in radians
def write_map(fname, data, box, fmt=None):
	if fmt == None:
		if   fname[-5:] == ".fits": fmt = "fits"
		elif fname[-4:] == ".hdf":  fmt = "hdf"
	if fmt == "fits":
		write_fits(fname, data, box)
	elif fmt == "hdf":
		write_hdf(fname, data, box)
	else:
		raise ValueError

def read_map(fname, fmt=None):
	if fmt == None:
		if   fname[-5:] == ".fits": fmt = "fits"
		elif fname[-4:] == ".hdf":  fmt = "hdf"
	if fmt == "fits":
		return read_fits(fname)
	elif fmt == "hdf":
		return read_hdf(fname)
	else:
		raise ValueError

# This one implements the general World Coordinate System
# standard for FITS.
def write_fits(fname, data, box):
	ra, dec = np.array(box.T)*180/np.pi
	naxis = data.shape[::-1]
	#if ra[1] - ra[0] > 180: ra = np.array([ra[1]-360,ra[0]])
	cardList = pyfits.CardList()
	cardList.append(pyfits.Card('NAXIS', data.ndim))

	# Insert the pixel axes first
	n, c, delta = naxis[0], naxis[0]/2, (ra[1]-ra[0])/(naxis[0]-1)
	cardList.append(pyfits.Card('NAXIS1', n))
	cardList.append(pyfits.Card('CRPIX1', c+1))
	cardList.append(pyfits.Card('CDELT1', delta))
	cardList.append(pyfits.Card('CRVAL1', ra[0]+delta*c))
	cardList.append(pyfits.Card('CTYPE1', 'RA'))
	cardList.append(pyfits.Card('CUNIT1', 'DEG'))
	
	n, c, delta = naxis[1], naxis[1]/2, (dec[1]-dec[0])/(naxis[1]-1)
	cardList.append(pyfits.Card('NAXIS2', n))
	cardList.append(pyfits.Card('CRPIX2', c+1))
	cardList.append(pyfits.Card('CDELT2', delta))
	cardList.append(pyfits.Card('CRVAL2', dec[0]+delta*c))
	cardList.append(pyfits.Card('CTYPE2', 'DEC'))
	cardList.append(pyfits.Card('CUNIT2', 'DEG'))

	# And then all the remaining dimensions, which we don't know
	# much about.
	for i in range(2,data.ndim):
		cardList.append(pyfits.Card('NAXIS%d' % (i+1), naxis[i]))
		cardList.append(pyfits.Card('CRPIX%d' % (i+1), naxis[i]/2+1))

	header = pyfits.Header(cards=cardList)

	# And finally write the data. PyFits will write it in row-major order,
	# which fits interprets as a transposition, which is what we want to
	# make the above keywords match.
	hdus = pyfits.HDUList([pyfits.PrimaryHDU(data, header)])
	hdus.writeto(fname, clobber=True)

def wcs2range(n, ref, step, center):
	return [center - (ref-1)*step, center + (n-ref+1)*step]

# This function assums that we are using flat equatorial coordinates.
# It can be generalized later if needed.
def read_fits(fname):
	hdu    = pyfits.open(fname)[0]
	header = hdu.header
	box    = np.array([wcs2range(header["NAXIS%d" % i], header["CRPIX%d" % i], header["CDELT%d" % i],
		header["CRVAL%d" % i]) for i in range(1,3)]).T*np.pi/180
	print box.shape, box[0], box[1]
	return np.array(hdu.data), box

# Sadly, I don't know of any WCS equivalent for HDF, so this
# one will just, dump the data and box.
def write_hdf(fname, data, box):
	hfile = h5.File(fname,"w")
	hfile["data"] = data
	hfile["box"]  = box[:,::-1]*180/np.pi
	hfile["system"] = "equ"
	hfile.close()

def read_hdf(fname):
	with h5.File(fname,"r") as hfile:
		return np.array(hfile["data"]), np.array(hfile["box"])[:,::-1]*np.pi/180
