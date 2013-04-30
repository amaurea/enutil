from pyslalib.slalib import *
import iers, pyfsla, math, numpy as np

arcsec = np.pi/180/60/60

# Transformations work in radians, with all angles as right-handed
# north polar coordinates (phi,theta), with theta=0 pointing in the
# z direction. Use astro2std to transform between these angles and
# the standard ones.
hor = 1; equ = 2; gal = 3

def astro2std(point, sys):
	res = np.array(point)
	off = int(res.shape[1]==3)
	p = res[:,off:]
	if sys == hor:
		p[:,0] = -p[:,0]
		p[:,1] = np.pi/2-p[:,1]
	elif sys == equ or sys == gal:
		p[:,1] = np.pi/2-p[:,1]
	return res

std2astro = astro2std

# Transform from apparent horizontal to equatorial (celestial) mean coordinates.
class hor2equ:
	def __init__(self, site):
		self.site = site
	def __call__(self, icoord):
		tcoord = np.array(icoord)
		mjd0 = icoord[0,0]
		info   = iers.lookup(mjd0)
		as2rad = math.pi/180/60/60
		ao = sla_aoppa(mjd0, info.dUT, self.site.lon, self.site.lat, self.site.alt,
			info.pmx*as2rad, info.pmy*as2rad, self.site.T, self.site.P, self.site.hum,
			299792.458/self.site.freq, 0.0065)
		am = sla_mappa(2000.0, mjd0)
		tcoord[:,1:] = pyfsla.aomulti(icoord[:,0], icoord[:,1:].T, ao, am).T
		return tcoord

class equ2gal:
	def __call__(self, icoord):
		tcoord = np.array(icoord)
		tcoord[:,1:] = pyfsla.equ2gal(icoord[:,1:].T).T
		return tcoord

class rotchain:
	def __init__(self, rots):
		self.rots = rots
	def __call__(self, icoord):
		out = np.array(icoord)
		for rot in self.rots:
			out = rot(out)
		return out

# This function performs rotation "rotation" on on coordinates icoord,
# but also computes the psi of the rotation this
# induces in the local coordinate system. This is done by duplicating
# every point into two closely offset points, transforming both, and
# measuring their relative orientation after the transformation.
# This is general (i.e. it will work independently of which underlying
# pointing library is used), but is a bit inefficient, and may suffer from
# numerical instability if the points are chosen to be too close.
class rot_polang:
	def __init__(self, rot, step=arcsec):
		self.rot  = rot
		self.step = step
	def __call__(self, icoord):
		icoord2 = np.array(icoord)
		icoord2[:,2] = icoord[:,2] - self.step
		# Handle pole-crossing
		inds = np.where(icoord2[:,2] < 0)
		icoord2[inds,2] = -icoord2[inds,2]
		icoord2[inds,1] = icoord2[inds,1]+np.pi
		# Transform both sets of coordinates
		ocoord  = self.rot(icoord)
		ocoord2 = self.rot(icoord2)
		# Make sure we don't have a 2pi skip between ocoord and ocoord2
		ocoord2 = ocoord + (ocoord2-ocoord+np.pi) % (2*np.pi) - np.pi
		# Get distance between pairs of points by haversine. These
		# should be very close to step.
		dists = haverdist(ocoord[:,1:],ocoord2[:,1:])
		# Calculate psi by the haversine formula
		hC  = (haversin(ocoord2[:,2])-haversin(ocoord[:,2]-dists))/(np.sin(ocoord[:,2])*np.sin(dists))
		# Numerical issues with very small angles may make hC negative:
		hC = np.maximum(0, hC)
		psi = 2*np.arcsin(hC**0.5)
		# Differentiate between left and right rotation:
		inds = np.where(ocoord2[:,1]>ocoord[:,1])
		psi[inds] = -psi[inds]
		res = np.zeros((icoord.shape[0],4))
		res[:,:3] = ocoord
		res[:,3]  = psi
		return res

def haversin(theta): return np.sin(theta/2)**2
def haverdist(c1, c2):
	lat1, lat2 = np.pi/2-c1[:,1], np.pi/2-c2[:,1]
	h = haversin(lat1-lat2) + np.cos(lat1)*np.cos(lat2)*haversin(c1[:,0]-c2[:,0])
	return 2*np.arcsin(h**0.5)

# Given a list of angles in radians, make sure the
# angles all change continuously, i.e. remove jumps due
# to [0:2pi]-bracketing.
def safe_angles(a):
	# First make sure all values are 0:2pi bracketed
	b = a % (2*np.pi)
	# Then remove all the jumps
	b[1:] -= np.cumsum(np.round((b[1:]-b[:-1])/(2*np.pi)))*2*np.pi
	return b

def safe_grid(a):
	# Make angles safe in each direction
	b = a % (2*np.pi)
	for d in range(a.ndim):
		c = np.swapaxes(b, d, -1)
		cflat = np.reshape(c, (np.prod(c.shape[:-1]), c.shape[-1]))
		for sub in cflat:
			sub[1:] -= np.cumsum(np.round((sub[1:]-sub[:-1])/(2*np.pi)))*2*np.pi
	b = b-np.round(np.mean(b)/(2*np.pi))*2*np.pi
	return b

# Storing full pointing is probably too much (i.e. too much time
# spent on I/O if re-reading, and not enough memory to hold it all
# in memoru if not. During map-making P size is Nscan*Nsamp*ndet*4 for T
# and Nscan*Nsamp*ndet*(4+4*3) for TQU. For comparison, the noise decorrelation
# step requires Nsamp*ndet*4*2 Hence, P dominates. By using
# interpolation, we can reduce P to Nscan*nsamp*4*2. As long as Nscan >> 1, this
# is a huge improvement.
#
# When using P, what we need is P(ndet,nsamp,[T,Q,U]). We can either store
# detectors as [daz,del,dpsi] and then use psi in the interpolation, followed
# by cos(2psi), sin(2psi), or we can store [daz,del],[T,Q,U], and interpolate
# in [Q,U].
#
# cos(a+b) = cos(a)*cos(b)-sin(a)*sin(b)
# sin(a+b) = cos(a)*sin(b)+sin(a)*cos(b)
#
# [cos(a+b)] = [ cos(b) -sin(b) ] * [cos(a)]
# [sin(a+b)]   [ sin(b)  cos(b) ]   [sin(a)]
#
# In our case, a = 2 dpsi and b = psi, so
# [T]'   [ 1      0          0     ]   [T]
# [Q]  = [ 0  cos(2psi) -sin(2psi) ] * [Q]
# [U]    [ 0  sin(2psi)  cos(2psi) ]   [U]
# This approach completely avoids trigonometric evaluations in the map-maker.
# The job of the interpolator is then to go from (mjd,az,el) -> (ra,dec,cos(2psi),sin(2psi)).
#
# Sievers' tiled interpolation is a good idea for constant elevation scans, so
# I'll implement that. The interpolation parameters will be stored in the level2-
# file, so the map-maker need only be concerned with using the interpolation,
# not creating it. In particular, the map-maker will not need to call slalib
# or other coordinate transformation libraries.
