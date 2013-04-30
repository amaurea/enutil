# Implementation of preconditioned conjugate gradients. More general than
# the scipy version, in that it does not assume that it knows how to perform
# the dot product. This makes it possible to use this for distributed x
# vectors. It is pretty much a straight port of the fortran cg solver in
# the quiet pipeline.

# According to
# http://etna.math.kent.edu/vol.13.2002/pp56-80.dir/pp56-80.pdf
# the error (As the A-norm of the difference between the current
# and true value) at step i can be approximated well as the sqrt of
# sum(j=0:d-1) alpha(i+j)*rz(i+j) for large enough d, where 4 is
# sufficient. This means that we can't estimtate the current error,
# but we can estimate the error 4 steps in the past, which is good
# enough.


import numpy as np

def default_M(x):     return np.copy(x)
def default_dot(a,b): return a.dot(b)

class CG:
	def __init__(self, A, b, x0=None, M=default_M, dot=default_dot):
		# Init parameters
		self.A   = A
		self.b   = b
		self.M   = M
		self.dot = dot
		self.x   = np.copy(x0)
		if self.x == None:
			self.x = np.zeros(b.shape)
		# Internal work variables
		n = b.size
		self.r   = b-self.A(self.x)
		self.z   = self.M(self.r)
		self.rz  = self.dot(self.r, self.z)
		self.rz0 = float(self.rz)
		self.p   = self.z
		self.err = np.inf
		self.d   = 4
		self.arz = []
		self.err_true = np.inf
		self.i   = 0
	def step(self):
		Ap = self.A(self.p)
		alpha = self.rz/self.dot(self.p, Ap)
		self.x += alpha*self.p
		self.r -= alpha*Ap
		self.z = self.M(self.r)
		next_rz = self.dot(self.r, self.z)
		self.err = next_rz/self.rz0
		beta = next_rz/self.rz
		self.rz = next_rz
		self.p = self.z + beta*self.p
		self.arz.append(self.rz*alpha)
		# Update proper error
		if len(self.arz) > self.d:
			# Good estimate of error d steps ago
			self.err_true = sum(self.arz[-self.d:])
		self.i += 1

def test():
	def A(x): return np.array([[4,1],[1,3]],dtype=float).dot(x)
	b = np.array([1.,2])
	cg = CG(A, b, x0=np.array([2.,1.]))
	while cg.err > 1e-4:
		cg.step()
		print cg.i, cg.err, cg.x
