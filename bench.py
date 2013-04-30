# Simple benchmarking utility. Compared to the fortran version, this is
# not fast enough to be used in tight inner loops, but it will serve my
# purposes. It keeps a full history of the bencharmking events as a flat
# list.
import time, numpy as np

class Bench:
	def __init__(self):
		self.entries = []
	def start(self, name):
		self.entries.append([0, time.time(), time.clock(), name])
	def stop(self, name):
		self.entries.append([1, time.time(), time.clock(), name])
	def dump(self, fname, fmt="%d %12.5f %12.5f %s"):
		f = open(fname, "w")
		for e in self.entries:
			print >> f, fmt % tuple(e)
		f.close()
#
#class Entry:
#	def __init__(self):
#		self.t  = 0.0
#		self.dt = 0.0
#		self.n  = 0
#
#class Bench:
#	def __init__(self):
#		self.steps = {}
#	def start(self, name):
#		if not name in self.steps:
#			self.steps[name] = Entry()
#		self.steps[name].t = time.time()
#	def stop(self, name):
#		# Finish previous step first
#		step = self.steps[name]
#		step.dt += time.time() - step.t
#		step.n  += 1
#	def add(self, other):
#		for k, v in other.steps.iteritems():
#			if not k in self.steps:
#				self.steps[k] = v
#			else:
#				self.steps[k].dt += v.dt
#				self.steps[k].n  += v.n
#	# Reduce collects all information in the
#	# root node, leaving the others empty.
#	def reduce(self, comm):
#		myid  = comm.Get_rank()
#		nproc = comm.Get_size()
#		if myid == 0:
#			for i in range(1, nproc):
#				self.add(comm.recv(source=i))
#		else:
#			comm.send(self, dest=0)
#			self.steps = {}
#	def to_str(self, fmt="%-16s %8d %15.7e %15.7e\n"):
#		keys = sorted(list(self.steps), key=str.lower)
#		res  = ""
#		for k in keys:
#			v = self.steps[k]
#			res += fmt % (k, v.n, v.dt, v.n and v.dt/v.n or np.nan)
#		return res
#	def prt(self, fmt="%-16s %8d %15.7e %15.7e\n"):
#		print self.to_str(fmt),
#	def dump(self, fname, fmt="%-16s %8d %15.7e %15.7e\n"):
#		f = open(fname,"w")
#		print >> f, self.to_str(fmt),
#		f.close()
#
#def btest(comm):
#	bench = Bench()
#	for j in range(10):
#		bench.start("A")
#		for i in range(1000): pass
#		bench.stop("A")
#		bench.start("B")
#		for i in range(10000): pass
#		bench.stop("B")
#		bench.start("C")
#		for i in range(100000): pass
#		bench.stop("C")
#		bench.start("D")
#		for i in range(1000000): pass
#		bench.stop("D")
#	bench.reduce(comm)
#	myid = comm.Get_rank()
#	if myid == 0:
#		bench.prt()
