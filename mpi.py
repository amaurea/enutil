# This module defines parallel iteration over sequences.
from mpi4py import MPI

comm  = MPI.COMM_WORLD
myid  = comm.Get_rank()
nproc = comm.Get_size()

# Given a sequence which defines len and [], and an mpi
# communicator, define an object which can be iterated over,
# and which provides len, [] and has members specifying
# the mapping between local and global indices
class mine:
	def __init__(self, seq, comm=comm):
		self.seq   = seq
		self.comm  = comm
		self.myid  = comm.Get_rank()
		self.nproc = comm.Get_size()
		self.N     = len(seq)
		self.n     = (len(seq)-self.myid+self.nproc-1)/self.nproc
	def __len__(self):
		return self.n
	def __getitem__(self, i):
		class Entry: pass
		ind = Entry()
		ind.i = i
		ind.I = i*self.nproc+self.myid
		ind.n = self.n
		ind.N = self.N
		ind.nproc = self.nproc
		ind.myid  = self.myid
		ind.comm  = self.comm
		return ind, self.seq[ind.I]
	def __iter__(self):
		for i in range(len(self)):
			yield self[i]
