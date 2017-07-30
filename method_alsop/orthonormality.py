#!/usr/bin/python

import sys
import itertools
import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt

import read_earth_io as reo

#################################################

class ortho():
	
	def __init__(self,infile,ps,dcon,usewt):

		reoobj=reo.read_egnfile_per(infile,ps)
		self.egnfnmat=reoobj.utmat
		self.dep=reoobj.dep
		self.dcon=dcon
		nmodes=self.egnfnmat.shape[1]
		self.matint=np.zeros((nmodes,nmodes)) # matrix of integrals
		self.kmode=reoobj.wavnum
		self.mu=reoobj.mu

		#**********************************************************
		# Before we integrate we need to add extra depth points to
		# the arrays in order to account for the discontinuity in
		# the integration
		
		eps=0.000001
		valtoins=[self.dcon-eps,self.dcon+eps]
		ind=np.searchsorted(self.dep,valtoins)
		#print "ind is ", ind
		# Get values of eigenfunctions AT discontinuity
		efdcon=self.egnfnmat[ind[0],:]
		#print efdcon
		# Allocate value to point just above discontinuity
		self.egnfnmat=np.insert(self.egnfnmat,ind[0],efdcon,axis=0)
		# and to point just below
		self.egnfnmat=np.insert(self.egnfnmat,ind[1]+1,efdcon,axis=0)
		#print self.egnfnmat[ind[0]-5:ind[0]+5,:]
		# Add extra points to the depth and mu arrays
		self.dep=np.insert(self.dep,ind,valtoins)
		valtoins=[self.mu[ind[0]-1],self.mu[ind[1]]]
		self.mu=np.insert(self.mu,ind,valtoins)
		#print "mu is: ", self.mu[ind[0]-5:ind[0]+5]
		#print "depth is: ", self.dep[ind[0]-5:ind[0]+5]
		# *********************************************************

		for ij in itertools.combinations_with_replacement(range(nmodes),2):
			i=ij[0]
			j=ij[1]
			if not usewt:
				wt=np.ones(len(self.dep))
			else:
				if i==j:
					wt=self.kmode[i]*self.mu
				else:
					wt=self.mu
			ans=self.integrate(self.egnfnmat[:,i],self.egnfnmat[:,j],wt)
			#self.matint[ij]=ans
			self.matint[j][i]=ans
			if __name__=='__main__':
				print "Mode number %d and %d: %f" %(i,j,ans)
		if usewt: 
			self.norms=np.diagonal(self.matint)
		

	def integrate(self,egf1,egf2,wfn):

		sid=np.where(self.dep==self.dcon)[0][0]
		#print "sid is ", sid
		# integration above discontinuity
		phi_ij=egf1[:sid]*egf2[:sid]
		prod=phi_ij*wfn[:sid]
		int_above=spi.simps(prod,self.dep[:sid])

		# integration below discontinuity
		phi_ij=egf1[sid+1:]*egf2[sid+1:]
		prod=phi_ij*wfn[sid+1:]
		int_below=spi.simps(prod,self.dep[sid+1:])

		#plt.plot(egf1,self.dep)
		#plt.plot(egf2,self.dep)
		#plt.show()

		return int_above + int_below

if __name__=='__main__':
	effile=sys.argv[1] # eigenfunction file
	per=float(sys.argv[2]) # period in seconds
	discon=float(raw_input("Enter depth of discontinuity: "))
	orobj=ortho(effile,per,discon,True)
	print orobj.matint
	print "Normalization factors are: ", orobj.norms
