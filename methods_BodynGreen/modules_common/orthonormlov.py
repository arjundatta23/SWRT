#!/usr/bin/python

import sys
import itertools
import numpy as np
import scipy.integrate as spi

import read_earth_io as reo

#################################################

class ortho():
	
	def __init__(self,infile,ps,dcon,usewt=True):

		reoobj=reo.read_egnfile_per(infile,ps)
		omega=2*np.pi/ps
		oegnfnmat=reoobj.utmat
		origdep=reoobj.dep
		self.dcon=dcon
		nmodes=oegnfnmat.shape[1]
		matint=np.zeros((nmodes,nmodes)) # matrix of integrals
		self.norms=np.zeros(nmodes) # vector of integrals
		self.kmode=reoobj.wavnum
		self.mu=reoobj.mu

		#**********************************************************
		# Before we integrate we need to add extra depth points to
		# the arrays in order to account for discontinuities in
		# the integration

		self.dep=origdep
		self.egnfnmat=oegnfnmat
	
		eps=0.000001
		for dc in self.dcon:
			valtoins=[dc-eps,dc+eps]
			ind=np.searchsorted(self.dep,valtoins)
			# Get values of eigenfunctions AT discontinuity
			efdcon=self.egnfnmat[ind[0],:]
			#print "For interface ", dc
			#print "ind is ", ind
			#print "eigenfunction before adding extra point: ", self.egnfnmat[ind[0]-5:ind[0]+5,:]
			#print "Value of eigenfunction at interface: ", efdcon
			# Allocate value to point just above discontinuity
			self.egnfnmat=np.insert(self.egnfnmat,ind[0],efdcon,axis=0)
			# and to point just below
			self.egnfnmat=np.insert(self.egnfnmat,ind[1]+1,efdcon,axis=0)
			# Add extra points to the depth and mu arrays
			self.dep=np.insert(self.dep,ind,valtoins)
			#print "LOOKING FOR ", ind[0]-1, ind[1], len(self.mu)
			if len(self.mu)>ind[1]:
				valtoins=[self.mu[ind[0]-1],self.mu[ind[1]]]
			else:
			#this happens when the very last depth point is itself an interface
				valtoins=[self.mu[ind[0]-1],self.mu[ind[1]-1]]
			self.mu=np.insert(self.mu,ind,valtoins)
			#print "eigenfunction after adding extra point: ", self.egnfnmat[ind[0]-5:ind[0]+5,:]
			#print "mu after extra: ", self.mu[ind[0]-5:ind[0]+5]
			#print "depth after extra: ", self.dep[ind[0]-5:ind[0]+5]
		# *********************************************************

		#for i in range(nmodes):
		for ij in itertools.combinations_with_replacement(range(nmodes),2):
                        i=ij[0]
                        j=ij[1]
			if not usewt:
				wt=np.ones(len(self.dep))
			else:
				wt=self.kmode[i]*self.mu
			ans=(self.integrate(self.egnfnmat[:,i],self.egnfnmat[:,j],wt))/omega
			""" Division by omega is necessary to get energy flux - see notes for details """
			if __name__=='__main__':
                                print "Mode number %d and %d: %f" %(i,j,ans)
			matint[j][i]=ans
			self.norms=np.diagonal(matint)

	def integrate(self,egf1,egf2,wfn):

		sumint=0.0
		for l in range(len(self.dcon)):
			sid=np.where(self.dep==self.dcon[l])[0][0]
			if l>0:
				sid_prev=np.where(self.dep==self.dcon[l-1])[0][0]
				top=sid_prev+1
			else:
				top=0
			if __name__=='__main__':
				print "Depth sample excluded from integration: ",self.dep[sid]
			# integration above each horizontal interface
			phi_ij=egf1[top:sid]*egf2[top:sid]
			prod=phi_ij*wfn[top:sid]
			int_above=spi.simps(prod,self.dep[top:sid])
			sumint += int_above

		# integration below deepest horizontal interface
		phi_ij=egf1[sid+1:]*egf2[sid+1:]
		prod=phi_ij*wfn[sid+1:]
		int_below=spi.simps(prod,self.dep[sid+1:])

		return sumint + int_below

if __name__=='__main__':
	effile=sys.argv[1] # eigenfunction file
	per=float(sys.argv[2]) # period in seconds
	dhif=raw_input("Enter depths of horizontal interfaces: ")
	discon=[float(d) for d in dhif.split()]
	orobj=ortho(effile,per,discon,True)
	print "Normalization factors are: ", orobj.norms
