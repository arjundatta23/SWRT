#!/usr/bin/python

import sys
import itertools
import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt

import read_earth_io as reo

#################################################

class ortho():
	
	def __init__(self,infile,ps,dcon):

		reoobj=reo.read_egnfile_per(infile,ps)
		omega=2*np.pi/ps
		self.orig_b1=reoobj.uzmat
		self.orig_b2=reoobj.urmat
		self.orig_b3=reoobj.tzmat
		self.orig_b4=reoobj.trmat
		self.origdep=reoobj.dep
		self.dcon=dcon
		nmodes=self.orig_b1.shape[1]
		matint=np.zeros((nmodes,nmodes)) # matrix of integrals
		self.norms=np.zeros(nmodes) # vector of integrals
		self.kmode=reoobj.wavnum.reshape(1,len(reoobj.wavnum))
		self.mu=reoobj.mu.reshape(len(reoobj.mu),1)
		self.lamda=reoobj.lamda.reshape(len(reoobj.mu),1)
		#print "shape is ", b1.shape
		#print "shapes are ", self.mu.shape, self.kmode.shape
		kmu=np.dot(self.mu,self.kmode)
		klamda=np.dot(self.lamda,self.kmode)
		d_b2_dz=(omega*self.orig_b4-np.multiply(kmu,self.orig_b1))/self.mu # numpy.multiply does element wise array multiplication
		d_b1_dz=(np.multiply(klamda,self.orig_b2)+omega*self.orig_b3)/(self.lamda+2*self.mu)
		dxz=np.gradient(self.orig_b2[:,0])
		dzz=np.gradient(self.orig_b1[:,0])
		if __name__=='__main__':
#		print "new shape is ", kmu.shape, d_b2_dz.shape
			#plt.plot(dxz,self.origdep,label='numerical')
			#plt.plot(d_b2_dz[:,0],self.origdep,label='theoretical')
			plt.plot(dzz,self.origdep,label='numerical')
			plt.plot(d_b1_dz[:,0],self.origdep,label='theoretical')
			plt.legend(loc='best')
			plt.ylim(100,0)
			plt.show()
			

		#**********************************************************
                # Before we integrate we need to add extra depth points to
                # the arrays in order to account for discontinuities in
                # the integration
		
		def modify_egn(ind,egn_in):
			# Get values of eigenfunctions AT discontinuity
			efdcon=egn_in[ind[0],:]
			#print "For interface ", dc
			#print "ind is ", ind
			#print "Value of eigenfunction at interface: ", efdcon
			# Allocate value to point just above discontinuity
			egn_out=np.insert(egn_in,ind[0],efdcon,axis=0)
			# and to point just below
			egn_out=np.insert(egn_out,ind[1]+1,efdcon,axis=0)
			return egn_out

		self.b1=self.orig_b1
		self.b2=self.orig_b2
		self.b4=self.orig_b4
		#print "shapes are: ", (klamda+2*kmu).shape, self.b2.shape, d_b1_dz.shape, self.lamda.shape
		self.psi_xx=((self.lamda*d_b1_dz)-(np.multiply((klamda+2*kmu),self.orig_b2)))/omega
                self.dep=self.origdep

		
		eps=0.000001
		#print "eigenfunction before adding extra point: ", self.b1[25:35,:]
		for dc in self.dcon:
                        valtoins=[dc-eps,dc+eps]
                        ins_ind=np.searchsorted(self.dep,valtoins)
			self.b1=modify_egn(ins_ind,self.b1)
			self.b2=modify_egn(ins_ind,self.b2)
			self.b4=modify_egn(ins_ind,self.b4)
			self.psi_xx=modify_egn(ins_ind,self.psi_xx)
			# Finally, add extra points to the depth array
                        self.dep=np.insert(self.dep,ins_ind,valtoins)
		#print "eigenfunction after adding extra point: ", self.b1[25:35,:]
		#print "Shapes of altered eigenfunctions: ", self.b1.shape, self.psi_xx.shape

		for ij in itertools.combinations_with_replacement(range(nmodes),2):
			i=ij[0]
			j=ij[1]
			ans=0.5*(self.integrate(self.b1[:,i],self.b1[:,j],self.b2[:,i],self.b2[:,j],self.b4[:,i],self.b4[:,j],self.psi_xx[:,i],self.psi_xx[:,j]))
			""" As opposed to the Love wave case, above expression ITSELF is equal to energy flux -- NO division by omega
				required -- see notes for details. """
			matint[j][i]=ans
			if __name__=='__main__':
				print "Mode number %d and %d: %f" %(i,j,ans)
			self.norms=np.diagonal(matint)

	def integrate(self,b1_m,b1_n,b2_m,b2_n,b4_m,b4_n,psi_m,psi_n):

		sumint=0.0
		for l in range(len(self.dcon)):
			sid=np.where(self.dep==self.dcon[l])[0][0]
			if l>0:
				sid_prev=np.where(self.dep==self.dcon[l-1])[0][0]
				top=sid_prev+1
			else:
				top=0
			#print "Depth sample excluded from integration: ",self.dep[sid]
			# integration above each horizontal interface
			integrand=b2_n[top:sid]*psi_m[top:sid] + b2_m[top:sid]*psi_n[top:sid] - b1_n[top:sid]*b4_m[top:sid] - b1_m[top:sid]*b4_n[top:sid]
			int_above=spi.simps(integrand,self.dep[top:sid])
                        sumint += int_above

                # integration below deepest horizontal interface
		integrand=b2_n[sid+1:]*psi_m[sid+1:] + b2_m[sid+1:]*psi_n[sid+1:] - b1_n[sid+1:]*b4_m[sid+1:] - b1_m[sid+1:]*b4_n[sid+1:]	
		int_below=spi.simps(integrand,self.dep[sid+1:])

		return abs(sumint + int_below)

if __name__=='__main__':
	effile=sys.argv[1] # eigenfunction file
	per=float(sys.argv[2]) # period in seconds
	dhif=raw_input("Enter depths of horizontal interfaces: ")
	discon=[float(d) for d in dhif.split()]
	orobj=ortho(effile,per,discon)
