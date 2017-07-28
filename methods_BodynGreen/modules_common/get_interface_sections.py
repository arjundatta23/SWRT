#!/usr/bin/python

import sys
import gzip
import numpy as np
import matplotlib.pyplot as plt

# modules written by me
import get_slice as gs

class discontinuities():
	
	def __init__(self,modfile,xdcon,flipsides):

		if modfile.endswith('.gz'):
			fhandle = gzip.GzipFile(modfile,'r')
		else:
			fhandle = open(modfile,'r')
		#self.sgrid=int(float(fhandle.readline().split()[-2])/1000.0)
		self.sgrid=float(fhandle.readline().split()[-2])/1000.0
		xpts=[xdcon-self.sgrid,xdcon+self.sgrid]
		fhandle.close()
		del fhandle
		for i,x in enumerate(xpts):
			print "Extracting model at x = ", x
			vsobj=gs.vertical_slice(modfile,x)
			if i==0:
				vs1=np.array(vsobj.vs)/1000.0
				rho1=np.array(vsobj.rho)/1000.0
			else:
				vs2=np.array(vsobj.vs)/1000.0
				rho2=np.array(vsobj.rho)/1000.0
			if len(vsobj.vs) != vsobj.points_z:
				sys.exit('Error reading model')
		vsdiff=np.abs(vs2-vs1)
		if np.unique(vsdiff)[0]==0:
			# means there is at least one layer (most likely the half-space) for which there is no lateral discontinuity
			uvals=np.unique(vsdiff)[1:] 
		else:
			uvals=np.unique(vsdiff)
		betalr=np.zeros((len(uvals),2))
		rholr=np.zeros((len(uvals),2))
		self.veld=[] # veld is Vertical Extent of Lateral Discontinuity
		#print uvals
		for j in range(betalr.shape[0]):
			ind=np.where(vsdiff==uvals[j])[0]
			self.veld.append(ind)
			thk=len(ind)*self.sgrid
			print thk, vs1[ind[0]], vs2[ind[0]]
			betalr[j,0]=vs1[ind[0]]
			betalr[j,1]=vs2[ind[0]]
			rholr[j,0]=rho1[ind[0]]
			rholr[j,1]=rho2[ind[0]]
		dep=range(len(vs1))
		self.veld=np.array(self.veld)
		#plt.plot(vs1,dep)
		#plt.plot(vs2,dep)
		#plt.plot(vsdiff,dep)
		#plt.xlim(0,6)
		#plt.ylim(100,dep[0])
		#plt.show()
		if flipsides:
			self.beta_it=np.fliplr(betalr)
			self.rho_it=np.fliplr(rholr)
			self.get_horizontal_interfaces(vs2,vs1)
		else:
			self.beta_it=betalr
			self.rho_it=rholr
			self.get_horizontal_interfaces(vs1,vs2)

	def get_horizontal_interfaces(self,vsi,vst):
		
		uvi,indl=np.unique(vsi,return_index=True)
		uvt,indr=np.unique(vst,return_index=True)
		""" NB: by default numpy unique returns SORTED unique values. This can mess things up if layer velocities
		are not all increasing with depth (eg. low velocity layer at depth). Hence before proceeding we must
		rearrange uvi and uvt to represent actual ordering of layers """  
		#print "unique values as returned by numpy:	", uvi
		uvi=[b for a,b in sorted(zip(indl,uvi))]
		uvt=[b for a,b in sorted(zip(indr,uvt))]
		#print "unique values in correct order:	", uvi
		thks=[]
		for j,ui in enumerate(uvi):
			ind=np.where(vsi==ui)[0]
			if j==0:
				prev=0
			else:
				prev=thks[-1]
			thks.append(prev+(len(ind)*self.sgrid))
		self.ishif=thks[:-1] # lshif is Incidence Side Horizontal InterFaces
		thks=[]
		for j,ut in enumerate(uvt):
			ind=np.where(vst==ut)[0]
			if j==0:
				prev=0
			else:
				prev=thks[-1]
			thks.append(prev+(len(ind)*self.sgrid))
		self.tshif=thks[:-1] # Transmission Side Horizontal Interfaces

if __name__=='__main__':
	infile=sys.argv[1]
	discon=float(raw_input("x-coordinate of vertical discontinuity in model: "))
	ichoice=raw_input("Consider incidence from 'right hand side' ? (y/n): ")
	if ichoice=='y':
		flip=True
	else:
		flip=False
	dconobj=discontinuities(infile,discon,flip)
	print "vs on either side: ", dconobj.beta_it
	print "rho on either side: ", dconobj.rho_it
	print "Incidence side horizontal interfaces: ", dconobj.ishif
	print "Transmission side horizontal interfaces: ", dconobj.tshif
	print "vertical indices of lateral discontinuity: ", dconobj.veld
