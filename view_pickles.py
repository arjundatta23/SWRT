#!/usr/bin/python

import sys
import pickle
import numpy as np
import matplotlib.pyplot as plt

#######################################################

def read_single(pfile):
	jar=open(pfile)
	print "Reading ", pfile
	cookie1 = pickle.load(jar)
	cookie2 = pickle.load(jar)
	try:
		cookie3 = pickle.load(jar)
	except EOFError:
		cookie3 = None
	jar.close()
	#print cookie1
	#print cookie2
	#print len(cookie2)
	return cookie1, cookie2, cookie3

#######################################################

nfiles=len(sys.argv)
pfnames=[]
compare_pub=False
for pf in range(1,nfiles):
	if sys.argv[pf].endswith('.pckl'):
		pfnames.append(sys.argv[pf])
	else:
		datfile=sys.argv[pf]
		pub_data=open(datfile)
		compare_pub=True
usrc=int(raw_input("Plotting results for different models (1) or different methods on the same model (2) ? : "))
iseng=False
usrc2=raw_input("Is this/are these energy pickle(s)? (y/n): ")
if usrc2=='y':
	iseng=True
for i,pkl in enumerate(pfnames):
	freq,values,verr=read_single(pkl)
	cols=['k','b','r','g','y','m','0.3','0.5']
	mname=['Model F','Model L','r','g','y','m']
	#mname=['dz 0.5','dz 1','r','g','y','m']
	#mname=['rhs 7','rhs 5','r','g','y','m']
	#mname=['FD-pert','FD-Bielak']
	#mname=['Alsop method', 'Body wave method', 'Green''s function method']
	#mname=['1D model (incidence side medium)', '2D model (transmission side)']
	#mname=['Model 1, 'r'$\alpha=0.1$','Model 2, 'r'$\alpha=0.4$','Model 3, 'r'$\alpha=0.7$']
	if usrc==1:
		cname="Model %d" %(i+1)
	#print "Lengths are: ", len(freq), len(values[0]), len(values[1])
		for mode in range(len(values)):
		#for mode in range(1):
			cname_temp="Mode %d" %(mode)
			if mode==0:
			#if (i==0 and mode==0) or (i==1 and mode==4):	# USED FOR FIGURE 4.4 in thesis
				plt.plot(freq,values[mode],'-o',color=cols[i],label=mname[i])
			else:
				plt.plot(freq,values[mode],'--',color=cols[i])
				#plt.plot(freq,values[mode],'--',color=cols[i],label=cname_temp)
				#plt.plot(freq,values[mode],'--.',color=cols[i])
	elif usrc==2:
		for mode in range(len(values)):
		#for mode in range(1):
			cname="Mode %d" %(mode)
			if i==0:
				#if mode==1 or mode==2:
				#	values[mode]=[-1*j for j in values[mode]]
				try:
					plt.plot(freq,values[mode],color=cols[mode],label=cname)
				except IndexError:
					# means no. of modes in pickle is > len(colors)
					plt.plot(freq,values[mode],label=cname)
				#plt.errorbar(freq,values[mode],yerr=verr[mode],ecolor='gray',label=cname)
			elif i==1:
				if verr==None:
					try:
						plt.plot(freq,values[mode],'--o',color=cols[mode],mec=cols[mode],mfc="None",ms=8)
					except IndexError:
						plt.plot(freq,values[mode],'-o')
				else:
					try:
						#plt.errorbar(freq,values[mode],yerr=verr[mode],fmt='.',color=cols[mode],ecolor='gray')
						plt.plot(freq,values[mode],'.',color=cols[mode],markersize=6)
					except IndexError:
						# means no. of modes in pickle is > len(colors)
						plt.errorbar(freq,values[mode],yerr=verr[mode],fmt='.',ecolor='gray',label=cname)
			elif i==2:
				try:
					plt.plot(freq,values[mode],'-*',color=cols[mode],mec=cols[mode],ms=8)
				except IndexError:
					plt.plot(freq,values[mode],'-*')
				
		if iseng and len(values)>1:
			if i==0:
				plt.plot(freq,np.nansum(values,axis=0),'--',color='k',label='Total Analytic')
			elif i==1:
				plt.plot(freq,np.nansum(values,axis=0),'-.',color='k',label='Total FD')
#plt.ylim(-0.16,0.15)
#plt.ylim(0.4,1.2)
#plt.ylim(0.7,0.8)
#plt.ylim(1.2,1.28)
#plt.ylim(0,1)
#plt.xlim(0,0.05)
#plt.text(0.02,0.2,r'Model: $\alpha=0.7$',fontsize=18)
#plt.text(0.1,0.6,r'Model: $\alpha=0.4$',fontsize=18)
#plt.text(0.08,1.3,r'M-step model: $\alpha=0.4$',fontsize=18)
#plt.text(0.12,0.95,'Model L',fontsize=18)
#plt.text(0.01,0.65,'Rayleigh waves',fontsize=18)
#plt.text(0.08,0.2,'Love waves',fontsize=18)
if usrc==1:
	plt.legend(loc='best')
	#plt.legend(loc='best',prop={'size':12})
if usrc==2 and len(values)>1:
	plt.legend(loc='best')
	#plt.legend(loc=9,prop={'size':12},ncol=4)
	#plt.legend(loc=3)
plt.xlabel("Frequency [Hz]")
if iseng:
	plt.ylabel("Fraction of incident energy")
else:
	plt.ylabel("Transmission surface ratio")
	#plt.ylabel("Mode Participation Factors")
#plt.ylabel("Reflection Coefficient")
#plt.grid(True,color='0.6')
if compare_pub:
	x=[]
	y=[]
	for line in pub_data:
		x.append(float(line.split(',')[0]))
		y.append(float(line.split(',')[1]))
	plt.plot(x,y,'o')
	pub_data.close()
	del pub_data
plt.show()
