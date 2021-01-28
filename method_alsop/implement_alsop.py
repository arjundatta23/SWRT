#!/usr/bin/python

# Available python modules
import os
import sys
import pickle
import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt

sys.path.append(os.path.expanduser('~/code_general/modules.python'))
# path to the "SW1D_earthsr" set of modules

# modules writted by me
import maineq as meq
import SW1D_earthsr.utils_pre_code as upre
import SW1D_earthsr.egnfunc_integrals_norms as ein

###################################################################################################

def do_single_freq(per):

	#****************************************************
	# Get normalization factors
	#****************************************************
	print("Computing normalization factors for period %f" %(per))
	omega=2*np.pi/per

	# medium 1
	oobj1 = ein.energy_integrals_lov(effile1, per, dep_pts_mod1, dcon1)
	dep1=oobj1.dep # depth sampling for medium 1
	mu1=oobj1.mu
	k1=oobj1.kmode
	efmat1=oobj1.l1 # matrix containing eigenfunctions for medium 1
	oobj1.orthogonality_products()
	Nmed1=oobj1.norms
	n=len(Nmed1) # no. of modes in medium 1 (at period per)
	Nmat1=np.sqrt(np.outer(Nmed1,Nmed1)) # matrix of normalizing factors for medium 1

	# medium 2
	oobj2 = ein.energy_integrals_lov(effile2, per, dep_pts_mod2, dcon2)
	dep2=oobj2.dep # depth sampling for medium 2
	# different from dep1 in general because additional depth points will have been added by the ort module
	# around the horizontal interfaces of the two media
	efmat2=oobj2.l1 # matrix containing eigenfunctions for medium 2
	oobj2.orthogonality_products()
	Nmed2=oobj2.norms
	m=len(Nmed2) # no. of modes in medium 2
	Nmat2=np.sqrt(np.outer(Nmed2,Nmed2)) # matrix of normalizing factors for medium 2

	# matrix of mixed (medium 1 & 2) normalizing factors - required for the integrals T and P
	N12=np.sqrt(np.outer(Nmed1,Nmed2))

	# ******************************************************************************************************
	# Compute required integrals - build P,S,T,V matrices
	# S and V matrices can be built by ort module as only 1 medium is involved hence the eigenfunctions have
	# same depth sampling.
	# P and T matrices must be built from scratch in this program as two different media are involved,
	# so depth sampling of eigenfunctions is different (in general), so ort cannot do the required integrations
	# *******************************************************************************************************
	print("Computing required integrals")
	S=np.zeros((n,n))
	V=np.zeros((m,m))
	T=np.zeros((n,m))
	P=np.zeros((n,m))

	# S matrix
	sobj=ein.energy_integrals_lov(effile1, per, dep_pts_mod1, dcon1, False)
	sobj.orthogonality_products()
	S=sobj.matint/Nmat1

	# V matrix
	vobj=ein.energy_integrals_lov(effile2, per, dep_pts_mod2, dcon2, False)
	vobj.orthogonality_products()
	V=vobj.matint/Nmat2

	# T and P matrices
	def integrate(ef1,ef2,z1,z2,wt=None):
		checking=False
		if wt is None:
			weight=np.ones(len(z1))
		else:
			weight=wt
		if len(np.setdiff1d(z1,z2))>0:
			# print(np.setdiff1d(z1,z2))
			# print(np.setdiff1d(z2,z1))
			""" Strategy will be to modify depth sampling of eigenfunction of medium 2 - remove extra points at discontinuity of medium 2 and add extra points at discontinuity of medium 1. This is because in the T and P integrals, only mu1 has a role, not mu2 """
			to_add=np.setdiff1d(z1,z2)
			to_remove=np.setdiff1d(z2,z1)
			addid=np.where(z1==to_add[0])[0][0]
			remid=np.where(z2==to_remove[0])[0][0]
			# remove appropriate points
			if checking:
				print("z2 before removal: ", z2[remid-5:remid+5])
				print("ef2 before removal: ", ef2[remid-5:remid+5])
				print("deleted samples ", z2[remid], z2[remid+2])
			z2=np.delete(z2,[remid,remid+2])
			ef2=np.delete(ef2,[remid,remid+2])
			if checking:
				print("z2 after removal ", z2[remid-5:remid+5])
				print(" ef2 after removal ", ef2[remid-5:remid+5])
			# add appropriate points
			z2=np.insert(z2,[addid,addid+1],[z1[addid],z1[addid+2]]) # careful with usage
			if checking:
				print("ef2 before additon: ", ef2[addid-5:addid+5])
			ef2=np.insert(ef2,[addid,addid+1],[ef2[addid],ef2[addid]])   # of numpy insert
			if checking:
				print("added samples ", z1[addid], z1[addid+2])
				print("z2 after addition ", z2[addid-5:addid+5])
				print("ef2 after addition ", ef2[addid-5:addid+5])
			if len(z2) != len(z1) or len(ef2) != len(ef1):
				sys.exit('Problem with depth sampling of eigenfunctions from the two media')

		sumint=0.0
		for l in range(len(dcon1)):
			sid=np.where(z1==dcon1[l])[0][0]
			if l>0:
				sid_prev=np.where(z1==dcon1[l-1])[0][0]
				top=sid_prev+1
			else:
				top=0

			# integration above each discontinuity
			phi_ij=ef1[top:sid]*ef2[top:sid]
			prod=phi_ij*weight[top:sid]
			int_above=spi.simps(prod,z1[top:sid])
			sumint += int_above

		# integration below deepest discontinuity
		phi_ij=ef1[sid+1:]*ef2[sid+1:]
		prod=phi_ij*weight[sid+1:]
		int_below=spi.simps(prod,z1[sid+1:])

		integral = sumint + int_below
		return integral

	for i in range(n):
		for j in range(m):
			wfn=k1[i]*mu1
			T[i,j]=integrate(efmat1[:,i],efmat2[:,j],dep1,dep2)/omega
			P[i,j]=integrate(efmat1[:,i],efmat2[:,j],dep1,dep2,wfn)/omega
			# NB: in the above, division by omega is not part of the written equations but
			# it is required for consistency with the ort module -- in the Love wave ort module,
			# I have divided the normalising integral by omega so that it equals the energy
			# flux (see notes).
	T=T/N12
	P=P/N12


	########################################################################
	# Set up and solve matrix equation
	########################################################################

	rc,tc = meq.do_main(P,S,T,V)
	if len(allper)==1:
		print("Matrix of normalization factors for medium 1:\n ", Nmat1)
		print("Norms for medium 1:\n ", Nmed1)
		print("Matrix of normalization factors for medium 2:\n ", Nmat2)
		print("Matrix of mixed normalization factors:\n ", N12)
		print("P matrix:\n ", P)
		print("S matrix:\n ", S)
		print("T matrix:\n ", T)
		print("V matrix:\n ", V)
	else:
		print("Done")
	rc=np.ndarray.flatten(rc)
	tc=np.ndarray.flatten(tc)
	#energy_ref = sum(rc**2)
	#energy_trans = sum(tc**2)
	energy_ref = rc**2
	energy_trans = tc**2
	print("Reflection coefficients (Alsop normalization): ", rc)
	print("Transmission coefficients (Alsop normalization): ", tc)
	print("Fraction of energy transmitted: ", energy_trans)
	rc_proper = rc*np.sqrt(Nmed1[0]/Nmed1)
	tc_proper = tc*np.sqrt(Nmed1[0]/Nmed2)
	# Note that because I've taken Nmed1[0], I'm assuming that incident mode
	# is the fundamental mode.
	print("Transmission coefficients (surface displacement ratio): ", tc_proper)
	return rc_proper,tc_proper, energy_ref, energy_trans

############################################## MAIN PROGRAM #####################################################

modfile1=sys.argv[1]
effile1=sys.argv[2]
modfile2=sys.argv[3]
effile2=sys.argv[4]

# Read incidence-side model
ureo1 = upre.model_1D(modfile1)
dep_pts_mod1 = ureo1.deps_all
dcon1 = ureo1.mod_hif
print("Depths of discontinuity for the incidence side medium: ", dcon1)
# dcon1=input("Depths of discontinuity for the incidence side medium: ")

# Read transmission-side model
ureo2 = upre.model_1D(modfile2)
dep_pts_mod2 = ureo2.deps_all
dcon2 = ureo2.mod_hif
print("Depths of discontinuity for the transmission side medium: ", dcon2)
# dcon2=input("Depths of discontinuity for the transmission side medium: ")
# dcon1=[float(i) for i in dcon1.split()]
# dcon2=[float(i) for i in dcon2.split()]

fstep=-0.005
frange=input('Enter frequency range: ')
fl=float(frange.split()[0])
fh=float(frange.split()[1])
freq=np.arange(fh,fl+fstep,fstep)
allper=1/freq
# for the plotting to work properly allper must be sorted in ascending order
rcoeff=list(range(len(allper)))
tcoeff=list(range(len(allper)))
eng_ref=list(range(len(allper)))
eng_trans=list(range(len(allper)))
for p,per in enumerate(allper):
	rcper,tcper,erefper,etper=do_single_freq(per)
	if p==0:
		# find the number of modes reflected/transmitted at the shortest period
		maxmr=len(rcper)
		maxmt=len(tcper)
	else:
		if len(rcper)<maxmr:
			missing=maxmr-len(rcper)
			rcper=np.append(rcper,np.array(missing*[np.nan]))
			erefper=np.append(erefper,np.array(missing*[np.nan]))
		if len(tcper)<maxmt:
			missing=maxmt-len(tcper)
			tcper=np.append(tcper,np.array(missing*[np.nan]))
			etper=np.append(etper,np.array(missing*[np.nan]))
	rcoeff[p] = rcper
	tcoeff[p] = tcper
	eng_ref[p] = erefper
	eng_trans[p] = etper
#print(eng_ref, eng_trans)
rcm=np.empty((maxmr,len(freq)))
erm=np.empty((maxmr,len(freq)))
tcm=np.empty((maxmt,len(freq)))
etm=np.empty((maxmt,len(freq)))
for i in range(maxmr):
	rcm[i,:]=[rc[i] for rc in rcoeff]
	erm[i,:]=[er[i] for er in eng_ref]
for j in range(maxmt):
	tcm[j,:]=[tc[j] for tc in tcoeff]
	etm[j,:]=[et[j] for et in eng_trans]

fig1=plt.figure()
fig2=plt.figure()
ax1=fig1.add_subplot(111)
ax2=fig2.add_subplot(111)

# plot the reflected/transmitted energy
#for i in range(maxmt):
#	ax1.plot(freq,etm[i],'o-')
#for i in range(maxmr):
#	ax2.plot(freq,erm[i],'o-')
#plt.ylim(0,1.1)
#plt.ylabel('Fraction of incident energy')

# plot the reflection/transmission surface ratios
for i in range(maxmt):
	ax1.plot(freq,tcm[i],'o-')
for i in range(maxmr):
	ax2.plot(freq,rcm[i],'o-')
ax1.set_ylabel('Transmission surface ratio')
ax2.set_ylabel('Reflection surface ratio')
plt.show()
usrc=input("Do you want to save the result ? (y/n) : ")
if usrc=='y':
	jarname="ref.pckl"
	jar=open(jarname,'wb')
	pickle.dump(freq,jar)
	pickle.dump(rcm,jar)
	jar.close()
	jarname2="trans.pckl"
	jar=open(jarname2,'wb')
	pickle.dump(freq,jar)
	pickle.dump(tcm,jar)
	jar.close()
	jarname3="eref.pckl"
	jar=open(jarname3,'wb')
	pickle.dump(freq,jar)
	pickle.dump(erm,jar)
	jar.close()
	jarname4="etrans.pckl"
	jar=open(jarname4,'wb')
	pickle.dump(freq,jar)
	pickle.dump(etm,jar)
	jar.close()
	print("Results stored in %s, %s, %s and %s" %(jarname,jarname2,jarname3,jarname4))
