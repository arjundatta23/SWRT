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

# modules written by me
import SW1D_earthsr.utils_pre_code as upre
import SW1D_earthsr.egnfunc_integrals_norms as ein

#######################################################################################

def layers_lateral_discontinuity():

	vs1 = ureo1.oreo.beta
	vs2 = ureo2.oreo.beta
	rho1 = ureo1.oreo.rho
	rho2 = ureo2.oreo.rho

	try:
		assert vs1.size==vs2.size
		assert rho1.size==rho2.size
	except AssertionError:
		raise SystemExit("The two models must have the same 'size'")

	vsdiff=np.abs(vs2-vs1)
	if np.unique(vsdiff)[0]==0:
		# means there is at least one layer (most likely the half-space) for which there is no lateral discontinuity
		uvals=np.unique(vsdiff)[1:]
	else:
		uvals=np.unique(vsdiff)

	# print("uvals: ", uvals)
	beta_it=np.zeros((uvals.size,2))
	rho_it=np.zeros((uvals.size,2))
	veld=[] # veld -> vertical_extent_of_lateral_discontinuity

	for j in range(beta_it.shape[0]):
		ind=np.where(vsdiff==uvals[j])[0]
		# print("IND: ", ind)
		veld.append(ind)
		beta_it[j,0]=vs1[ind[0]]
		beta_it[j,1]=vs2[ind[0]]
		rho_it[j,0]=rho1[ind[0]]
		rho_it[j,1]=rho2[ind[0]]

	# print("veld: ", veld)
	return beta_it, rho_it, veld

#######################################################################################

def do_single_freq(per):

	print("Computing normalization factors for period %f" %(per))
	omega=2*np.pi/per

	# medium 1
	oobj1 = ein.energy_integrals_lov(effile_is, per, dep_pts_is, ishif)
	dep1=oobj1.dep # depth sampling with necessary extra points for medium 1
	mu1=oobj1.mu
	k1=oobj1.kmode
	efmat1=oobj1.l1 # matrix containing eigenfunctions for medium 1
	oobj1.orthogonality_products()
	Nmed1=oobj1.norms
	n=len(Nmed1) # no. of modes in medium 1 (at period per)

	# medium 2
	oobj2 = ein.energy_integrals_lov(effile_ts, per, dep_pts_ts, tshif)
	dep2=oobj2.dep # depth sampling with necessary extra points correspoding to horizontal interfaces of medium 2
	efmat2=oobj2.l1 # matrix containing eigenfunctions for medium 2
	mu2=oobj2.mu
	k2=oobj2.kmode
	oobj2.orthogonality_products()
	Nmed2=oobj2.norms
	m=len(Nmed2) # no. of modes in medium 2

	c=2*np.pi/(k1*per)[0] # means we consider only the fundamental mode incident
	c2=2*np.pi/(k2*per)[0]
	print("Phase velocities: ", c,c2)

	# **************************************************************************
	def get_reftrans_coeffs():
		# for each vertical section of the lateral discontinuity
		rcoeff=np.zeros(numsec)
		tcoeff=np.zeros(numsec)
		for sec in range(numsec):
			root = np.sqrt(1/c**2 + 1/(beta12[sec,1])**2 - 1/(beta12[sec,0])**2)
			numer = mu12[sec,0]/c - (mu12[sec,1]*root)
			denom = mu12[sec,0]/c + (mu12[sec,1]*root)
			rcoeff[sec]=numer/denom
			tcoeff[sec]=1+rcoeff[sec]
		print("ref coeff: ", rcoeff)
		print("trans coeff: ", tcoeff)
		return rcoeff, tcoeff
	# **************************************************************************

	rc,tc=get_reftrans_coeffs()

	# get reflected & transmitted SH displacement on vertical interface
	tc_alldep=np.ones(len(dep1))
	rc_alldep=np.ones(len(dep1))
	# print("length of dep1 is: ", len(dep1))
	foundhif=False
	for sec in range(numsec):
		# start=gisobj.veld[sec][0]
		# end=gisobj.veld[sec][-1]
		start=ldve[sec][0]
		end=ldve[sec][-1]
		# print("start and end original: ", start, end)
		# now 'start' and 'end' are depth indices corresponding to the original depth samples
		# i.e. without the extra depth points added by the ein module
		# they must be modified to account for the horizontal interfaces of the incidence side medium
		# so that the vector tc_alldep is built correctly for multiplication with the incident eigenfunctions
		if foundhif:
			start+=1
			end+=2
			foundhif=False
		if (end+1)*dz in ishif:
			end+=1
			foundhif=True
		# print("start and end modified: ", start, end)
		tc_alldep[start:end+1]=tc[sec]
		rc_alldep[start:end+1]=rc[sec]
	# print(tc_alldep[:100])
	inclov=efmat1[:,0] # means the incident Love wave mode is the fundamental
	shtrans=inclov*tc_alldep
	shref=inclov*rc_alldep
	# plt.plot(inclov,dep1)
	# plt.plot(shtrans,dep1)
	# plt.ylim(500,0)
	# plt.show()

	# Compute the coupling coefficient (integral)
	def integrate(phish,philov,z1,z2,wfn):

		checking=False
		to_add=np.setdiff1d(z2,z1)
		to_remove=np.setdiff1d(z1,z2)
		if len(to_add)>0:
			addid=np.where(z2==to_add[0])[0][0]
			remid=np.where(z1==to_remove[0])[0][0]
			# remove appropriate points
			if checking:
				print("z1 before removal: ", z1[remid-5:remid+5])
				print("phish before removal: ", phish[remid-5:remid+5])
				print("deleted samples ", z1[remid], z1[remid+2])
			z1=np.delete(z1,[remid,remid+2])
			phish=np.delete(phish,[remid,remid+2])
			if checking:
				print("z1 after removal ", z1[remid-5:remid+5])
				print(" phish after removal ", phish[remid-5:remid+5])
			# add appropriate points
			z1=np.insert(z1,[addid,addid+1],[z2[addid],z2[addid+2]]) # careful with usage
			if checking:
				print("phish before additon: ", phish[addid-5:addid+5])
			phish=np.insert(phish,[addid,addid+1],[phish[addid],phish[addid]])   # of numpy insert
			if checking:
				print("added samples ", z2[addid], z2[addid+2])
				print("z1 after addition ", z1[addid-5:addid+5])
				print("phish after addition ", phish[addid-5:addid+5])
		if len(z2) != len(z1) or len(philov) != len(phish):
			sys.exit('Problem with depth sampling of eigenfunctions from the 2 media')

		# now that depth samples are appropriate, do the integration
		sumint=0.0
		for j in range(len(tshif)):
			sid=np.where(dep2==tshif[j])[0][0]
			if j>0:
				sid_prev=np.where(dep2==tshif[j-1])[0][0]
				top=sid_prev+1
			else:
				top=0
			# print("Depth sample excluded from integration: ", dep2[sid])
			# integration above each horizontal interface
			phi_12=philov[top:sid]*phish[top:sid]
			prod=phi_12*wfn[top:sid]
			int_above=spi.simps(prod,dep2[top:sid])
			sumint+=int_above

		# integration below deepest horizontal interface
		phi_12=philov[sid+1:]*phish[sid+1:]
		prod=phi_12*wfn[sid+1:]
		int_below=spi.simps(prod,dep2[sid+1:])

		return sumint + int_below

	cctrans=np.zeros(len(k2))
	ccref=np.zeros(len(k1))
	for tmode in range(len(k2)):
		wt=mu2*(k2[tmode]+k2[tmode])/2
		cctrans[tmode]=(integrate(shtrans,efmat2[:,tmode],dep1,dep2,wt))/omega
		# NB: division of these integrals by omega is not part of the written equations but
		# it is required for consistency with the ein module -- in the Love wave ein module,
		# I have divided the normalising integral by omega so that it equals the energy
		# flux (see notes).
	for tmode in range(len(k1)):
		wt=mu1*(k1[tmode]+k1[tmode])/2
		ccref[tmode]=(integrate(shref,efmat1[:,tmode],dep1,dep1,wt))/omega

	# For energy ref/trans. coefficients we need the amplitude coefficients w.r.t. normalized modes
	tcoef_norm = cctrans/(np.sqrt(Nmed1[0])*np.sqrt(Nmed2))
	etrans=tcoef_norm**2
	# Note we have taken Nmed1[0] which means it is implicit that incident field is ONLY the fundamental mode
	rcoef_norm = ccref/(np.sqrt(Nmed1[0])*np.sqrt(Nmed1))
	eref=rcoef_norm**2
	print("Fraction of energy transmitted: ", sum(etrans))
	# Note that in above line, the sum is sum over modes

	# For surface ratio of displacements, we actually don't need "tcoef_norm" and "rcoef_norm"
	srtrans=cctrans/Nmed2
	srref=ccref/Nmed1
	return etrans, eref, srtrans, ccref

############################################## MAIN PROGRAM #####################################################
# Get user inputs
modfile_is=sys.argv[1]
effile_is=sys.argv[2]
modfile_ts=sys.argv[3]
effile_ts=sys.argv[4]
# is -> incidence side
# ts -> transmission side

#***********************************************************************
# Get all necessary information on model interfaces
#***********************************************************************
print("Identifying model interfaces..")

# Incidence-side model
ureo1 = upre.model_1D(modfile_is)
dep_pts_is = ureo1.deps_all
ishif = ureo1.mod_hif
dz = dep_pts_is[1] - dep_pts_is[0]

# Transmission-side model
ureo2 = upre.model_1D(modfile_ts)
dep_pts_ts = ureo2.deps_all
tshif = ureo2.mod_hif
dz2 = dep_pts_ts[1] - dep_pts_ts[0]

try:
	assert dz==dz2
except AssertionError:
	raise SystemExit("The two models must have the same depth sampling")
else:
	print("Horizontal interfaces on incidence side: ", ishif)
	print("Horizontal interfaces on transmission side: ", tshif)

beta12, rho12, ldve = layers_lateral_discontinuity()
mu12=rho12*(beta12)**2
print("mu on either side of vertical interfaces: ", mu12)
numsec=mu12.shape[0]

#***********************************************************************
fstep=-0.005
frange=input('Enter frequency range: ')
fl=float(frange.split()[0])
fh=float(frange.split()[1])
freq=np.arange(fh,fl+fstep,fstep)
#freq=np.array([0.04])
allper=1/freq
# for the plotting to work properly allper must be sorted in ascending order
# NB: The following VARIABLE NAMES ARE MISLEADING. They sound like reflection and transmission coefficients
# whereas actually they are reflection coefficient and transmission surface ratio.
rcoeff=list(range(len(allper)))
tcoeff=list(range(len(allper)))
eng_trans=list(range(len(allper)))
eng_ref=list(range(len(allper)))

for p,per in enumerate(allper):
	etper,erper,tcper,rcper=do_single_freq(per)
	if p==0:
		maxmr=len(rcper)
		maxmt=len(tcper)
	else:
		if len(rcper)<maxmr:
			missing=maxmr-len(rcper)
			rcper=np.append(rcper,np.array(missing*[np.nan]))
			erper=np.append(erper,np.array(missing*[np.nan]))
		if len(tcper)<maxmt:
			missing=maxmt-len(tcper)
			tcper=np.append(tcper,np.array(missing*[np.nan]))
			etper=np.append(etper,np.array(missing*[np.nan]))
	rcoeff[p] = rcper
	tcoeff[p] = tcper
	eng_trans[p] = etper
	eng_ref[p] = erper

rcm=np.empty((maxmr,len(freq)))
erm=np.empty((maxmr,len(freq)))
tcm=np.empty((maxmt,len(freq)))
etm=np.empty((maxmt,len(freq)))

for i in range(maxmr):
	rcm[i,:]=np.array([rc[i] for rc in rcoeff])
	erm[i,:]=np.array([er[i] for er in eng_ref])
for j in range(maxmt):
	tcm[j,:]=np.array([tc[j] for tc in tcoeff])
	etm[j,:]=np.array([et[j] for et in eng_trans])

#************* All done, now plot and pickle results ***************
fig1=plt.figure()
#first figure is for the transmission surface ratio
fig2=plt.figure()
#second figure is for the energy transmission coefficients
fig3=plt.figure()
#third figure is for the energy reflection coefficients
ax1=fig1.add_subplot(111)
ax2=fig2.add_subplot(111)
ax3=fig3.add_subplot(111)
for i in range(maxmt):
	mnum='Mode %d' %(i)
	ax1.plot(freq,tcm[i],'o-',label=mnum)
	ax2.plot(freq,etm[i,:],label=mnum)

ax2.plot(freq,np.nansum(etm,axis=0),'o-',label='Total')

for i in range(maxmr):
	mnum='Mode %d' %(i)
	ax3.plot(freq,erm[i,:],label=mnum)

ax3.plot(freq,np.nansum(erm,axis=0),'o-',label='Total')
ax1.set_ylabel('Transmission surface ratio')
ax2.set_ylabel('Fraction of energy transmitted')
ax3.set_ylabel('Fraction of energy reflected')
ax1.legend(loc='best')
ax2.legend(loc='best')
ax3.legend(loc='best')
ax2.axhline(y=1,linestyle='--',color='k')
ax2.set_ylim(0,1.1)
plt.show()
usrc=input("Do you want to save the result ? (y/n) : ")
if usrc=='y':
	jarnamet="gregal1974_trans.pckl"
	jar=open(jarnamet,'wb')
	pickle.dump(freq,jar)
	pickle.dump(tcm,jar)
	jar.close()
	jarnamee="gregal1974_etrans.pckl"
	jar2=open(jarnamee,'wb')
	pickle.dump(freq,jar2)
	pickle.dump(etm,jar2)
	jar2.close()
	jarnamer="gregal1974_eref.pckl"
	jar3=open(jarnamer,'wb')
	pickle.dump(freq,jar3)
	pickle.dump(erm,jar3)
	jar3.close()
	print("Results stored in %s, %s and %s" %(jarnamet,jarnamee,jarnamer))
