#!/usr/bin/python

# Available python modules
import sys
import pickle
import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt

# modules written by me
import orthonormray as ortr
import get_interface_sections as gis

def do_single_freq(per):
	
	omega=2*np.pi/per

	print "Computing normalization factors for period %f" %(per)

        # medium 1
        oobj1=ortr.ortho(effile1,per,ishif)
        dep1=oobj1.dep # depth sampling with necessary extra points for medium 1
	dep1_orig=oobj1.origdep	# depth sampling without extra points
        mu1=oobj1.mu
        k1=oobj1.kmode
        efmat1_b1=oobj1.orig_b1 # matrix containing eigenfunctions for medium 1
        efmat1_b4=oobj1.orig_b4 # matrix containing eigenfunctions for medium 1
        efmat1_b2=oobj1.b2 
        efmat1_psixx=oobj1.psi_xx 

        Nmed1=oobj1.norms
        n=len(Nmed1) # no. of modes in medium 1 (at period per)

        # medium 2
        oobj2=ortr.ortho(effile2,per,tshif)
        dep2=oobj2.dep # depth sampling with necessary extra points corresponding to horizontal interfaces of medium 2 
	dep2_orig=oobj2.origdep	# depth sampling without extra points
        efmat2_b1=oobj2.orig_b1 # matrix containing eigenfunctions for medium 2
        efmat2_b4=oobj2.orig_b4 # matrix containing eigenfunctions for medium 2
        efmat2_b2=oobj2.b2 
        efmat2_psixx=oobj2.psi_xx 
        mu2=oobj2.mu
        k2=oobj2.kmode
        Nmed2=oobj2.norms
        m=len(Nmed2) # no. of modes in medium 2

	def alter_egn(phi1,phi2,z1,z2,alterwhich):

		checking=False
		if alterwhich==1:
			to_add=np.setdiff1d(z2,z1)
			to_remove=np.setdiff1d(z1,z2)
			if len(to_add)>0:
				addid=np.where(z2==to_add[0])[0][0]
				remid=np.where(z1==to_remove[0])[0][0]
				# remove appropriate points
				if checking:
					print "z1 before removal: ", z1[remid-5:remid+5]
					print "phi1 before removal: ", phi1[remid-5:remid+5]
					print "deleted samples ", z1[remid], z1[remid+2]
				z1=np.delete(z1,[remid,remid+2])
				phi1=np.delete(phi1,[remid,remid+2])
				if checking:
					print "z1 after removal ", z1[remid-5:remid+5]
					print " phi1 after removal ", phi1[remid-5:remid+5]
				# add appropriate points
				z1=np.insert(z1,[addid,addid+1],[z2[addid],z2[addid+2]]) # careful with usage
				if checking:
					print "phi1 before additon: ", phi1[addid-5:addid+5]
				phi1=np.insert(phi1,[addid,addid+1],[phi1[addid],phi1[addid]])   # of numpy insert
				if checking:
					print "added samples ", z2[addid], z2[addid+2]
					print "z1 after addition ", z1[addid-5:addid+5]
					print "phi1 after addition ", phi1[addid-5:addid+5]
				if len(z2) != len(z1) or len(phi2) != len(phi1):
					sys.exit('Problem with depth sampling of eigenfunctions from the 2 media')
				return phi1

		elif alterwhich==2:
			to_add=np.setdiff1d(z1,z2)
			to_remove=np.setdiff1d(z2,z1)
			if len(to_add)>0:
				addid=np.where(z1==to_add[0])[0][0]
				remid=np.where(z2==to_remove[0])[0][0]
				# remove appropriate points
				if checking:
					print "z2 before removal: ", z2[remid-5:remid+5]
					print "phi2 before removal: ", phi2[remid-5:remid+5]
					print "deleted samples ", z2[remid], z2[remid+2]
				z2=np.delete(z2,[remid,remid+2])
				phi2=np.delete(phi2,[remid,remid+2])
				if checking:
					print "z2 after removal ", z2[remid-5:remid+5]
					print " phi2 after removal ", phi2[remid-5:remid+5]
				# add appropriate points
				z2=np.insert(z2,[addid,addid+1],[z1[addid],z1[addid+2]]) # careful with usage
				if checking:
					print "phi2 before additon: ", phi2[addid-5:addid+5]
				phi2=np.insert(phi2,[addid,addid+1],[phi2[addid],phi2[addid]])   # of numpy insert
				if checking:
					print "added samples ", z1[addid], z1[addid+2]
					print "z2 after addition ", z2[addid-5:addid+5]
					print "phi2 after addition ", phi2[addid-5:addid+5]

				if len(z2) != len(z1) or len(phi2) != len(phi1):
					sys.exit('Problem with depth sampling of eigenfunctions from the 2 media')
				return phi2
		
		
	
	def integrate(phi1,phi2,z,discon):
		if discon:
			sumint=0.0	
			for j in range(len(tshif)):
				sid=np.where(dep2==tshif[j])[0][0]
				if j>0:
					sid_prev=np.where(dep2==tshif[j-1])[0][0]
					top=sid_prev+1
				else:
					top=0
				#print "Depth sample excluded from integration: ", dep2[sid]
				# integration above each horizontal interface
				prod=phi2[top:sid]*phi1[top:sid]
				#prod=phi_12*wfn[top:sid]
				int_above=spi.simps(prod,dep2[top:sid])
				sumint+=int_above

			# integration below deepest horizontal interface
			prod=phi2[sid+1:]*phi1[sid+1:]
			#prod=phi_12*wfn[sid+1:]
			int_below=spi.simps(prod,dep2[sid+1:])

		        integral = sumint + int_below
		else:	
			prod=phi1*phi2
			integral=spi.simps(prod,z)

		return integral
	
	b1_med1=range(len(Nmed1)); b2_med1=range(len(Nmed1)); b4_med1=range(len(Nmed1)); psixx_med1=range(len(Nmed1))
	b1_med2=range(len(Nmed2)); b2_med2=range(len(Nmed2)); b4_med2=range(len(Nmed2)); psixx_med2=range(len(Nmed2))
	
	print "Length of Nmed1 and 2: ", len(Nmed1), len(Nmed2)
	hm=max(len(Nmed1),len(Nmed2))
	P=np.zeros((hm,hm)); S0=np.zeros(hm)
	# have chosen the above sizes for P and S, despite the fact that one of the media may have less modes than hm,
	# because I find it easier and more logical to build a 'full size' system of equations
	# (as if both media have the same number of modes) - the reflection/tranmission coefficients
	# for modes that do not exist, will turn out to be 0 anyway, because of the zeros in P and/or S. 
	
	for i in range(len(Nmed1)):
		b1_med1[i]=efmat1_b1[:,i]
		b4_med1[i]=efmat1_b4[:,i]
		b2_med1[i]=efmat1_b2[:,i]
		psixx_med1[i]=efmat1_psixx[:,i]
	for j in range(len(Nmed2)):
		b1_med2[j]=efmat2_b1[:,j]
		b4_med2[j]=efmat2_b4[:,j]
		b2_med2[j]=efmat2_b2[:,j]
		psixx_med2[j]=efmat2_psixx[:,j]

	# NB: b1 and b4 are the original eigenfunctions from ortr, but psixx (because it involves elastic constants)
	# and b2 (because it goes with psixx) have the depth sampling corresponding to the horizontal 
	# discontinuities of the medium to which they belong.

	# Therefore, integration of b1 and b4 can be done directly (with dep1_orig, which is equal to dep2_orig),
	# whereas that of psixx and b2 needs to be done more carefully.

	for m in range(len(Nmed1)):
		for n in range(len(Nmed2)):
			# first, make sure that b2 from medium a has the right depth sampling
			# to multiply with psixx from medium b
			b2_1to2=alter_egn(b2_med1[m],psixx_med2[0],dep1,dep2,1)
			b2_2to1=alter_egn(b2_med2[n],psixx_med1[0],dep2,dep1,1)
	
			I_b2psixx_12=integrate(b2_1to2,psixx_med2[n],dep2,True)
			I_b2psixx_21=integrate(b2_2to1,psixx_med1[m],dep1,True)
			I_b1b4_12=integrate(b1_med1[m],b4_med2[n],dep1_orig,False)
			I_b1b4_21=integrate(b1_med2[n],b4_med1[m],dep1_orig,False)
			# NB: division of these integrals by omega not required in this case
			# because the Rayleigh wave ort module does not do it.

			P[m][n]=abs(I_b2psixx_12 - I_b2psixx_21 + I_b1b4_12 - I_b1b4_21)/(2*Nmed1[m])
			if m==0:
			# remember whilst we compute the entire P[m][n] matrix, we only need the
			# S[0][n] values, for the case of the fundamental mode being incident.
				S0[n]=abs(I_b2psixx_12 + I_b2psixx_21 - I_b1b4_12 - I_b1b4_21)/(2*Nmed2[n])

	#tcoef=S0[0]/(1+(P[0][0])**2)
	#rcoef=-tcoef*(P[0][0])
	#this is the solution in case of coupling of the fundamental mode to itself only

	# now build the linear system of equations which needs to be solved to obtain the complete solution
	# see comment above about first building 'full size' (hypothetical) equation matrices
	lhs = np.zeros((2*hm,2*hm)); rhs = np.zeros(2*hm)
 
	rhs[hm:]=S0
	# the top part (or top half, if same number of modes exist in both media)
	# of the RHS column vector is all zeros.

	lhs[0:P.shape[0],0:P.shape[1]]=P
	# upper left sub-matrix is simply the P-matrix
	lhs[P.shape[0]:,P.shape[1]:]=-np.transpose(P)
	# lower right sub-matrix is the conjugate transpose of the P-matrix
	# (in this code P is real valued so transpose is sufficient)
	lhs[0:P.shape[0],P.shape[1]:]=np.eye(P.shape[1])
	lhs[P.shape[0]:,0:P.shape[1]]=np.eye(P.shape[1])
	# upper right and lower left sub-matrices are identity matrices

	#print "LHS of equation: ", np.matrix(lhs)
	#print "RHS of equation: ", rhs

	soln=np.linalg.solve(lhs,rhs)
	tcoef=soln[:hm]
	rcoef=soln[hm:]
	# top half of the solution vector is transmission coefficients, bottom half is reflection coefficients.

	# a simple check (necessary but not sufficient) to test that everything has worked correctly
	if np.count_nonzero(rcoef) != len(Nmed1) or np.count_nonzero(tcoef) != len(Nmed2):
		# zero-valued elements in rcoef or tcoef will correspond to modes that do not exist
		sys.exit('You have a problem. Quitting.')

	# for modes that do not exist, change the ref/trans coef values from 0 to nan
	if len(Nmed1) != len(Nmed2):
		extra = (hm - len(Nmed2)) if len(Nmed1)>len(Nmed2) else (hm-len(Nmed1))
		if len(Nmed1)>len(Nmed2):
		# number of transmission coefficients should be less than number of reflection coefficients
			tcoef[-extra:]=np.nan
		elif len(Nmed2)>len(Nmed1):
		# number of reflection coefficients should be less than number of transmission coefficients
			rcoef[-extra:]=np.nan
		
	# Unlike my related Love wave codes, in this code the transmission surface ratio IS EQUAL
	# to the transmission coefficient, because the full norm has been used in derving the expressions,
	# rather than assuming unit norm
	srtrans=tcoef #*np.sqrt(Nmed1[0]/Nmed2)
	print "srtrans is: ", srtrans
	return rcoef, srtrans
####################### Main program begins ###########################
# Get user inputs
mod2dfile=sys.argv[1]
effile1=sys.argv[2]
effile2=sys.argv[3]
vdxloc=float(raw_input("x-location of vertical discontinuity in model: "))
ichoice=raw_input("Consider incidence from 'right-hand side' ? (y/n): ")
if ichoice=='y':
	print "WARNING: eigenfunction files (command line arguments) must be entered in the correct order"
	print "i.e. first incidence side file, then transmission side file"
	flip=True
else:
	flip=False

########################################################################
# Get all necessary information on model interfaces
########################################################################
print "Identifying model interfaces.."
gisobj=gis.discontinuities(mod2dfile,vdxloc,flip)
ishif=gisobj.ishif
tshif=gisobj.tshif
print "Horizontal interfaces on incidence side: ", ishif
print "Horizontal interfaces on transmission side: ", tshif
beta12=gisobj.beta_it
rho12=gisobj.rho_it
mu12=rho12*(beta12)**2
print "mu on either side of vertical interfaces: ", mu12
numsec=mu12.shape[0]

########################################################################
fstep=-0.005
frange=raw_input('Enter frequency range: ')
fl=float(frange.split()[0])
fh=float(frange.split()[1])
freq=np.arange(fh,fl+fstep,fstep)
#freq=np.array([0.02])
allper=1/freq
# for the plotting to work properly, 'allper' must be sorted in ascending order
#per=allper[0]
rcoeff=range(len(allper))
tcoeff=range(len(allper))
for p,per in enumerate(allper):
	rcper,tcper=do_single_freq(per)
	if p==0:
	# shortest period
                maxmr=len(rcper)
                maxmt=len(tcper)
        else:
                if len(rcper)<maxmr:
                        missing=maxmr-len(rcper)
                        rcper=np.append(rcper,np.array(missing*[np.nan]))
                if len(tcper)<maxmt:
                        missing=maxmt-len(tcper)
                        tcper=np.append(tcper,np.array(missing*[np.nan]))
        rcoeff[p] = rcper
        tcoeff[p] = tcper

rcm=np.empty((maxmr,len(freq)))
tcm=np.empty((maxmt,len(freq)))
for i in range(maxmr):
        rcm[i,:]=np.array([rc[i] for rc in rcoeff])
for j in range(maxmt):
        tcm[j,:]=np.array([tc[j] for tc in tcoeff])	

for i in range(maxmt):
	mnum='Mode %d' %(i)
        plt.plot(freq,tcm[i],'o-',label=mnum)
plt.show()

usrc=raw_input("Do you want to save the result ? (y/n) : ")
if usrc=='y':
        jarnamet="green_trans.pckl"
        jar=open(jarnamet,'w')
        pickle.dump(freq,jar)
        pickle.dump(tcm,jar)
        jar.close()
        print "Results stored in %s" %(jarnamet)
