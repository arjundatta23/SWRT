#!/usr/bin/python

import numpy as np

def do_main(P,S,T,V):
	""" sets up and solves the matrix equation Ca=D
	    where 'a' is the vector of reflection coefficients.
            Transmission coefficients 'b' are obtained from 'a'
	"""

	n=P.shape[0]
	m=P.shape[1]
	Pt=np.transpose(P)
	Tt=np.transpose(T)

	mc1=np.dot(P,Tt)
	mc2=np.dot(T,Pt)
#	mc3=np.dot(P,np.dot(V,Pt))
	mc3=np.dot(np.dot(P,V),Pt)

	C = S + mc1 + mc2 + mc3

	#j=int(raw_input("Incident mode number (fundamental is 0: )"))
	j=0
	sd=S[:,j].reshape(n,1)
	md1=mc1[:,j].reshape(n,1)
	md2=mc2[:,j].reshape(n,1)
	md3=mc3[:,j].reshape(n,1)

	D = md2 + md3 - md1 - sd

	a=np.linalg.solve(C,D)

	b = Pt[:,j].reshape(m,1) - np.dot(Pt,a)
	# print("C matrix:\n ", C)
	# print("D matrix:\n ", D)

	return a,b
