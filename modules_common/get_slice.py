#!/usr/bin/python

import os
import sys
import gzip
import numpy as np
import matplotlib.pyplot as plt

###################################################################################################################

def usage():
	print "Usage-"
	print "This script takes either two or three arguments:"
	print "1. Input file in Steve's Gprofile section format - this can be a model parameter file or a spectral"
	print "		domain wavefield file. Both these file types have the same format but differ in no. of columns"
	print "2. x coordinate along path (in km) - distance from model origin at which you want to extract values"
	print "3. z coordinate - optional argument. If this is present, the script assumes that the input file is of"
	print "		wavefield type & extracts the value of the wavefield at the point (x,z)."
	print "		If this argument is absent, the script assumes that the input file is of model type and extracts"
	print "		the model at all depths for given x."
	print "		You only need to supply the file for one model parameter (hard coded to be Vs at the moment)"
	print "		the script automatically looks for matching Vp and Rho files in the same directory & thus"
	print "		extracts Vp, Vs and Rho." 
	print "Output of the script is a file formatted as an input file to 'earthsr'"
	print "Attenuation (both Qp and Qs) values are extracted from the *.grid file produced by mod1d2ben "
	print "This is because I don't modify attenuation at all when I perturb a 1D model to make it 2D"
	sys.exit()		

###########################################################################################################################

def write_earth_format(wavetype,zpoints,thkns,velp,vels,dens):
	outfile=open(out_name,'w')
	if wavetype=='rayleigh':
		indicate_type = 1 # In 'earthsr', 1 is for rayleigh
		binout_name = binout_name_ray+'\n'
	elif wavetype=='love':
		indicate_type = 2 # and anything <> 1 is for love
		binout_name = binout_name_lov+'\n'
	flattening=0
	ref_per=baseper
	top_line = "%d %d %.2f\n" %(zpoints,flattening,ref_per)
	outfile.write(top_line)
	for j in range(thkns.size):
		row_file = "%11.8f %11.8f %10.8f %10.8f %.6E %.6E\n" %(thkns[j],velp[j]/1000.0,vels[j]/1000.0,dens[j]/1000.0,qp[j],qs[j])
		outfile.write(row_file)
	line_bot1 = "%d\n" %(indicate_type) 
	outfile.write(line_bot1)
	line_bot2 = binout_name # name of output binary from earth
	outfile.write(line_bot2)
	line_bot3 = "%.6f %.6f %d %d\n" %(0,0,0,30) # phase velocity range and modes
	outfile.write(line_bot3)
	line_bot4 = "%d %d %f %f  ! frequency table and nsrc\n" %(1,fsamples,deltaf,f0) # number of sources, frequency samples, freq interval & starting frequency
	outfile.write(line_bot4)
	line_bot5 = "%f  !   source depth\n" %(sd) # source depth
	outfile.write(line_bot5)
	line_bot6 = "%f  !   receiver depth\n" %(0) # receiver depth
	outfile.write(line_bot6)
	line_bot7 = "%d" %(0) # not sure what this is
	outfile.write(line_bot7)
	outfile.close()
	del outfile

##############################################################################################################################

def get_attenuation(file_1dmod):
		qpp=[]
		qss=[]
		inpfile=open(file_1dmod)
		inpfile.readline() # first line describes the grid
		for rl in inpfile:
			qpp.append(float(rl.split()[4]))
			qss.append(float(rl.split()[5]))
		inpfile.close()
		del inpfile
		return qpp, qss
####################################################################################################################

class vertical_slice():

	def __init__(self,filein,x1,x2=None):

		# Get the model size and grid spacing
		if filein.endswith('.gz'):
			ascfile = gzip.GzipFile(filein,'r')
		else:
			ascfile = open(filein,'r')
		linelist = ascfile.readlines()
		header_one = linelist[0]
		gs_z = (float(header_one.split()[-1])/1000.0) # input is in m
		gs_x = (float(header_one.split()[-2])/1000.0) # input is in m
		header_two = linelist[1]
		gp_x  = int(header_two.split()[0])
		gp_z = int(header_two.split()[1])
		
		if len(linelist[2].split())==1:
			self.typefile=1
		elif len(linelist[2].split())==6:
			self.typefile=2
		else:	
			print "Unrecognized file type !!"
			sys.exit()
		xnum1=int(x1/gs_x)
		self.xcoord_start = xnum1*gs_x
		try:
			xnum2=int(x2/gs_x)
			xcoord_end = xnum2*gs_x
			self.xpts=np.arange(self.xcoord_start,xcoord_end+gs_x,gs_x)
		except TypeError:
			self.xpts=np.array([self.xcoord_start])
		self.points_x = gp_x
		self.points_z = gp_z
		self.deltaz = gs_z
		del linelist[0:2]
		calc_depth = np.zeros(gp_z)
		# the above attributes or variables exist regardless of file type
		if self.typefile==1:
			file2 = filein.replace('.vs.','.vp.')
			file3 = filein.replace('.vs.','.rho.')
			if not os.path.exists(file2) or not os.path.exists(file3):
				print "Problem finding corresponding vp and rho files"
				print "Cannot proceed with model extraction without the 3 files for vp, vs & rho"
				sys.exit()
			self.tkns = np.zeros(gp_z)
			if file2.endswith('.gz'):
				ascfile_vp = gzip.GzipFile(file2,'r')
				ascfile_rho = gzip.GzipFile(file3,'r')
			else:
				ascfile_vp = open(file2,'r')
				ascfile_rho = open(file3,'r')
			list_lines_vp = ascfile_vp.readlines()
			list_lines_rho = ascfile_rho.readlines()
			del list_lines_vp[0:2]
			del list_lines_rho[0:2]
			self.vs = self.do_read(xnum1,linelist,calc_depth)
			self.depths = calc_depth
			self.vp = self.do_read(xnum1,list_lines_vp,calc_depth)
			self.rho = self.do_read(xnum1,list_lines_rho,calc_depth)
			#print "Depth	Vp	  Vs      Rho"
			for j in range(calc_depth.size):
				if j < (gp_z-1):
					self.tkns[j] = calc_depth[j+1] - calc_depth[j]
				if j == (gp_z-1):
					self.tkns[j] = 0
				#row_screen = "%6.1f %9.6f %8.6f %8.6f" %(calc_depth[j],vp[j],vs[j],rho[j])
				#print row_screen	
			ascfile_vp.close()
			ascfile_rho.close()
			del ascfile_vp, ascfile_rho
		elif self.typefile==2:
			if __name__=='__main__':
				print "Extracting vertical slice between %.1f and %.1f km" %(self.xpts[0],self.xpts[-1])
			xpnint=[int(j/gs_x) for j in self.xpts]
			[ucomp,vcomp,wcomp] = self.do_read(xpnint,linelist,calc_depth)
			self.depths = calc_depth
			self.phaseu = np.angle(ucomp)
			self.phasev = np.angle(vcomp)
			self.phasew = np.angle(wcomp)
			self.ampu = np.abs(ucomp)
			self.ampv = np.abs(vcomp)
			self.ampw = np.abs(wcomp)
			self.realu = ucomp.real
			self.realv = vcomp.real
			self.realw = wcomp.real
			self.imagu = ucomp.imag
			self.imagv = vcomp.imag
			self.imagw = wcomp.imag
			self.uwhole = ucomp
			self.vwhole = vcomp
			self.wwhole = wcomp			
		
		ascfile.close()
		del ascfile

	def do_read(self,xpn,linelist,depth_cal):

		""" Currently this function works differently depending on whether it is reading a model file or a wavefield file.
		    It is designed to extract a vertical slice of the model at a single x-location, or of the wavefield at any number
		    of x-locations """

		z_done = 0
		if self.typefile==1:
			modval = range(self.points_z)
			for i,line in enumerate(linelist):
				lnum = xpn + z_done*self.points_x
				if i==lnum:
					depth_cal[z_done] = z_done*self.deltaz
					modval[z_done] = float(line.split()[0])
					z_done += 1
			return modval
		elif self.typefile==2:
			k=0
			uvalues = np.zeros((self.points_z,len(xpn)), dtype=complex)
			vvalues = np.zeros((self.points_z,len(xpn)), dtype=complex)
			wvalues = np.zeros((self.points_z,len(xpn)), dtype=complex)		
			for i,line in enumerate(linelist):
				lnum1 = xpn[0] + z_done*self.points_x
				lnum2 = xpn[-1] + z_done*self.points_x
				if i in range(lnum1,lnum2+1):
					k+=1
					if k==1:
						depth_cal[z_done] = z_done*self.deltaz
					[u_rp,u_ip] = [float(line.split()[0]), float(line.split()[1])]
					[v_rp,v_ip] = [float(line.split()[2]), float(line.split()[3])]
					[w_rp,w_ip] = [float(line.split()[4]), float(line.split()[5])]
					uvalues[z_done,k-1] = complex(u_rp,u_ip)
					vvalues[z_done,k-1] = complex(v_rp,v_ip)
					wvalues[z_done,k-1] = complex(w_rp,w_ip)
				if i==lnum2:		
					z_done += 1
					k=0
				if z_done > self.points_z:
					sys.exit("Problem extracting vertical slice. Aborting...")
			return uvalues, vvalues, wvalues
	
###############################################################################################################################

class depth_slice():

	def __init__(self,ascfile,z,x1,x2):
		
		# Get the model size and grid spacing
		fc = open(ascfile,'r')
		linelist = fc.readlines()
		header_one = linelist[0]
		gs_z = (float(header_one.split()[-1])/1000.0) # input is in m
		gs_x = (float(header_one.split()[-2])/1000.0) # input is in m
		header_two = linelist[1]
		gp_x  = int(header_two.split()[0])
		gp_z = int(header_two.split()[1])
		if len(linelist[2].split())==1:
			typefile=1
		elif len(linelist[2].split())==6:
			typefile=2
		else:	
			print "Unrecognized file type !!"
			sys.exit()
		del linelist[0:2]
		if x2>x1:
			xnum2=int(x2/gs_x)
		else:
			xnum2=gp_x
		self.xcoord_end = xnum2*gs_x
		xnum1=int(x1/gs_x)
		self.xcoord_start = xnum1*gs_x
		znum = int(z/gs_z)
		self.zcoord = znum*gs_z
		self.xdist=range(xnum2-xnum1)
		if __name__=='__main__':
			print "Extracting values at z = %.1f km, between x = %.1f km and x = %.1f km" %(self.zcoord,self.xcoord_start,self.xcoord_end)
		lnum_start = znum*gp_x + xnum1
		lnum_end = znum*gp_x + xnum2
		points_done=0
		if typefile==2:
			self.phaseu=[];self.phasev=[];self.phasew=[];self.ampu=[];self.ampv=[];self.ampw=[]
			self.realu=[];self.realv=[];self.realw=[]
			self.uwhole=[]; self.vwhole= [] ; self.wwhole=[]
			for i,line in enumerate(linelist):
				if (i>=lnum_start and i<lnum_end):
					self.xdist[points_done] = self.xcoord_start + points_done*gs_x
					[u_rp,u_ip] = [float(line.split()[0]), float(line.split()[1])]
					[v_rp,v_ip] = [float(line.split()[2]), float(line.split()[3])]
					[w_rp,w_ip] = [float(line.split()[4]), float(line.split()[5])]
					#print u_rp, u_ip, v_rp, v_ip, w_rp, w_ip
					self.realu.append(u_rp)
					self.realv.append(v_rp)
					self.realw.append(w_rp)
					self.uwhole.append(complex(u_rp,u_ip))
					self.vwhole.append(complex(v_rp,v_ip))
					self.wwhole.append(complex(w_rp,w_ip))
					self.phaseu.append(np.angle(self.uwhole[-1]))
					self.phasev.append(np.angle(self.vwhole[-1]))
					self.phasew.append(np.angle(self.wwhole[-1]))
					self.ampu.append(np.abs(self.uwhole[-1]))
					self.ampv.append(np.abs(self.vwhole[-1]))
					self.ampw.append(np.abs(self.wwhole[-1]))
					points_done += 1
		elif typefile==1:
			self.vel=[]
			#print "lstart and lend are: ", lnum_start, lnum_end
			for i,line in enumerate(linelist):
				if(i>=lnum_start and i<lnum_end):
					self.xdist[points_done] = self.xcoord_start + points_done*gs_x
					v=float(line)
					vexact=float("%.1f" %(v))
					self.vel.append(vexact)
					points_done += 1

###############################################################################################################################
if __name__ == '__main__':
	if len(sys.argv) <3 or len(sys.argv) > 5:
		usage()
		sys.exit()
	infile = sys.argv[1]
		
	if len(sys.argv) < 5:
		xstart = float(sys.argv[2])
		try:
			xend=float(sys.argv[3])
		except IndexError:
			xend=None
		vertical_extraction = True
		lateral_extraction = False
	else:
		zvalue = float(sys.argv[2])
		xstart = float(sys.argv[3])
		xend = float(sys.argv[4])
		vertical_extraction = False
		lateral_extraction = True	
	
	if vertical_extraction:
		vsobj = vertical_slice(infile,xstart,xend)
		if vsobj.typefile==1:
			baseper = float(raw_input("Reference period for earthsr: "))
			sd = float(raw_input("Enter depth of source: "))
			deltaf = float(raw_input("Enter frequency increment: "))
			fsamples = int(raw_input("Enter no. of freq. samples: "))
			f0 = float(raw_input("Enter starting frequency: "))
			binout_name_ray = raw_input("Name of binary output of excitations (Rayleigh): ")
			binout_name_lov = raw_input("Name of binary output of excitations (Love): ")
			[qp, qs] = get_attenuation('homomodel.txt.grid')
			out_name = 'mod.xdist.%d.ray' %(int(vsobj.xcoord_start))
			write_earth_format('rayleigh',vsobj.points_z,vsobj.tkns,vsobj.vp,vsobj.vs,vsobj.rho)
			print "For Rayleigh, model written to file: ", out_name
			out_name = 'mod.xdist.%d.lov' %(int(vsobj.xcoord_start))
			write_earth_format('love',vsobj.points_z,vsobj.tkns,vsobj.vp,vsobj.vs,vsobj.rho)
			print "For Love, model written to file: ", out_name
		elif vsobj.typefile==2:
			if xend==None:
				plt.plot(vsobj.realw,vsobj.depths,label='Displacement')
				plt.plot(vsobj.ampw,vsobj.depths,label='Amplitude')
				plt.ylim([vsobj.depths[-1],0])
				plt.ylabel('Depth [Km]')
				plt.grid(True)
				plt.title('Location %.1f km' %(xstart))
				plt.legend(loc=4)
				plt.show()
	elif lateral_extraction:
		dsobj = depth_slice(infile,zvalue,xstart,xend)
		print "Distances are: ", dsobj.xdist
		print "Actual displacement (vertical component): ", dsobj.realw
		print "Phase values (vertical component) are: ", dsobj.phasew
		print "Amplitudes (vertical component) are: ", dsobj.ampw
		plt.plot(dsobj.xdist,dsobj.phasew,'-o')
		plt.show()
	
