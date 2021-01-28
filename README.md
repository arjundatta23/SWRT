# SWRT
SWRT: Surface Wave Refelction and Transmission. This is a suite of codes for semi-analytic solutions of seismic surface wave propgation across sharp lateral discontinuities in 2-D media.  It accompanies Datta (2018):

Datta, A.: SWRT: A package for semi-analytical solutions of surface wave propagation, including mode conversion, across transversely aligned vertical discontinuities, Geosci. Instrum. Method. Data Syst., 7, 101-112, https://doi.org/10.5194/gi-7-101-2018, 2018.

Refer to the above for abbreviations/acronyms and references used in this README.

**********************************************************************************************
A. PACKAGE CONTENTS AND OVERVIEW

SWRT consists of four main programs, implementing the three algorithms discussed in Datta (2018):

1. "method\_alsop/implemen\t_alsop.py" - implements the Alsop method (Love waves only)
2. "methods\_BodynGreen/implement\_bodymethod.py" - implements the GregAl method (Love waves only)
3. "methods\_BodynGreen/implement\_green\_lov.py" - implements the Green's Function (GF) method for Love waves
4. "methods\_BodynGreen/implement\_green\_ray.py" - implements the GF method for Rayleigh waves

All programs rely on the 'SW1D_earthsr' set of modules, which is available as a separate respository:
https://github.com/arjundatta23/SW1D_earthsr

**********************************************************************************************
B. BASIC CODE USAGE

SWRT requires NumPy, SciPy and matplotlib (for visualization). All four programs can be run from the command line.

Apart from the command line arguments, the user is required to input the frequency range (lower and upper bounds of frequency in Hz) in which calculations are to be performed.

=======
Simple command line usage:

python <code_name> <mod\_file\_1> <eigen\_file\_1> <mod\_file\_2> <eigen\_file\_2>
=======
<mod\_file\_1> and <mod\_file\_2> are ASCII text files containing the layered-Earth model on either side of the vertical interface in the 2-D medium, corresponding to the incidence side and transmission side media respectively.
<eigen\_file\_1> and <eigen\_file\_2> are ASCII files containing the corresponding (local) Love or Rayleigh wave eigenfunctions on the two sides.

If you want to use the code as is, all input files are in 'earthsr' format (see example dir) and are therefore read by modules in 'SW1D_earthsr'.

Prompted user input: frequency range (lower and upper bounds in Hz) in which calculations are to be performed.

For example, if we consider the example of propagation from left to right in the demo model in Datta (2018, Figure 1),
to run "implement\_alsop" on this model in the frequency range 0.01 - 0.1 Hz, you would do:

############  
python implement_alsop.py <mod\_file\_1> <eigen\_file\_1> <mod\_file\_2> <eigen\_file\_2>
<SOME CODE OUTPUT>
Enter frequency range: 0.01 0.1  
############  

**********************************************************************************************
C. RUNNING THE EXAMPLE "example\_modL\_rhslyr7".

SWRT contains an example which corresponds to Figures 1 and 2 of Datta (2018). For Love wave propagation in the forward direction, <eigen\_file\_1> and <eigen\_file\_2> are "eigen.xdist.0.lov.gz" and "eigen.xdist.400.lov.gz" respectively. The GAl and GF methods can be run on this example (for Love waves) with the following commands:

For forward direction (medium 1 -> medium 2):

###############  
python ../methods_BodynGreen/implement_bodymethod.py mod.xdist.0.lov eigen.xdist.0.lov.gz mod.xdist.400.lov eigen.xdist.400.lov.gz  
python ../methods_BodynGreen/implement_green_lov.py mod.xdist.0.lov eigen.xdist.0.lov.gz mod.xdist.400.lov eigen.xdist.400.lov.gz  
###############  

For backward direction (medium 1 <- medium 2):

###############  
python ../methods_BodynGreen/implement_bodymethod.py mod.xdist.400.lov eigen.xdist.400.lov.gz mod.xdist.0.lov eigen.xdist.0.lov.gz  
python ../methods_BodynGreen/implement_green_lov.py mod.xdist.400.lov eigen.xdist.400.lov.gz mod.xdist.0.lov eigen.xdist.0.lov.gz  
###############    

After entering the frequency range, this should produce the results needed to reproduce Figure 2 of Datta 2018.  
NB: To do the Rayleigh case, simply replace the Love wave eigenfunction files ("eigen.xdist.0.lov.gz", "eigen.xdist.400.lov.gz") with Rayleigh ones ("eigen.xdist.0.ray.gz", "eigen.xdist.400.ray.gz").

**********************************************************************************************
D. VISUALIATION/PLOTTING SCRIPT

SWRT makes use of the Python pickle module to store the results of any run of a program as a "pickle" which can be loaded later for visualization etc. The script "view_pickles.py" is provided for this purpose. If the result of any of the main programs is stored as {pickle name}, figures such as those in Datta (2018) can be made using:

############  
python "view_pickles.py" {pickle name 1} {pickle name 2} .... upto any number of stored pickles.  
############
