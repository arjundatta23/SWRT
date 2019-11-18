# SWRT
SWRT: Surface Wave Refelction and Transmission. This is a package of code for semi-analytic solutions of seismic surface wave propgation across sharp lateral discontinuities in 2-D media.  It accompanies Datta (2018): 

Datta, A.: SWRT: A package for semi-analytical solutions of surface wave propagation, including mode conversion, across transversely aligned vertical discontinuities, Geosci. Instrum. Method. Data Syst., 7, 101-112, https://doi.org/10.5194/gi-7-101-2018, 2018.

Refer to the above for abbreviations/acronyms used in this README.

**********************************************************************************************
A. PACKAGE CONTENTS AND OVERVIEW

SWRT consists of four main programs, implementing the three algorithms discussed in Datta (2018):

1. "method\_alsop/implemen\t_alsop.py" - implements the Alsop method (Love waves only)
2. "methods\_BodynGreen/implement\_bodymethod.py" - implements the GregAl method (Love waves only)
3. "methods\_BodynGreen/implement\_green\_lov.py" - implements the Green's Function (GF) method for Love waves
4. "methods\_BodynGreen/implement\_green\_ray.py" - implements the GF method for Rayleigh waves

Modules which are used by all or some of these programs are in the directory "modules_common".

**********************************************************************************************
B. BASIC CODE USAGE

SWRT should run on any system with a standard Python installation. NumPy and SciPy may have to be installed separately, while matplotlib is required for visualization. All the four programs can be run from the command line.

Simple command line usage:

python {program name} {ef1} {ef2}

{ef1} and {ef2} are the two files containing Love or Rayleigh wave eigenfunctions on either side of the vertical interface in the model, corresponding to the incidence side and transmission side media respectively. I call these medium 1 and medium 2.

Prompted user inputs:

1. Horzontal discontinuities in medium 1 and medium 2. These need to be entered as a sequence of numbers (depths in km), for each of the two media separately. 
2. Frequency range (lower and upper bounds in Hz) in which calculations are to be performed.

For example, if we consider the example of propagation from left to right in the demo model in Datta (2018, Figure 1), this
model has horizontal discontunities at depths of 30, 50, 150, 215 km in medium 1, and 7, 50, 150, 215 km in medium 2. 
To run "implement\_alsop" on this model in the frequency range 0.01 - 0.1 Hz, you would do:

############  
python implement_alsop.py {ef1} {ef2}

Enter depth of horizonatl interfaces for medium 1: 30 50 150 215  
Enter depth of horizonatl interfaces for medium 2: 7 50 150 215  
Enter frequency range: 0.01 0.1  
############  

**********************************************************************************************
C. RUNNING THE EXAMPLE "example\_modL\_rhslyr7".

SWRT contains an example which corresponds to Figures 1 and 2 of Datta (2018). For Love wave propagation in the forward direction, {ef1} and {ef2} are "eigen.xdist.0.lov.gz" and "eigen.xdist.400.lov.gz" respectively. The GAl and GF methods can be run on this example (for Love waves) with the following commands:

For forward direction:

###############  
python methods_BodynGreen/implement_bodymethod.py m0theor.vs.ascii.gz eigen.xdist.0.lov.gz eigen.xdist.400.lov.gz  
python methods_BodynGreen/implement_green_lov.py m0theor.vs.ascii.gz eigen.xdist.0.lov.gz eigen.xdist.400.lov.gz  
###############  

For backward direction:

###############  
python methods_BodynGreen/implement_bodymethod.py m0theor.vs.ascii.gz eigen.xdist.400.lov.gz eigen.xdist.0.lov.gz  
python methods_BodynGreen/implement_green_lov.py m0theor.vs.ascii.gz eigen.xdist.400.lov.gz eigen.xdist.0.lov.gz  
###############  

NB: For backward propagation, you will also need to answer yes '(y)' to the prompted question
"Consider incidence from right hand side (y/n)?: "

NB: in the above commands there is an addtional argument "m0theor.vs.ascii.gz". This is an ascii file containing the Vs structure for the model, in a format used by Roecker et al. (2010). Since I worked with this file format, for my convenience I wrote a module "get\_interface\_sections.py", used by the "methods\_BodynGreen" programs, which goes through the ascii model file and extracts the horizontal interfaces in medium 1 and medium 2, provided the user supplies the lateral location (in km) of the vertical interface. However this preliminary task is completely extraneous to the implemented algorithms and can be done away with with a trivial modification of the code. The horizontal interfaces can be entered manually as with the "implement\_alsop.py" program.

With the code in its current form, simply enter the location of the vertical discontinuity, which happens to be 250 km in the example dir.

###########
x-location of vertical discontinuity in model: 250
###########

After enetring the frequency range, this should produce the results needed to reproduce Figure 2 of Datta 2018.  
NB: To do the Rayleigh case, simply replace the Love wave eigenfunction files ("eigen.xdist.0.lov.gz", "eigen.xdist.400.lov.gz") with Rayleigh ones ("eigen.xdist.0.ray.gz", "eigen.xdist.400.ray.gz").

**********************************************************************************************
D. VISUALIATION/PLOTTING SCRIPT

SWRT makes use of the Python pickle module to store the results of any run of a program as a "pickle" which can be loaded later for visualization etc. The script "view_pickles.py" is provided for this purpose. If the result of any of the main programs is stored as {pickle name}, figures such as those in Datta (2018) can be made using:

############
python "view_pickles.py" {pickle name 1} {pickle name 2} .... upto any number of stored pickles.
############
