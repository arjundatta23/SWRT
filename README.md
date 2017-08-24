# sw_reftrans
Code accompanying Datta 2017: "SWRT: A package for semi-analytical solutions of surface wave propagation, including mode conversion, across transversely aligned vertical discontinuities", submitted by Arjun Datta to Computers and Geosciences, 2017

Refer to the above for abbreviations/acronyms used in this README. 

SWRT consists of four programs, implementing the three algorithms discussed in Datta 2017:

1. "method_alsop/implement_alsop.py" - implements the Alsop method (Love waves only)
2. "methods_BodynGreen/implement_bodymethod.py" - implements the GregAl method (Love waves only)
3. "methods_BodynGreen/implement_green_lov.py" - implements the GF method for Love waves
4. "methods_BodynGreen/implement_green_ray.py" - implements the GF method for Rayleigh waves

Modules which are used by all or some of these programs are in the directory "modules_common".

BASIC CODE USAGE

SWRT should run on any system with a standard Python installation. NumPy and SciPy may have to be installed separately, while
matplotlib is required for visualization. All the above programs require the same input (command line) arguments, i.e. two 
files containing the Love or Rayleigh eigenfunctions on either side of the vertical interface. I denote these <ef1> and
<ef2>, corresponding to the incidence side and transmission side media respectively. I call these medium 1 and medium 2.

Beyond the command line arguments, there are two (prompted) user inputs -- first, the horzontal discontinuities in medium 1
and medium 2. These need to be entered (depths in km) as a sequence of numbers, for each of the two media separately. Second, 
the user needs to enter the frequency range (lower and upper bounds of frequency in Hz) in which calculations are to be performed.

For instance, if we consider the example of propagation from left to right in the demo model in Datta 2017 (Figure 1), this
model has horizontal discontunities at depths of 30, 50, 150, 215 km in medium 1, and 7, 50, 150, 215 km in medium 2. To run, say "implement alsop" on this model, do

############
python implement_alsop.py <ef1> <ef2>

Enter depth of horizonatl interfaces for medium 1: 30 50 150 215
Enter depth of horizonatl interfaces for medium 2: 7 50 150 215
Enter frequency range: 0.01 0.15
############

RUNNING THE EXAMPLE "example_modL_rhslyr7".

SWRT contains an example which corresponds to Figures 1 and 2 of Datta 2017. For Love wave propagation in the forward direction, <ef1> and <ef2> are eigen.xdist.0.lov.gz and eigen.xdist.400.lov.gz respectively. The GAl and GF methods can be run on this example (for Love waves) with the following commands:

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

Notice that in the above four commands there is an addtional argument "m0theor.vs.ascii.gz". This is an ascii file containing the Vs structure for the model, in a format used by Roecker et al. (2010). Since I worked with this file format, for my convenience I wrote a module "get_interface_sections.py", used by the "methods_BodynGreen" programs, which goes through the ascii model file and extracts the horizontal interfaces in medium 1 and medium 2, provided the user supplies the lateral location (in km) of the vertical interface. However this preliminary task is completely extraneous to the implemented algorithms and can be done away with with a trivial modification of the code. The horizontal interfaces can be entered manually as with the "implement_alsop.py" program.

With the code in its current form, simply enter the location of the vertical discontinuity, which happens to be 250 km in the example dir.

###########
x-location of vertical discontinuity in model: 250
###########

After enetring the frequency range, this should produce the results needed to reproduce Figure 2 of Datta 2017.
NB: To do the Rayleigh case, simply replace the Love wave eigenfunction files (eigen.xdist.0.lov.gz, eigen.xdist.400.lov.gz) with Rayleigh ones (eigen.xdist.0.ray.gz eigen.xdist.400.ray.gz).

VISUALIATION/PLOTTING SCRIPT

SWRT makes use of the Python pickle module to store the results of any run of a program as a "pickle" which can be loaded later for visualization etc. The script "view_pickles.py" is provided for this purpose. If the result of any of the four programs described above is stored as <pickle name>, figures such as those in Datta 2017 can be made using:

############
python "view_pickles.py" <pickle name 1> <pickle name 2> .... upto any number of stored pickles.
############
