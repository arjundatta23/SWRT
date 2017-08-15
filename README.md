# sw_reftrans
Code accompanying Datta 2017: "SWRT: A package for semi-analytical solutions of surface wave propagation, including mode conversion,
across transversely aligned vertical discontinuities", submitted by Arjun Datta to Computers and Geosciences, 2017

SWRT consists of four programs, implementing the three algorithms discussed in Datta 2017:

1. "implement_alsop.py" - implements the Alsop method (Love waves only)
2. "implement_body_method.py" - implements the GregAl method (Love waves only)
3. "implement_green_lov.py" - implements the GF method for Love waves
4. "implement_green_ray.py" - implemrnts the GF method for Rayleigh waves

Modules which are used by all or some of these programs are in the directory "common_modules".

How to use the code:

SWRT should run on any system with a standard Python installation. NumPy and SciPy may have to be installed separately, while
matplotlib is required for visualization. All the above programs require the same input (command line) arguments, i.e. two 
files containing the Love or Rayleigh eigenfunctions on either side of the vertical interface. I denote these <ef1> and
<ef2>, corresponding to the incidence side and transmission side media respectively. I call these medium 1 and medium 2.

Beyond the command line arguments, there is are two (prompted) user inputs -- first, the horzontal discontinuities in medium 1
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
