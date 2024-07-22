OS Name: Ubuntu 20.04.6 LTS
OS Type: 64-bit
GNOME Version: 3.36.8

Programming language: Python 3.8.10
All codes were run using IPython 7.3.10
This software comes pre-installed in most of the Linux distributions.

##########################################
All the data is obtained by solving the Allen-Cahn equation described in the supplementary material with the parameters specified in the text. We use periodic boundary conditions and a fast Fourier transform algorithm to discretize the differential equation. The time evolution is performed using the Adams-Bashfort method for the nonlinear terms and Crank-Nicolson for the differential operators, which are treated implicitly (see Eq. 2 in Ascher, Uri M., Steven J. Ruuth, and Brian TR Wetton. "Implicit-explicit methods for time-dependent partial differential equations." SIAM Journal on Numerical Analysis 32.3 (1995): 797-823.). The numerical solver can be provided upon reasonable request. 
##########################################

#############################################
The datasets and codes to analyze the results from the numerical solver are described below.


####################################################################################
Fig4 folder contents:
ExtensionAndTotalMass
####################################################################################

ExtensionAndTotalMass contents:

DropExtensionAnalysis.py takes "TrialData_Speed_StrongerNuc_ForPaper_"+str(j_aux)+".npy" where j_aux is a number defining the data set containing the time evolution of the concentration field in 2D.
With this script we follow the edge of the condensate and measure its extension with time. The output of running

python3 DropExtensionAnalysis.py 

is the values for the extension given by "Extension_"+str(j_aux)+".txt" and time points in which this is evaluated "ExtensionTime.txt". 

LuminiscenceAnalysis.py takes "TrialData_Speed_StrongerNuc_ForPaper_"+str(j_aux)+".npy" where j_aux is a number defining the data set containing the time evolution of the concentration field in 2D.
With this script we sum all the concentration variables within a condensate and how it evolves with time. The output of running

python3 LuminiscenceAnalysis.py 

provides the values of the total condensate material amount and they are given in "NormalizedIntensity"+str(j_aux)+".txt" and time points in which this is evaluated "TimeLuminiscence.txt". 
####################################################################################


####################################################################################
FigS5 folder contents:
MultiDropAnalysis
NumericalCalculationOfVelocity
####################################################################################

####################################################################################
MultiDropAnalysis contents:

DropExtensionAnalysis.py takes as an input "MultipleDrops.npy", which is the data set containing the time evolution of multiply nucleated droplets.
With this script we follow the edges of multiply nucleated condensates and follow their extension with time. This is done until all condensates have fused and extended until reaching the calculation box-size. 

The output of running
python3 DropExtensionAnalysis.py 
is the time arrays tracking the extension of a fixed number of drops until two of them fuse and the extension of each drop for this time. The time arrays are given by "TimeForNumOfDrops"+str(j)+".txt", where j keeps track of the changes in number of drops. Other output files contain the extension values and are of the form "AllExtensionsPerDrop_"+str(j)+"_"+str(k)+".txt" where str(k) identifies each of the drops in a certain time window. 

LuminiscenceAnalysis.py takes as an input "MultipleDrops.npy".
With this script we follow the total amount of material of multiply nucleated condensates and follow it over time. This is done until all condensates have fused and extended until reaching the calculation box-size. 

The output of running
python3 LuminiscenceAnalysis.py 
is the time arrays tracking the extension of a fixed number of drops until two of them fuse and the extension of each drop for this time. The time arrays are given by "TimeForNumOfDropsLum_"+str(j)+".txt", where j keeps track of the changes in number of drops. Other output files contain the total amount of material in each drop and are of the form "AllLuminiscencesPerDrop_"+str(j)+"_"+str(k)+".txt" where str(k) identifies each of the drops in a certain time window. 

####################################################################################

####################################################################################
NumericalCalculationOfVelocity contents:
DropExtensionAnalysis.py takes "TrialData_Speed_StrongerNuc_ForPaper_"+str(j_aux)+".npy" where j_aux is a number defining the data set containing the time evolution of the concentration field in 2D.
With this script we follow the edge of the condensate and measure its extension with time, and then use this information to make a linear fit to extract the velocity of growth condensate. 

The output of running 
python3 DropExtensionAnalysis.py
are "Velocity.txt" that contains the velocity for values of the relative binding affinity that are stored in "BindingAffinities.txt". The errors for each fit are saved as "FittingErrors.txt". 

####################################################################################
