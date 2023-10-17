What is InbredMetapop?

InbredMetapops is a spatially and genetically explicit individual-based model to simulate eco-evolutionary dynamics in metapopulations, explicitly accounting for accumulation of genetic load and dispersal evolution.
The model has been created to test the extent to which strong population structure, the build-up of deleterious mutations (i.e. genetic load) and the evolution of dispersal can act synergistically to facilitate reduction of the genetic load in metapopulations and postpone extinctions through mutational meltdowns. For more details, please see https://doi.org/10.1111/evo.14620. 

Parameters

All parameter of the model are defined in the parameter.h file. Parameter values can be changed in parameters.cpp before producing an executable. 
The default parameters are those used in the main text https://doi.org/10.1111/evo.14620.


How to compile InbredMetapop

Simulations were conducted on a high performance computer (Linux). The executable was compiled in Release mode with the compiler gcc (GCC) 4.8.5 20150623 (Red Hat 4.8.5-36).
To compile the code on a Linux system, the macro LINUX in PolyDisp.h must be set to LINUX 1 (the default). To compile on a Windows system, LINUX in PolyDisp.h must be set to LINUX 0.

One way of creating an executable on a Linux system is:
Collect all cpp and h files from this repository into a folder named src. 
Add the Makefile to the directory containing src.
While in this directory, run the command "make all CFG=Release". 
This will create a folder called bin.Release which will contain your executable. 
