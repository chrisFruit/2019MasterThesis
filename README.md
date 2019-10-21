# 2019MasterThesis
Source code for OpenFOAM modifications

>>> CompositionModel.C:
- The Cp() function which is responsible for calculating particle phase cp values is found between lines #435-483

>>> ReactingMultiphaseParcel.C:
- This file contains both the original OpenFOAM implementation in addition to the Kirov model.
- The use of which model to use in the CpEff() function is dictated by the 'define' statement found on line #59.
- The OpenFOAM cp model is defined between lines #62-67.
- The Kirov Cp model is defined between lines #70-109.

>>> coalFoam/:
- This directory contains the user copy of the 'coalChemistryFoam' solver

>>> DualCompetingRateDevolatilisation/:
- Contains the source and header file for the dual-competing rate solver. 
- The actual model implementation is found in the source file on line #139.

>>> MATLAB_Optimization/:
- Contains the full optimization routine for the kinetic parameters
- The mass loss sub-routine is described between lines #70-133
