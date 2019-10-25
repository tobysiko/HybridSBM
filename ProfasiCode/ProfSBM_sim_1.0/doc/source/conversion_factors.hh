//Contents of this file are required only for documentation

/**
 \page convfact Conversion factors between PROFASI units and standard units
1 PROFASI internal temperature unit = 675.32 Kelvin \n
1 PROFASI internal energy unit = 1.342 kcal/mol \n

Using gas constant R = 1.9872159 cal/mol/K \n\n

The conversion factors are calculated by identifying the observed melting temperature of Trp-cage, 1L2Y, in PROFASI and in experiment. For the release version of PROFASI 1.1.2, this melting temperature (obtained by fitting a two state description with the average helix content seen in simulations) is 0.466446 internal PROFASI temperature units. This sets the temperature and energy scales for \e all simulations with PROFASI. There are no free parameters to tune to get the "right" energy or temperature scales any more. 

The scale is subject to change when geometrical properties of some amino acids are changed or modifications are made to the energy function. All but the most experienced users are best advised to use the conversion factors as they are given here, without attempting to modify them in the code. 

*/

