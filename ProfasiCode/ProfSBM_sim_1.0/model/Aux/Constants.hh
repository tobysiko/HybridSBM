/*******************************************************************************
    PROFASI: Protein Folding and Aggregation Simulator, Version 1.5
    Copyright (C) (2012)  Anders Irback and Sandipan Mohanty
    Email: profasi@thep.lu.se
    Home Page: http://cbbp.thep.lu.se/activities/profasi/
    Version control (git) : https://trac.version.fz-juelich.de/PROFASI

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License
    (see PROFASI/gpl.txt).

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
********************************************************************************/

#ifndef Constants_HH
#define Constants_HH

//! Set of frequently used constants

namespace UnivConstants
{
    //! pi
    const double pi=3.1415926535897932385;
    //! 2*pi
    const double pi2=6.2831853071795864769;
    const double twoPi=6.2831853071795864769;
    //! 0.5*pi
    const double pid2=1.5707963267948966192;
    //! 2*pi/3
    const double twoPid3=2.0943951023931954923;
    //! 1 radian in degrees
    const double radian_in_degrees=180/pi;

    //! Temperature conversion factor: 1 profasi units in Kelvin.
    const double pru_in_kelvin=2000.0/3;
    //! Temperature conversion factor: 1 Kelvin in profasi units
    const double kelvin_in_pru=1.0/pru_in_kelvin;
    //! Energy conversion factor: 1 profasi unit in kcal/mol
    const double prf_energy_in_kcal_per_mol=1.9858*pru_in_kelvin/1000;
    //! Energy conversion factor: 1 kcal/mol in profasi unit
    const double kcal_per_mol_in_prf_energy= 1.0/prf_energy_in_kcal_per_mol ;  // added by TS

    //! Minimum Ramachandran phi angle for helix
    const double helix_phimin=-90*pi/180.0  ;
    //! Maximum Ramachandran phi angle for helix
    const double helix_phimax=-30*pi/180.0  ;
    //! Minimum Ramachandran psi angle for helix
    const double helix_psimin=-77*pi/180.0  ;
    //! Maximum Ramachandran psi angel for helix
    const double helix_psimax=-17*pi/180.0  ;

    //! Minimum Ramachandran phi angle for beta strand
    const double sheet_phimin=-150*pi/180.0 ;
    //! Maximum Ramachandran phi angle for beta strand
    const double sheet_phimax=-90*pi/180.0  ;
    //! Minimum Ramachandran psi angle for beta strand
    const double sheet_psimin=90*pi/180.0   ;
    //! Maximum Ramachandran psi angel for beta strand
    const double sheet_psimax=150*pi/180.0  ;
    
    //const double kBoltz = 1.3806488e-23 ;// in J/K
    const double kBoltz = 0.0019872041; // kcal/mol/K     // added by TS
    //const double kBoltzProfasi = kBoltz / kcal_per_mol_in_prf_energy * kelvin_in_pru ; 
}

#endif
