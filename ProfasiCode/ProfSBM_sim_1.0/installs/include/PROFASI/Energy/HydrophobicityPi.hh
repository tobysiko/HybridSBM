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

#ifndef HydrophobicityPi_HH
#define HydrophobicityPi_HH
#include "Energy.hh"
#include <map>


namespace prf
{
//! Effective hydrophobic attraction between non-polar side chains
    /**
     * This class represents the hydrophobicity interaction in the model. In its
     * current version, hydrophobicity is a simple pair-wise additive potential.
     * A degree of contact between two hydrophobic side chains is calculated
     * based on the proximity of a pre-determined set of atoms from each
     * side chain. This degree of contact is multiplied by a strength assigned
     * to pairs of side chain, which are in rough correspondence with their size.
     * Like other terms, the implementation of the class takes care of the
     * optimization. Contribution from each pair of hydrophobic side-chains is
     * saved in a contribution matrix. For each update, the affected
     * contributions are determined, and reevaluated, and compared with the
     * value stored in the matrix. Contributions from hydrophobic pairs that
     * are nearest and second-nearest neighbours in sequence are supressed by
     * factors 0.0 and 0.5 respectively.
     * \ingroup profasi_energies
     */

    class HydrophobicityPi:public Energy
    {
    public:
        HydrophobicityPi();
        ~HydrophobicityPi();
        void init();
        double evaluate();
        double gradientXYZ(std::valarray<double> &gx);
        double deltaE(Update *);
        void Accept(Update *);
        void DisplayMatrix();
        void rangeEstimate(double &x1, double &x2);
        inline void SetNNStrength(double n1, double n2) {
            strnn = n1;
            strnnn = n2;
        }

        inline void ScaleStrength(double x) {
            sclfct = fabs(x);
        }

        inline double contribution(int i, int j) const {
            return sclfct*Mehp[NHAA * i + j];
        }

        inline int nhaa() const {
            return NHAA;
        }
        inline void setfocalRes(int res){focalRes=res;}
        inline void setfocalRes2(int res){focalRes2=res;}
        inline void setfocalResScale(double factor){focalResScale=factor;}
        
        double interhp();
        //! i and j refer to the list of hydrophobic residues in the whole system.
        double contact_frac(int i, int j,double(*distf)(int,int));
        //! i and j refer to the list of hydrophobic residues in the whole system.
        double arpi(int i, int j, bool &success);
        double lookupPi(std::string key, double r, double theta, double phi, bool &success);
        //! iaa and jaa refer to the ligand indices relative to the whole system
        double hp_contact_frac(int iaa, int jaa,double(*distf)(int,int));
        double InterChain(int ich, int jch);
        
        void initPiInteractionMap();
        inline void setCatPiStrength(double val){arpistrength = val;};
        inline void setCatPiRepulsion(double val){cprep = val;};
        inline void setMode(std::string m){mode = m;}
        inline std::string getMode(){return mode;}
    private:
    	int focalRes,focalRes2;
    	double focalResScale,arpistrength, cprep;
        int htyp(OneLetterCode cod);
        int ptyp(OneLetterCode cod);
        int ctyp(OneLetterCode cod);
        int groupOf(OneLetterCode cod);
        double hp_pair(int i, int j,double(*distf)(int,int));
        double dhp_pair(int i, int j, bool perhint,
                        std::valarray<double> &gx);
        double getIndex(double num, double min, double factor);
        Vector3 centroid(Ligand* lig);
        Vector3 normal(Ligand* lig);
        //double centroidDistance(int i, int j);
        //double theta(int i, int j);
        //double phi(int i, int j);
        std::valarray < double > Mehp, Vehp;
        std::valarray<int> changed;
        int NHAA, N, nchanges;
        double a2, b2, d2,slpe;
        double strnn, strnnn, sclfct,max_at_dist2,max_cmp_dst2;
        std::valarray<double> strength;
        std::valarray<int> atom;
        std::valarray<int> natom, ih, chainof, soften, hlig;
        double r2min[12];
        int r2minpartner[12];
        std::valarray <OneLetterCode> seq,holc;
        std::valarray <bool> ispi, iscat;
        std::valarray <Vector3> centroids, normals;
        static const int nHPgrp;
        static const double hpstr[];
        static std::string hatms[][6];
        
        static const int nPPgrp, nCPgrp;
        static const double ppstr[],cpstr[];
        static std::string patms[][8];
        static std::string natms[][7];
        static std::string cenatms[][3];
        static std::string catatms[][2];
        static double pcharge[][8];
        static double ncharge[][7];
        
        double rstep,thetastep,phistep;
        double rmin,rmax,thetamin,thetamax,phimin,phimax;
        std::string key4Pair(Ligand* lig1, Ligand* lig2);
        std::map< std::string,std::map< std::string, double >* > AA2PiMap;
        static const std::string pipimap_filenames[12];
        static const int nPiFilenames;
        std::vector<std::string> pipikeys;
        bool debug;
        Vector3 displacePoint(Vector3 point, Vector3 vect, double distance);
        std::map<std::string, bool> validCatPi;
        //std::Vector<std::string> bannedKeys;
        std::string mode;
    };
}

#endif
