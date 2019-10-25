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

#ifndef ExVol_HH
#define ExVol_HH
#include "Energy.hh"
#include "ExVolBase.hh"
#include "../Aux/ConnectionsMatrix.hh"

namespace prf
{
//! The excluded volume term
    /**
     * \ingroup profasi_energies
     * This class represents the excluded volume interactions and is the most
     * painstakingly optimized part of PROFASI. A backup matrix, as used for
     * other energy terms would not be very effective, because almost every
     * pair of atoms can interact with this term. The optimization is achieved
     * using the following ...
     *
     * <ol>
     * <li> Since the model assumes fixed bond lengths and bond angles, there are
     * many pairs of atoms with fixed separation. Excluded volume contribution
     * from such pairs is dynamically uninteresting and is set to zero. </li>
     * <li> Atom pairs separated by 3 covalent bonds are treated separately in the
     * LocExVol class </li>
     * <li> The excluded volume potential in the model is a strong repulsion with
     * a finite range. So, the periodic box for the whole system is divided into
     * cells with linear dimensions greater than or equal to the range of the
     * excluded volume potential. Each atom is assigned a cell. Atom pairs which
     * are neither in the same cell nor in neighbouring cells would have distances
     * greater than the range of the potential, and hence wont contribute. So,
     * during evaluation, such pairs are ignored. For each atom, contribution from
     * other atoms in the same cell and in the neighbouring cells is calculated.
     * For each update, the cell-assignment is carefully updated, so that only the
     * moved atoms are re-assigned. Energy contributions are recalculated only on
     * "active" cells, in which something happened: atoms moved in or moved out of
     * the cell or one of its neighbours. </li>
     * </ol>
     */

    class ExVol:public ExVolBase, public Energy
    {
    public:
        ExVol();
        ~ExVol();
        void init();
        void rangeEstimate(double &x1, double &x2);
        double evaluate();
        double gradientXYZ(std::valarray<double> &ans);
        double deltaE(Update *);
        double deltaEwithlimit(Update * updt, double emax);
        void Accept(Update *);
        void Revert(Update *);
        double sacfull();
        inline double sqr(double xx) {
            return xx * xx;
        }

        inline double boxcrdx(int i) {
            return (AtomCoordinates::xval(i) - neginf) -
                   bxl *((int)((AtomCoordinates::xval(i) - neginf) / bxl));
        }

        inline double boxcrdy(int i) {
            return (AtomCoordinates::yval(i) - neginf) -
                   bxl *((int)((AtomCoordinates::yval(i) - neginf) / bxl));
        }

        inline double boxcrdz(int i) {
            return (AtomCoordinates::zval(i) - neginf) -
                   bxl *((int)((AtomCoordinates::zval(i) - neginf) / bxl));
        }

        inline double dist2a(int i, int j) {
            i*=3;j*=3;
            return sqr(lxyz[i] - lxyz[j]) + sqr(lxyz[i+1] - lxyz[j+1]) +
                   sqr(lxyz[i+2] - lxyz[j+2]);
        }

        inline double dist2b(int i, int j) {
            i*=3;j*=3;
            double xtmp=fabs(lxyz[i]-lxyz[j]);
            double ytmp=fabs(lxyz[i+1]-lxyz[j+1]);
            double ztmp=fabs(lxyz[i+2] - lxyz[j+2]);
            return sqr((xtmp<hbxl) ? xtmp : (bxl - xtmp)) +
                   sqr((ytmp<hbxl) ? ytmp : (bxl - ytmp)) +
                   sqr((ztmp<hbxl) ? ztmp : (bxl - ztmp));
        }

        inline Vector3 vecitoj_a(int i, int j)
        {
            i*=3;j*=3;
            return Vector3(lxyz[j]-lxyz[i],lxyz[j+1]-lxyz[i+1],
                           lxyz[j+2]-lxyz[i+2]);
        }

        inline Vector3 vecitoj_b(int i, int j)
        {
            i*=3;j*=3;
            double xtmp=lxyz[j]-lxyz[i];
            double ytmp=lxyz[j+1]-lxyz[i+1];
            double ztmp=lxyz[j+2] - lxyz[i+2];
            if (xtmp<-hbxl) xtmp+=bxl;
            if (xtmp>hbxl) xtmp-=bxl;
            if (ytmp<-hbxl) ytmp+=bxl;
            if (ytmp>hbxl) ytmp-=bxl;
            if (ztmp<-hbxl) ztmp+=bxl;
            if (ztmp>hbxl) ztmp-=bxl;
            return Vector3(xtmp,ytmp,ztmp);
        }
        double calc_esa_with_grd(std::valarray<double> &gx);
        double cellcalc_with_grd(int inc,int nl,int &incnr,
                                    std::valarray<double> &gx);
        inline double dist_cutoff() const {return cutg;}

        int initcells();

    private:
        double LAMBDA, cutg, cutg2;
        std::valarray < double > sig2, asa, bsa;
        ConnectionsMatrix connected;
        int MAXNGB, NTO, nx, ny, nz, nc, xpitch, ypitch, zpitch;
        int useperiod;
        int imprequest; //impossible deltaE requested
        double cellx, celly, cellz, cellxc, cellyc, cellzc;
        double xmin, ymin, zmin, xtmp;
        double xmax, ymax, zmax, bxl, hbxl, neginf;
        std::valarray < int > listt, pnt; /* pointers to atoms        */
        std::valarray < int > cell;         /* division into cells      */
        std::valarray < double > lxyz, blxyz;
        std::valarray < int > atomloc; /* cell in which atom i is located   */
        std::valarray < int > atlcbk;         /* backup for atom locations */
        std::valarray < int > va_ngbdisp; /* valarray for neighbour cells */
        std::valarray < char > atom_stat_bkp, atom_stat, atom_used;
        //! Different kinds of neighbours for an atom
        /**
         * Under an update an atom could be unchanged, be a part of a subsystem
         * which moves rigidly, or end up in a sybsystem which changes shape.
         * Excluded volume contribution need not be calculated inside the
         * fixed part, or internally for the part that moves rigidly. The
         * different parts are tagged after an update, and such terms which
         * are guaranteed to be fixed are left out in the calculations of
         * deltaE.
         */
        std::valarray < int > unmoved_ngb, rigid_ngb, flxbl_ngb;
        int inc_unmoved, inc_rigid, inc_flxbl, n_unmoved, n_flxbl,
        n_rigid;
        int cellstatus;
        double enew;
    public:
//! Categorizes different kinds of pairs of atoms
        /**
         * If NA is the number of types of atoms (5 in our case, H, C, N, O, S),
         * there are NA*NA kinds of pairs. This information is used to pre-calculate
         * the unchanging parts of the functional form of excluded volume, and only
         * insert the distance information for a new configuration when required. So,
         * the function PairType(i,j) here, returns numbers between -1 and 25. 0-24
         * for all kinds of pairs. -1 means "don't care", and is returned when there
         * is less than or equal to three covalent bonds between them. With 1 or 2
         * bonds, the distance between the pair can not change in the model. For 3
         * covalent bonds, the excluded volume is calculated separately in the
         * LocExVol class.
         */
        inline int PairType(int i, int j) {
            return connected(i,j) ? -1 :p->PairType(i,j);
        }

    private:
        double cellcalc(int inc, int nl, int &incnr);
        double sac(Update * updt, double emax);
//! Updates cell assignments for atoms with unique_ids between istr, and iend
        int updatecells(int istr, int iend);
        void mark_moved_atoms(int i, int j, int k, int l);
        void restore_markings(int i, int j);
        void reset_markings(int i, int j);
        /**
         * \brief Partial evaluation of excluded volume.
         *
         * Calculates the contribution to excluded volume by atoms istrt..iend because
         * of their interactions among themselves and all other atoms. Aborts when
         * value exceeds an upper-limit etest. etest can be calculated in advance.
         * For every update, all other energy terms are calculated before ExVol. The
         * random number used in the Monte Carlo question about acceptance is also
         * generated before ExVol. With this one can calculate what is the maximum
         * value the contribution from this range of atoms can be, so that the event
         * would be accepted. If the contribution exceeds this, the event and the
         * new conformation will be rejected, so there is no point continuing the
         * excluded volume calculation any further. All contributions to ExVol are
         * positive definite, which allows for this simplification.
         */
        double partial_esa(int istrt, int iend, double etest);
//! Finds neighbours for atom i
        void select_neighbours(int i);
        double calc_esa(double);
//! Restores cell structure for atoms istr to iend in case an event is rejected
        int restorecells(int istr, int iend);
        int checksubsets(Update * updt);

//! Implements the functional form of excluded volume
        /**
         * Both Vexva and Vexvb implement the functional form of excluded volume.
         * Vexva is used for all atoms which are not located in a cell at the box
         * boundary. Vexvb is used for atoms in cells at the box boundary, for which
         * part of the neighbouring cells may lie on the opposite side or cornor of the
         * periodic box.
         */

        inline double Vexva(int i, int j) {
            double r2,r6;
            int a;

            if ((r2=dist2a(i,j))>cutg2) return 0;

            if ((a=PairType(i,j))<0) return 0;

            r6=sig2[a]/r2;

            r6*=r6*r6;

            return r6*r6+asa[a]+bsa[a]*r2;
        }

//! Implements the functional form of excluded volume
        inline double Vexvb(int i, int j) {
            double r2,r6;
            int a;

            if ((r2=dist2b(i,j))>cutg2) return 0;

            if ((a=PairType(i,j))<0) return 0;

            r6=sig2[a]/r2;

            r6*=r6*r6;

            return r6*r6+asa[a]+bsa[a]*r2;
        }

    };

}

#endif
