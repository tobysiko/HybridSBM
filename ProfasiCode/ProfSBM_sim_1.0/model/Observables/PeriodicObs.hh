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

#ifndef PeriodicObs_HH
#define PeriodicObs_HH
#include "Observable.hh"

namespace prf
{
    // Forward declaration of ObsHandler
    class ObsHandler;

    //! Minimize another Observable with respect to the periodoc box
    /**
     * Measurements such as the radius of gyration present a problem for
     * for multi-chain systems. If there are two chains touching each other
     * they might end up near opposite walls of the periodic box. Their
     * centre of mass would be at the centre of the box, so translating to
     * the centre of mass would not change anything. The system would appear
     * have a very large radius of gyration although as far as the force
     * field goes, it is a very compact system. Sometimes, one wants to
     * move the periodic box a bit, so that the chains appear at the middle
     * of the box, and radius of gyration evaluates to a small value. This
     * is the purpose of this class.
     *
     * In this preliminary version, there is no smart algorithm to do this.
     * This class works by searching. It moves the periodic box around on
     * a mesh, and tries to find the position that minimizes something. This
     * something has to be another Observable, known to this class through
     * a pointer.
     *
     * It is worth noting that since this class works with a search, it should
     * be avoided during a simulation, and used only to extract information
     * during a post-run analysis. So, it should be added to the settings
     * file only for the extract_props program.
     *
     * \ingroup profasi_observables
     */

    class PeriodicObs : public Observable
    {
    public:
        explicit PeriodicObs(prf::ObsHandler *hdl);
        ~PeriodicObs();
        int init_obs();
        //! attach to two atoms given by their unique ids.
        inline void set_mesh_size(int i1) {mshsz=i1;}

        double evaluate();

        void rangeEstimate(double &xmin,double &xmax);

    private:
        prf::ObsHandler *handler;
        Observable *theObs;
        int mshsz;
        std::vector<prf::Vector3> mesh_disp, bkpcrds;
    };
}

#endif
