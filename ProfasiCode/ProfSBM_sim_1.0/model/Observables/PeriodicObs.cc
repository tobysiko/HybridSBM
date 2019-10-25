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

#include "PeriodicObs.hh"
#include "ObsHandler.hh"
using std::string;
using std::vector;

/**
* \page opt_PeriodicObs PeriodicObs
\section options Available options
<b>of</b> <i>of BBRMSD</i><br>
Evaluate the Observable BBRMSD with periodic box centered on the mesh points,
and return the minimum value. The argument passed to "of" is the user specified
alias for the other Observable. <br>
<b>mesh_size</b><i>mesh_size 5</i><br>
During initialisation, PeriodicObs creates a mesh of 5X5X5 (in this case) points
inside the current periodic box. During evaluation, it translates the box centre
to each of these points, applies periodic boundary conditions and evaluates the
given Observable. The return value is the minimum value found.
\section examples Examples
new_obs PeriodicObs prg of rg; mesh_size 4 <br>
If rg is the alias for one of the Observables tracked by ObsHandler, this creates
a new Observable prg, that is the minimum value of rg when the box centre is moved
to 4X4X4 grid points.
\sa prf::PeriodicObs
*/

namespace prf
{
    PeriodicObs::PeriodicObs(prf::ObsHandler *hdl) : handler(hdl),
            theObs(NULL), mshsz(5) {Name("PeriodicObs");}

    PeriodicObs::~PeriodicObs() {}

    int PeriodicObs::init_obs()
    {
        if (handler==NULL or Observable::init_obs()==0) return 0;
        bkpcrds.resize(p->NumberOfAtoms());

        theObs=NULL;

        vector<string> parts;

        for (size_t i =0;i<usrcmd.size();++i) {
            parts.clear();
            split(usrcmd[i],parts);

            if (parts[0]==string("of")) {
                if (parts.size()<2) {
                    prf::cerr<<Name()<<"> Syntax error in command \"of\". Usage:\n"
                    <<"of valid_observable_alias\n";
                    continue;
                }

                std::string obsname=parts[1];

                if ((theObs=handler->get_obs(obsname))==NULL) {
                    prf::cerr<<Name()<<"> The observable handler could not connect me "
                    <<"to any Observable called "<<obsname<<". Initialisation failed. \n";
                    return 0;
                }
            }

            if (parts[0]==string("mesh_size")) {
                if (parts.size()<2) {
                    prf::cerr<<Name()<<"> Syntax error in command \"mesh_size\". Usage: \n"
                    <<"mesh_size some_small_positive_integer\n";
                } else {
                    set_mesh_size(std::atoi(parts[1].c_str()));

                    if (mshsz<3) {
                        prf::cerr<<Name()<<"You asked for a mesh size of "<<mshsz<<". The "
                        <<"minimum allowed value is 3, which will be used instead of"
                        <<" your value.\n";
                        mshsz=3;
                    }
                }
            }
        }

        if (theObs==NULL) {
            prf::cerr<<Name()<<"> No observable set. Initialisaiton failed. \n";
            return 0;
        }

        prf::Vector3 xcap(1,0,0),ycap(0,1,0),zcap(0,0,1);

        double displegt=(1.0/(mshsz-1))*AtomCoordinates::boxL();
        double offset=0.5*AtomCoordinates::boxL();

        mesh_disp.resize(mshsz*mshsz*mshsz);

        for (int i=0;i<mshsz;++i) {
            for (int j=0;j<mshsz;++j) {
                for (int k=0;k<mshsz;++k) {
                    mesh_disp[mshsz*mshsz*i+mshsz*j+k]=
                        (i*displegt-offset)*xcap+
                        (j*displegt-offset)*ycap+
                        (k*displegt-offset)*zcap;
                }
            }
        }

        Logger()(log_thres)<<Name()<<"> Tracking Observable "<<theObs->Name()

        <<" with "<<mshsz<<" mesh points on x, y and z axes\n";
        return 1;
    }

    double PeriodicObs::evaluate()
    {
        double obvl=theObs->evaluate();
//        prf::cout<<"Original value = "<<obsval<<"\n";

        for (size_t i=0;i<bkpcrds.size();++i) bkpcrds[i]=AtomCoordinates::vec(i);

        for (size_t i=0;i<mesh_disp.size();++i) {
            AtomCoordinates::BlockTranslate(mesh_disp[i],0,AtomCoordinates::numberOfAtoms());
            p->EnforceBC();
            double tmpvl=theObs->evaluate();

            if (tmpvl<obvl) {
                obvl=tmpvl;
//      prf::cout<<"New value "<<obsval<<" for translation "<<mesh_disp[i]<<"\n";
            }

            for (size_t j=0;j<bkpcrds.size();++j)
                AtomCoordinates::vec(j,bkpcrds[j].x(),bkpcrds[j].y(),bkpcrds[j].z());
        }
        return obvl;
    }

    void PeriodicObs::rangeEstimate(double &xmin,double &xmax)
    {
        theObs->rangeEstimate(xmin,xmax);
        if (!userbinsz) xbin0=theObs->his_bin_size();
    }

}
