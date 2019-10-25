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

#include "Minimizer.hh"
#include "../Aux/profasi_io.hh"

Minimizer::Minimizer() : maxcyc(-1), updttyp(1), ftol(1.0e-7),
goodenough(-1.0e7),knowgoodenough(false),absscale(false) {}

Minimizer::~Minimizer() {}

void Minimizer::set_function(NDFunctional *gf)
{
    if (gf==NULL) return;
    f=gf;
    ncgreset=ndof=gf->n_dim();
    if (maxcyc<0) maxcyc=10*ndof;
    lp.set_base_function(f);
}

double Minimizer::minimize()
{
    prf::Logger blog(10);
    std::valarray<double> x;
    x.resize(ndof,0);
    f->get_position(x);
    const double eps=1.0e-7;
    double fp=f->value();
    double fret=fp;
    std::valarray<double> g(0.0,ndof),h(0.0,ndof),xi(0.0,ndof);
    double dgg, gam, gg, lmstpsize=0.01;
    f->gradient(xi);

    for (int i=0;i<ndof;++i) {
        xi[i]=h[i]=g[i]=-xi[i];
    }
    blog<<"Conjugate gradient minimization : starting function value is "
            <<fp<<"\n";
    for (int i=0;i<maxcyc;++i) {
        lp.set_origin(x);
        lp.set_direction(xi);
        double t=0;
        fret=lp.nearest_minimum(t,lmstpsize);
        blog(20)<<"Minimization cycle "<<i<<" found a minimum value "
                <<fret<<"\n";
        lp.get_last_position(x);
        f->set_position(x);
        if ((absscale && fabs(fret-fp)<ftol)) {
            blog(15)<<"Ending CG minimization loop because we are using "
                    <<"absolute scales and the difference from the value "
                    <<"at the end of the last iteration "<<(fret-fp)
                    <<" is smaller than the tolerance "<<ftol<<"\n";
            break;
        } else if ((!absscale)&&2*fabs(fret-fp)<ftol*(fabs(fret)+fabs(fp)+eps)) {
            blog(15)<<"Ending CG minimization loop because the minimum  "
                    <<fret<<" obtained in the current iteration is close to "
                    <<fp<<", where the last iteration ended.\n";
            break;
        } else if (knowgoodenough and fret<goodenough) {
            prf::cout<<"Ending CG minimization loop because the latest "
                    <<"minimum value "<<fret<<" is below the \"good enough\""
                    <<"threshold "<<goodenough<<"\n";
            break;
        }
        fp=f->value();
        f->gradient(xi);
        dgg=gg=0;
        if (updttyp==0) {
            //Fletcher-Reeves
            for (int j=0;j<ndof;++j) {
                gg+=g[j]*g[j];
                dgg+=xi[j]*xi[j];
            }
        } else {
            //Polak-Ribiere
            for (int j=0;j<ndof;++j) {
                gg+=g[j]*g[j];
                dgg+=(xi[j]+g[j])*xi[j];
            }
        }
        if (gg==0) {
            blog(15)<<"Ending CG minimization loop because gradient is 0\n";
            blog<<"Final function value is "<<fp<<"\n";
            break;
        }
        gam=dgg/gg;
        if ((i+1)%ncgreset==0) gam=0;
        for (int j=0;j<ndof;++j) {
            g[j]=-xi[j];
            h[j]=g[j]+gam*h[j];
            xi[j]=h[j];
        }
    }
    blog(10)<<"Conjugate gradient minimization finished with function value "
            <<fret<<"\n";
    return fret;
}
