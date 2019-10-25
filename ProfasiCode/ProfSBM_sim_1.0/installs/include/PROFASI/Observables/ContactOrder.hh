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

#ifndef ContactOrder_HH
#define ContactOrder_HH
#include "Observable.hh"
#include "ContactFunctions.hh"
#include "../Elements/Population.hh"

//! Relative contact order
/**
* This observable makes sense for individual chains only. It is a
* measure of the complexity of the topology of a protein fold, which
* is useful for proteins large enough to have tertiary structure. The
* definition is as in :<br>
*
* "Contact Order, Transition State Placement and the Refolding Rates
* of Single Domain Proteins", Kevin W. Plaxco, Kim T. Simons and David Baker,
* J. Mol. Biol. (1998) 277, 985-994<br>
*
* Qualitatively, contact order is the mean sequence separation between
* residues in physical contact, normalized by the chain length.
* A purely helical structure will have relatively small contact order,
* as the residues in contact are close in sequence. In a &#946;-sheet,
* the  contacting residues could be  very far in sequence, leading to
* high contact orders. This class returns the contact order in percentage,
* which is conventional.
*
* For the evaluation of contact order, the following definition of a
* contact is used: two residues are in contact if any non-hydrogen
* atom of one is within D Angstroms of a non-hydrogen atom of the
* other. The default value of cutoff D is 6, but can be changed. There
* is no minimum sequence separation for a contact for this observable.
* So, neighbouring residues are always in contact, and contribute to the
* average.
*
* As a note on using this class, we point out that for a 100 residue
* protein with 1000 heavy atoms, in the worst case scenario, about
* half a million distances would have to be evaluated. This would be
* very expensive. It is therefore not recommended to use this Observable
* as one of the very frequently evaluated variable, like for instance,
* secondary structure content. Using the "cell" information from
* the excluded volume energy term can speed it up. But that is not
* available to Observable classes at the moment. Until further notice,
* we recommend that this class is only used for post run analysis
* programs.
* \ingroup profasi_observables
*/

class ContactOrder : public Observable
{

public:
    ContactOrder();
    ~ContactOrder();
    //! Which chain. Contact order only works on one chain. Default is 0.
    inline void of_chain(int i) {ich=i;}

    //! Which chain is currently in use
    inline int of_chain() const {return ich;}

    //! Currently used cutoff to define contacts
    inline double get_cutoff() const {return dcut;}

    //! Define a cutoff distance to be used in the definition of a contact
    inline void set_cutoff(double x) {contact_cutoff_distance(x);}

    inline void contact_cutoff_distance(double x) {dcut=x;cf.set_cutoff(dcut);}

    double evaluate();
    int init_obs();
    void rangeEstimate(double &x1, double &x2);

private:
    int ich;
    double dcut;
    CaContact cf;
};

#endif
