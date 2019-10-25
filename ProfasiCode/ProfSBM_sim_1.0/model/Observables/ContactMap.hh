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

#ifndef ContactMap_HH
#define ContactMap_HH

#ifndef PRFTARGETSIZE
#define PRFTARGETSIZE 32
#endif
#ifndef PRFMAXCONTACTS
#define PRFMAXCONTACTS 65536
#endif

#include <bitset>
#include <vector>
#include <string>
#include "Observable.hh"
#include "ContactFunctions.hh"

namespace prf
{
    //! Utility to help keep track of a set of contacts
    /**
     * Contact map is an abstract class to keep track of a set of "Contacts".
     * A contact might be a hydrogen bond, or hydrophobic contact, or mere
     * proximity ... Typically we will have a function bool inContact(int,int),
     * taking two integer arguments representing objects and telling us if they
     * are in contact, like two amino acids. The function may or may not be
     * symmetric.
     *
     * What we want to do is to read in contact information from a file, and
     * make a set of Contacts we want to monitor. At any given point of time,
     * one might want to scan over the set of specified contacts to see which
     * ones are on, which ones are off. This information can then be
     * conveniently output as one or more integers, whose bits would represent
     * individual contacts. It should also be possible to invert this set of
     * integers to recreate the contacts they represent.
     *
     * The origin of this class was a project in which such an enermous amount
     * of contact information had to be saved that size of the files on disk
     * became an issue. A contact, on the other hand, is either present or
     * absent. So, it can be stored in a bit. This allows one to store
     * the state of 32 different contacts in one single 32-bit integer.
     *
     * As an Observable, its value is the number of intact contacts. The
     * function refresh is left as virtual so that it can be over-ridden
     * again in subsequent inheritted classes.
     * \ingroup profasi_observables
     */

    class ContactMap : public Observable
    {
    public:
        ContactMap();
        virtual ~ContactMap();
        double evaluate();
        int init_obs();
        void write_snapshot(Output &op);
        void write_rtkey(Output &op);
        void avg_reset();
        void avg_fill(int itmp);
        void avg_write(Output &op);
        void WriteMaps();
        void init();
        void ReadMap(std::string mapfile);
        void track_all_contacts();

        inline void setContactFunction(ContactFunction *cf)
        {contact_function=cf;}

        inline int numContacts() const {return n_monitored_contacts;}

        inline int repSize() const {return representation_size;}

        void repSize(int i);
        inline int Count() {return mainmap.count();}


        void rangeEstimate(double &x0,double &x1);

        inline bool StateOf(int i) {return (bool) mainmap[i];}

        inline bool Intact(int i) {return (bool) mainmap[i];}

        inline bool Broken(int i) {return ~(bool) mainmap[i];}

        inline std::bitset<PRFMAXCONTACTS> MainMap() const {return mainmap;}

        /* The following functions should only be used when the number of
           contacts monitored is less than the number of bits in unsigned
           long. Else an exception will be thrown. Remember that PRFMAXCONTACTS
           is a larger number. */

        inline unsigned long State() {return mainmap.to_ulong();}

        inline void State(unsigned long val) {
            mainmap=std::bitset<PRFMAXCONTACTS>(val);
        }

        //If the number of contacts monitored is larger, use the follwing set
        // of functions.
        void DivideToGroups();
        inline unsigned long StateGroup(int i=0) {return cmap[i].to_ulong();}

        inline void StateGroup(int i, unsigned long val)
        {cmap[i]=std::bitset<PRFTARGETSIZE>(val);}

        void CombineGroups();
        inline Contact getContact(int i) {return contact[i];}

        inline int nTotNatBond() const {return numContacts();}

        inline Matrix<double> & ContactFrequencyTable() {return prof;}

        // Overwrite base class function to not let the range and nbins change
        void his_range(double xmn, double xmx);
        void his_nbins(int n);
    protected:
        void parse_CF_opts(std::vector<std::string> &ops, ContactFunction *f);
        std::bitset<PRFMAXCONTACTS> mainmap;
        std::vector<std::bitset<PRFTARGETSIZE> > cmap;
        int n_monitored_contacts,representation_size;
        std::vector<Contact> contact;
        std::vector<int> dndt;
        Matrix<double> prof;
        ContactFunction *contact_function;
        bool mycf;
    };
}

#endif
