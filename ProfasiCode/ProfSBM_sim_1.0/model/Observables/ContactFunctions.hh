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

#ifndef ContactFunctions_HH
#define ContactFunctions_HH
#include "Contact.hh"
#include "../Energy/HydrogenBond.hh"
#include "../Energy/Hydrophobicity.hh"


namespace prf
{
    //! A boolean function object class on Contacts
    /**
     * A contact function is a boolean function, taking either a Contact
     * or two integers representing a contact as argument, that returns true
     * if there is really a contact between those points. For instance, let's
     * consider a graph with three points {1,2,3}. Let there be edges between
     * (1,2), and (1,3), but not between (2,3). A contact function describing
     * this system should then return true for Contact objects (1,2) and (1,3)
     * while returning false for (2,3). The points 1,2,3 could be 3 amino
     * acids, and we may be interested in whether two amino acids are close
     * in space. Then the contact function represents the concept that this
     * condition is satisfied. It is true if the given amino acids are close,
     * false otherwise.
     *
     * With increasing number of different types of contact functions, such as
     * CaContact, HPContact, Proximity, HBContact, it seems useful to introduce
     * a convention for what the two integers stored in a Contact mean. With
     * PROFASI v. 1.1, for ContactFunction objects representing contacts between
     * residues or end-groups, these integers represent the ligand indices in
     * the entire population. For instance, if the population consists of
     * 3 copies of the peptide &lt; ACE * KLVFFAE * NH2 &gt;, a contact between
     * the K and A of the second chain would be denoted with integers (10,15).
     * If the population consists of 3 copies of &lt;  * KLVFFAE * &gt;, a contact
     * between the same residues will be denoted by (7,12).
     */

    class ContactFunction
    {
    public:
        ContactFunction();
        virtual ~ContactFunction();
        inline void min_seq_sep(int i) {minseqsep=i;}
        inline bool is_symmetric() {return issym; }

        virtual int init(Population *p);
        virtual void set_cutoff(double xx);
        virtual double get_cutoff() const;
        // This base class version (always true) is expected to be overridden
        virtual bool operator()(int arg1,int arg2) const ;
        virtual bool operator()(const Contact & c) const ;
        virtual double operator()(int arg1,int arg2, bool dummy) const ; //NEW! -TS
    protected:
        int minseqsep;
        bool issym;
    };

    //! A simple contact function based on Calpha separation
    /**
     * This contact function merely checks if the Calpha atoms of two
     * residues are within a certain given distance of each other. If they
     * are, the residues are regarded as being in contact. The default
     * for the distance cut-off is 6.0 Angstroms. It can be changed
     * with the function set_cutoff(double).
     */

    class CaContact : public ContactFunction
    {

    public:
        CaContact();
        CaContact(double dtc);
        int init(Population *p);
        bool operator()(const Contact &c) const ;
        bool operator()(int i1,int i2) const ;
        double operator()(int i1, int i2, bool dummy) const ; //NEW! -TS
        //! Set and get cutoff for C-alpha contacts
        void set_cutoff(double xx);
        inline void setCut(double xx) {distcut=xx;distcut2=distcut*distcut;}

        double get_cutoff() const;
    private:
        std::vector<int> indca;
        double distcut,distcut2;
    };

    //! A ContactFunction based on spatial proximity
    /**
    * This contact function returns true if the two concerned amino acids
    * are close in space. The closeness is defined in this way: If a heavy
    * atom in amino acid i is closer than a certain given cutoff distance
    * to any heavy atom in amino acid j, there is said to be a link between
    * the two residues. The residues are said to be in contact if there is
    * at least two (or some other number specified by the minimum_links(int)
    * function) links between them. The cutoff distance is tunable. This
    * contact function needs an initialization with a Population * as
    * the argument.
    */

    class Proximity : public ContactFunction
    {

    public:
        Proximity();
        ~Proximity();

        //! Minimum number of close pairs of heavy atoms for a contact, default=2
        /**
         If you set minimum links to 1, this contact function will work like a
        simple heavy atom contact function.
        */
        inline void minimum_links(int nlnk) {min_links=nlnk;}

        inline void distCutoff(double xx) {natCut=xx;natCut2=natCut*natCut;}

        void set_cutoff(double xx);
        double get_cutoff() const;
        int init(Population *p);
        bool operator()(const Contact &c) const ;
        bool operator()(int i,int j) const ;
    private:
        std::vector<std::vector<int> > hatoms;
        std::vector<int> natoms,nhatoms;
        double natCut,natCut2;
        int min_links;
    };

    //! A ContactFunction based on hydrogen bonds
    /**
     * This contact function checks if there is a hydrogen bond between two
     * given amino acids, represented by the contact. That is, for a given
     * Contact object with elements i and j, it returns true if there is a
     * hydrogen bond with a donor (NH) on i and an acceptor (CO) on j. The
     * cutoff to decide when to regard a hydrogen bond as formed is tunable.
     * The function looks only at the backbone backbone hydrogen bonds. The
     * indices i and j are to be interpreted as ligand indices relative to
     * the entire population, starting from 0. This is for consistency with
     * other ContactFunctions.
     */

    class HBContact : public ContactFunction
    {

    public:
        HBContact();
        HBContact(HBMM *gohbmm,double hbc=1.03);
        ~HBContact();
        int init(Population *p);
        bool operator()(const Contact &c) const;
        bool operator()(int i1,int i2) const;
        inline void connectHB(HBMM *hbm) {ohbmm=hbm;}

        inline void setHBCut(double xx) {hbcut=xx;}

        void set_cutoff(double xx);
        double get_cutoff() const;
    private:
        HBMM *ohbmm;
        double hbcut;
        bool myhb;
    };

    typedef HBContact HBCheck;

    //! A Contact function for Hydrophobic contacts
    /**
     * This contact function checks if the given residues are in
     * hydrophobic contact. It is not based on the hydrophobic energy,
     * but simply on how many hydrophobic atoms of one non-polar
     * residue are in contact with one or more such atoms of the other.
     * This number is then scaled by its maximum value to obtain a
     * fraction between 0 and 1. When this number is greater than a
     * certain lower cutoff, by default set to 0.1, the residues are
     * regarded as being in contact. The minimum contact fraction can
     * of course be changed with the function set_cutoff(double).
     *
     * The hydrophobic contact makes sense for amino acids as well as
     * end groups.
     */

    class HPContact : public ContactFunction
    {

    public:
        HPContact();
        HPContact(Hydrophobicity *gohpf,double ftc=0.1);
        ~HPContact();
        inline void connectHP(Hydrophobicity *hpf) {ohpf=hpf;}

        int init(Population *p);
        bool operator()(const Contact &c);
        bool operator()(const Contact &c) const ;
        bool operator()(int i1,int i2) const ;
        inline void setCut(double xx) {fraccut=xx;}

        //! Change contact fraction cutoff for hydrophobic contacts
        void set_cutoff(double xx);
        double get_cutoff() const;

    private:
        Hydrophobicity *ohpf;
        double fraccut;
        bool myhp;
    };

    class HBContactChains : public ContactFunction
    {

    public:
        HBContactChains();
        HBContactChains(HBMM *gohbmm,double hbc=1.03);
        ~HBContactChains();
        int init(Population *p);
        bool operator()(const Contact &c) const;
        bool operator()(int i1,int i2) const;
        inline void connectHB(HBMM *hbm) {ohbmm=hbm;}

        inline void setHBCut(double xx) {hbcut=xx;}

        void set_cutoff(double xx);
        double get_cutoff() const;
    private:
        HBMM *ohbmm;
        double hbcut;
        bool myhb;
    };
}

#endif
