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

#ifndef Contact_HH
#define Contact_HH

namespace prf
{
    //! An abstraction of a contact
    /**
     * This represents a link or a bond between two objects which are
     * represented by two integers. It could, for instance, be an "edge" of
     * graph theory, between two "nodes" labeled as integers. An object of
     * this type merely stores those integers. It can be used to store
     * information about one hydrogen bond or one native contact.
     */

    class Contact
    {
    public:
        Contact();
        ~Contact();
        //! Create a Contact object, between integers i1 and i2
        Contact(int i1,int i2);
        //! Copy a contact
        Contact(const Contact &c);

        //! Assignment
        Contact &operator=(const Contact &c);

        //! Comparison to see if two contact objects represent the same contact
        inline bool operator==(const Contact &c) {
            return(first==c.first&&second==c.second);
        }

        //! "Less than" operation is useful for sorting
        inline bool operator<(const Contact &c) {
            return(first==c.first?second<c.second:first<c.first);
        }

        inline bool operator>(const Contact &c) {
            return(first==c.first?second>c.second:first>c.first);
        }

        inline bool operator>=(const Contact &c) {
            return(first==c.first?second>=c.second:first>=c.first);
        }

        inline bool operator<=(const Contact &c) {
            return(first==c.first?second<=c.second:first<=c.first);
        }

        int first,second;
    };

}

#endif
