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

#ifndef Named_HH
#define Named_HH
#include <string>

namespace prf
{
    //! Anything that has a name
    /**
     * This little class does nothing other than provide an interface to
     * assign and retrieve the name of an object.
     */

    class Named
    {
    public:
        //! Create an object with name "unnamed"
        Named();
        //! Create an object with name st given as a C string
        Named(const char *st);
        //! Create an object with a name given as a string
        Named(std::string st);
        virtual ~Named();
        //! Retrive the name of an object
        inline std::string Name() const {return nm;}

        //! Assign a new name to an object
        inline void Name(std::string gnm) {nm=gnm;}

    protected:
        std::string nm;
    };
}

#endif
