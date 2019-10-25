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

#ifndef PRF_RAW_CONF_HH
#define PRF_RAW_CONF_HH
#include "TrajSnapshot.hh"
#include "BinStream.hh"
#include "prf_xml.hh"
#include <fstream>
#include <deque>

namespace prf_traj
{
    //! Reader for PROFASI's raw binary data format
    /**
      This class interprets the binary format in which PROFASI usually
      stores its trajectory information. It is not meant to be instantiated
      directly but only as a back-end for prf_traj::TrajSeg class, which is
      created in the background by the prf_traj::Trajectory class. This class
      also serves as a "template" on which readers for alternative binary
      formats should be written, so that they are directly usable as
      interpreters for parts of PROFASI trajectories.

      \ingroup prf_trajectory
      \sa prf_traj::Trajectory, prf_traj::TrajSeg
      */
    class prf_raw_conf
    {
    public:
        prf_raw_conf();
        ~prf_raw_conf();
        prf_raw_conf(const prf_raw_conf &);
        prf_raw_conf &operator=(const prf_raw_conf &);
        //! Set data file
        void set_file(std::string st);
        //! Retrieve data file name
        inline std::string get_filename() {return datafl;}
        //! Parse headers and file metadata
        /**
          This function reads the header section of the file and determines
          information about the snapshots stored, such as byte order, size
          of different kinds of numbers in bytes, but also sequences of
          the chains in the system. It then determines the starting and
          ending MC cycles in the file and the number of snapshots stored.
          */
        bool init();
        //! Number of blocks or snapshots in the file
        inline unsigned n_blocks() {return nblk;}
        //! Number of DOFs stored in each snapshots
        inline unsigned n_coordinates() const { return ncrd; }
        //! Retrieve the iblk'th snapshot in the file
        bool get_block(unsigned iblk, TrajSnapshot &snp);
        //! Close file handle but keep all information about the file
        bool sleep();
        //! (Re-)open file and perform initialisation if necessary
        bool wake_up();
        //! Return an XML map of information about all snapshots
        prf_xml::XML_Node *xml_map() { return mymap; }
        //! Set "info" file containing header
        /**
          A separate info file is no longer required in PROFASI since version
          1.4. To interpret older configurations without inline headers it
          normally suffices to look for an info file by the name "conf.info"
          in the same directory as the data file. If a different info file is
          desired for some reason, this function can be used.
          */
        inline void set_info_file(std::string st) {infofl=st;}
    private:
        void init_xml_map();
        bool get_header();
        bool read_inline_header();
        bool read_info_file();
        bool parse_header();
        bool check_sizes();

        // File names and other string info
        std::string datafl, infofl, headertxt, byte_order_of_file;
        // sizes
        unsigned blocksz, nblk, ncrd, header_offset;
        unsigned long nbytes;
        // Number of bytes for different kinds of variables
        std::vector<unsigned int> sizes_file,sizes_host;
        // Byte buffer to hold a chunk of unformatted data from the binary file
        char *buf;
        // A sequence of VarType representing the types of variable was written
        std::deque<prf_utils::VarType> signature;
        // XML map of the population and the data layout
        prf_xml::XML_Node *mymap;
        BinStream bin;
        std::ifstream myfl;
        bool awake,initialized;
    };
}
#endif // PRF_RAW_CONF_HH
