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

#ifndef TRAJECTORY_HH
#define TRAJECTORY_HH
#include "TrajSnapshot.hh"
#include "TrajSeg.hh"
#include "profasi_io.hh"
#include <deque>

/**
  \page traj_files PROFASI trajectory files
    A "trajectory" in PROFASI is a sequence of snapshots stored at regular
    intervals during a Monte Carlo run. Trajectory files are generated by
    PROFASI simulation programs such as BasicMCRun, SimTempRun and ParTempRun.

    \section intro Multi-segment trajectories
    On modern computing systems, the run time for a single job is normally
    limited, quite often to one or a few days. This is not enough for a
    production run of a decent protein system. So, one needs several restarts.
    This is how PROFASI handles such multi-part runs:

    When you start a run, tell the program how much time it has. For instance,
    if you request 24 hours for a parallel job on 64 cores from a PBS queue
    manager, put the following line in the PROFASI settings file for the job:
    \verbatim
    available_time 23:45:0
    \endverbatim

    A lot of things need to be saved when the job ends. So, it is normally
    better to leave a bit of time at the end of job: if you have 24 hours,
    tell PROFASI that you have 23 hours 45 minutes. If the job is not finished
    when the available time expires, the simulation "suspends". To restart
    from where it left, just start the same job again in the same directory.
    PROFASI detects that there was an unfinished run and continues the run
    instead of starting a new one. Imagine you had to restart a run 17 times
    to finish 10 million MC sweeps for a system. These 18 stages of the run will
    be called segments.

    \li Each segment of the run creates a new data file to save snapshots
    during that run
    \li The segment files contain headers with information about the
    layout of the data in them and about the machine where that segment
    of the simulation was made.
    \li A meta-data file called "traj" is created when a run starts for the
    first time. The "traj" file contains the names of the segment files which
    make up the trajectory. It also contains the starting MC cycle, interval
    between successive writes to the trajectory and the total number of
    snapshots in each segment.
    \li In case a run "resumes" from the middle of a previous run and traces
    a different trajectory, the original trajectory file is backed up and
    a new one started containing the part of the older trajectory prior to
    the restart point. Only a new meta data file is created. The two "traj"
    files share the data files of the common segments.
    \li Different segments can run on different machines with different
    architectures. While browsing a trajectory, this detail is not visible.

    To extract a particular snapshot with the extract_snapshot program or
    analyze properties of the trajectory with the extract_props program,
    the "traj" files are used as input. But the "traj" files alone are useless
    since they don't contain the actual data. <div style="color:red"> It is not
    enough to copy the tiny "traj" files back from the cluster or supercomputer
    running the simulations. All conf.data... files must be copied.</div>

    The prf_traj::Trajectory class parses the "traj" files and provides an
    interface as if it was one smooth series of snapshots. Only one file handle
    is kept open corresponding to the segment of the snapshot last accessed.
    This is called the "active" segment. If access is attempted to a snapshot on
    a different segment of the run than the active segment, the active segment
    is closed, and the necessary segment is opened and the correct snapshot
    returned. To the outside function, all these details are hidden. They
    iterate over the trajectory and access elements when they want to.
    \sa \ref prf_trajectory, extract_snapshot, extract_props, \ref traj_gen
  */
/**
    \defgroup prf_trajectory PROFASI trajectory module
    \ingroup utilities
    @brief Interface to trajectory files in PROFASI
    This module  contains a set of classes to manage interpretation of
    trajectories (\ref traj_files) .

    The class prf_traj::Trajectory provides an interface to such a trajectory
    saved from a previous run. Each point along the trajectory is called a
    snapshot, represented by the prf_traj::TrajSnapshot class. One can think of
    a Trajectory object as if it were a very restricted version of the standard
    library list with TrajSnapshot elements. The access is read-only, with
    sequencial iterators given by the TrajIterator class. Functions begin(),
    end() and find() are used to browse trajectories. Behind the scene, the
    Trajectory object handles a lot of complications.

    The segments of the trajectory are internally handled by objects of the
    class prf_traj::TrajSeg . The TrajSeg objects again depend on lower level
    structures to actually interpret the binary data in the segment data files.
    But they insulate the Trajectory class from all such low level details. The
    dirty work of reading binary bytes and interpreting them as energy values or
    coordinate values is done by prf_raw_conf . It uses knowledge of the data
    layout to interpret the bytes. In principle, the prf_raw_conf can be
    replaced by an optional alternative format interpreter, without changing
    much in the Trajectory class. In the future, alternatives to PROFASI's own
    binary format, such as the HDF5 format for binary data storage may be
    preferred.

    \sa \ref traj_files
    */
namespace prf_traj {
    class Trajectory;
    //! Iterator to browse trajectories
    /**
      TrajIterator provides an iterator for Trajectory. One can go back and
      forth along the trajectory with increment and decrement operators. Only
      when the value corresponding the the iterator is accessed with the
      "*" operator or the "->" operator, the trajectory is requested to focus
      on that particular snapshot and retrieve the data.
      \ingroup prf_trajectory
      */
    class TrajIterator
    {
    public:
        TrajIterator();
        TrajIterator(const TrajIterator &t);
        TrajIterator(Trajectory *tr, unsigned long wh);
        ~TrajIterator();
        TrajIterator &operator=(const TrajIterator &);
        bool operator==(const TrajIterator &);
        bool operator!=(const TrajIterator &);
        bool operator<(const TrajIterator &);
        bool operator>(const TrajIterator &);
        bool operator<=(const TrajIterator &);
        bool operator>=(const TrajIterator &);
        int operator+=(int n);
        TrajIterator & operator++();
        TrajIterator operator++(int);
        TrajIterator & operator--();
        TrajIterator operator--(int);
        TrajSnapshot & operator*();
        TrajSnapshot * operator->();
        inline bool good() const {return isgood;}
        void refresh();
    private:
        unsigned long idx;
        bool isgood,needsrefresh;
        TrajSnapshot snpsh;
        Trajectory *mytraj;
    };

    //! PROFASI (multi-segment) trajectories
    /**
      This class provides an STL list like interface to a PROFASI trajectory.
      Each element of this fictitious list is a TrajSnapshot. The list can
      be browsed using a TrajIterator from "begin()" to "end()". There is
      also a "find()" function to locate a snapshot with a particular MC time.

      Internally, this class is a collection of TrajSeg objects representing
      different segments of a run. It parses a "traj" file to decide what
      segments make up the trajectory, and the cycle limits assigned to each
      segment. Typical use:

      \verbatim
      prf_traj::Trajectory traj;
      traj.parse("somedirectory/traj");
      if (traj.init()) {
          for (prf_traj::Trajectory::iterator it=traj.begin();
               it!=traj.end();++it) {
               prf::cout<<"Time = "<<it->MC_time()
                        <<", energy = "<<it->energy()<<"\n";
               //Or do something else with the snapshot represented by
               //(*it)
          }
      }
      \endverbatim
      \ingroup prf_trajectory
      \sa prf_traj::TrajSnapshot, prf_traj::TrajSeg
      */
    class Trajectory
    {
    public:
        typedef TrajIterator iterator;
        Trajectory();
        ~Trajectory();
        inline iterator begin() { return bgptr; }
        inline iterator end() { return ndptr; }
        inline unsigned long min_cycle() const { return cycmin; }
        inline unsigned long max_cycle() const { return cycmax; }
        iterator find(unsigned long mct);
        bool focus(unsigned iblk, prf_traj::TrajSnapshot &snpsh);
        int parse(std::string trjfl);
        int append_list(std::string refdir, std::deque<std::string> &fllist);
        int append_segment(std::string sgfl, unsigned long st, unsigned long iv,
                           unsigned long nbl);
        int append_segment(std::string sgfl);
        int init();
        int auto_adjust();
        int status();
        inline void save() { save(myfile); }
        void save(std::string flnm);
        void save_up_to(unsigned long icyc, std::string flnm);
        inline prf_xml::XML_Node *xml_map() {return mymap;}
    private:
        unsigned seg_of_block(unsigned blkid);
        unsigned seg_of_time(unsigned long mct);
        unsigned curblk,curseg;
        unsigned long cycmin,cycmax;
        iterator bgptr, ndptr;
        std::deque<TrajSeg *> segment;
        prf_xml::XML_Node *mymap;
        std::string myfile;
    };
}

#endif // TRAJECTORY_HH
