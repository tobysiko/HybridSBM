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

#ifndef ProgUtils_HH
#define ProgUtils_HH
#include <string>
#include <vector>
#include <deque>
#include <map>
#include <list>
#include "InstructionString.hh"

namespace prf_utils
{
    //! Base class for program switches and options
    class OptBase
    {
    public:
        OptBase();
        virtual ~OptBase();
        OptBase(std::string longnm, std::string shortnm, std::string shorthlp);
        OptBase(const OptBase &ob);
        OptBase &operator=(const OptBase &ob);
        //! Set a name for an option created with the default constructor
        inline void Name(std::string st) {myname=st;}
        //! Retrieve the option name.
        inline std::string Name() const {return myname;}

        //! Retrieve long name
        inline std::string LongName() const {return longname;}
        //! Set a long name
        inline void LongName(std::string st) {longname=st;}

        //! Retrieve short help-text/info
        inline std::string short_help() const {return shhlp;}
        //! Set short help text
        inline void short_help(std::string st) {shhlp=st;}

        //! Help text
        inline std::string help() const { return hlp; }
        //! Set help text
        inline void set_help_text(std::string st) { hlp=st; }
        //! Append to existing help text
        inline void append_help_text(std::string st) { hlp+=st; }

        //! Activate the option
        /**
        * Just because a program accepts 15 options does not mean that
        one would need to specify all of them every time one runs it!
        So, an option is created inactive. Only when the user specifies
        it on does one "activate" it. Activated or not-
        activated describes states of an option. But the flag is not
        used in this class for any internal purpose.
        */
        inline void set_given() {given=true;}

        //! Is the option active ?
        inline bool state() const {return given;}
    protected:
        bool given;
        std::string myname,longname,shhlp,hlp;
    };

    //! A program switch
    /**
      * It's a boolean parameter altering the program behaviour. If it
      * is declared as a ProgSwitch, instead of a plain bool variable,
      * it can be used in connection with ProgArgs and will be automatically
      * interfaced through settings files and the command line. A program
      * switch can be given on the command line in positive or negative
      * form. For instance, a flag called "green" makes something happen.
      * It can be given on the command line as "--green". The explicit
      * negation of the flag can be given as "--no-green". In a settings
      * file, these two take the form "green on " and "green off".
      */
    class ProgSwitch : public OptBase
    {
    public:
        ProgSwitch();
        ~ProgSwitch();
        //! Constructor with long and short names, default state and short help
        ProgSwitch(std::string longname,std::string shortname, bool defstate,
                              std::string somedesc="");
        ProgSwitch(const ProgSwitch &ps);
        ProgSwitch &operator=(const ProgSwitch &ps);
        //! Query if the switch is on
        inline bool on() const { return val; }
        //! Query if the switch is off
        inline bool off() const { return !val; }
        //! Turn the switch on
        inline void turn_on() { val=true; }
        //! Turn the switch off
        inline void turn_off() { val=false; }
    private:
        bool val;
    };

    //! One option that a program may handle
    /**
    * This class is mostly intended to be used along with ProgArgs. It
    * represents one of the many options which might be handled by
    * ProgArgs.
    * \ingroup utilities
    */

    class ProgOpt : public OptBase
    {
    public:
        ProgOpt();
        //! Construct with name, number of arguments and a description
        /**
         This is likely the most useful way of constructing these. A
        long name is a descriptive full name of the option. On the
        command line one would specify an option with the long name by
        preceding it with two hyphens. Like, for an option with a long
        name "output_file", one would write, for example, <br>
         a.out --output_file somefile.dat <br>
         Short name is the abbreviation for the same option, specified
        on hte command line with a single hyphen. The integer argument
        numargus specifies how many values are passed to this option.
        Like, one value for the output_file option above. The somedesc
        argument is a some short information about the option. It could,
        for instance, specify that an integer argument actually can take
        only 3 legitimate values.
        */
        ProgOpt(std::string longname,std::string shortname,int numargus=1,
                std::string somedesc="");
        //! Copy constructor
        ProgOpt(const ProgOpt &);
        ~ProgOpt();

        ProgOpt &operator=(const ProgOpt &po);

        //! Set/change the number of expected arguments.
        inline void nargs(int i) {na=i;vals.resize(na,"");}

        //! The i'th value passed to this option.
        inline std::string value(int i) {return vals[i];}

        //! Set the i'th value for this option.

        inline void value(int i, std::string st) {if (i<(int)vals.size())vals[i]=st;}


        //! Retrieve number of expected arguments.
        inline int max_args() const {return na;}
    public:
        int na;
        std::vector<std::string> vals;
    };

    //! A utility to help manage program parameters
    /**
    * \ingroup utilities
    This utility class handles command line and settings file parameters
    for PROFASI programs. It converts all parameters passed through
    the command line, primary and secondary settings files, into a list
    of InstructionString objects.

    An application program would need to define a ProgArgs
    object, tell it what options the application will accept,  and
    how many arguments each of these options should expect. Then
    the application can ask  the ProgArgs object to parse the command
    line and the settings files. Once this is done, the entire set
    of user requests can be obtained as a list of InstructionString
    objects using the function "get_options()". The application
    can then process the InstructionStrings in as appropriate.

    On the command line, a string passed is an option if it starts
    with "-", unless it is a  number. Options can be specified using
    long names (like "--output_file_name"), or with a  short name
    ("-o").  ProgArgs does not discard any arguments found in the
    command line that it does not  recognize
    as an option. All such arguments are stored as "spare_args", and
    can be retrieved as such by the application. If a program uses
    ProgArgs to parse its command line options, they can be passed
    in any order.

    All options of a program can also be set using a settings file.
    The command line option "--settings-file" or "-st" can be used
    to change the name of the settings file. The default is
    "settings.cnf". This option, if present, is executed immediately.
    The other command line options are first collected and stored in
    a list. Settings file options are appended to that list.

    The main program decides the priority of command line vs settings
    file. This class provides methods to parse command line and to
    read a settings file. Which ever comes last will have priority,
    in case a command appears twice. For some options, such as
    "add_chain" for programs like BasicMCRun, a second instance of
    the option does not replace the first, but simply adds one more
    chain.

    The old ProFASi way of using ProgArgs was to access what value a user
    passed as a command line argument like in the following example.
    <br>
    \verbatim
    ProgArgs opts;
    opts.option("output_file_name","o", 1);
    opts.analyze(argc,argv);

    std::string ofile="generic.output";
    if (opts.option_given("o")) ofile=opts.option("o");
    \endverbatim

    This will still work. But this will not retrieve instructions
    from the settings files. Avoid this. Do the following instead.
    <br>
    \verbatim
    ProgArgs opts;
    opts.option("output_file_name","o", 1);
    opts.init_options(argc,argv);
    std::string ofile="generic.output";
    std::list<InstructionString> cmd=opts.get_options();
    for (std::list<InstructionString>::iterator it=cmd.begin();
        it!=cmd.end();++it) {
        if (it->head()=="output_file_name") ofile=it->tail().str();
    }
    \endverbatim
    Although the above seems longer than the old usage shown earlier,
    the new usage parses both the command line and settings file and
    deals with multiple instances of the same command in a sensible
    way. The developer needs to write only one command handling
    function taking care of the InstructionString, and not a separate
    one for command line and settings file. The end-user needs to
    remember only one set of options which work both with the
    settings file and the command line.

    This  class may some day be removed from PROFASI, and replaced with
    an alternative from a standard package (such as boost) providing
    the same functionality. For the present these few hundred lines
    of code help keep PROFASI self contained.
    */

    class ProgArgs
    {
    public:
        //! Default constructor
        ProgArgs();
        ~ProgArgs();
        //! Declare a new option with some specified properties
        /**
        When a new option is added for handling, it has to be
        identified with a long and a short name. The number of
        arguments to be expected after this option (for command line
        input) must be specified, and optionally a small info text
        clarifying the option when needed. For instance, an option
        might expect only integers between 0 and 5. The help text
        could then be "{0--5}".

        A special case arises when we extend the command line
        parser to process what used to be handled exclusively
        in a settings file. Settings file commands such as add_chain
        or new_obs take an indefinite number of arguments:

        \verbatim
        add_chain 1 < ACE * NYSDFRIKLK * NH2 >
        new_obs ProteinRMSD rmsd using +BB ; struc1 abc.pdb:5:A,2,11 ; struc2 $::A
        new_obs Rg rg
        \endverbatim

        To set up such commands for ProgArgs, pretend that they
        are commands which take one or two string arguments.
        \verbatim
        ProgArgs opts;
        opts.option("add_chain","add_chain", 2,"(number of chains and sequence)");
        opts.option("new_obs","new_obs", 1, "(add a new observable)");
        opts.init_options(argc,argv);
        \endverbatim

        This way, add chain will expect two arguments and new_obs only 1. They
        can be used as follows:
        \verbatim
        $ some_program --add_chain 1 "< ACE * NZSDFRIKLK * NH2 >" --new_obs "Rg rg"
        \endverbatim

        Remember that every command line option can now be written in the
        settings file, without the "--" and the quotation marks in the
        settings file.
        */
        void option(std::string longname,std::string shortname, int nargus,
                    std::string helptext="");
        //! Create a new program switch
        /**
          * A switch is like an option, but it does not take any arguments.
          * It may be turned on or off on the command line or the settings
          * file for some effect on the program. For example, imagine that
          * a switch called "verbose", or "v" is created with this function.
          * In the program, one can now use the command line option "-v",
          * to have the intended effect of "verbose". If the program is
          * so written that the default mode is verbose, one would be able
          * to turn this switch off by using the option "--no-verbose". The
          * role of the function parameter "defstate" is to specify if the
          * switch is on or off by default. This default value is used,
          * when the user does not explicitly pass it on the command line.
          */
        void new_switch(std::string longname, std::string shortname,
                        bool defstate=true, std::string helptext="");

        //! Get options from command line and settings files
        /**
          This function first acquires the command line options. If the special
          option "--settings_file" is given, that instruction is executed
          immediately, to change the name of the settings file. After that,
          instructions are acquired from the primary and secondary settings
          files. The command line instructions are appended at the end, so that
          there is one coherent list of instructions which can be executed to
          set up the application.

          \anchor secondary_settings

          Under the hood, some additional work is done to speed up processing
          of settings files. When a parallel application has to read the
          settings file, it is inefficient to have each and every process
          separately open and close the file. Now, only the rank 0 opens the
          primary settings file, and broadcasts the contents to the other
          ranks. In case there is a "secondary_settings on" instruction,
          each process proceeds to open its secondary settings file. There is
          also a new way to provide rank specific instructions in the primary
          settings file, so that the secondary settings is not really needed.
          For example, to set log_level to 100 on rank 0 and to 3 on every
          other rank, you would write this in the primary settings file:

          \verbatim
          log_level 3
          for_rank 0 log_level 3
          \endverbatim
          The rank specific commands are executed after the rank independent
          settings. There is also a "for_rank_range i j" version for the above
          command. The advantage of this method over the "secondary_settings"
          mechanism is that only one settings file is read by one process during
          start up. This will be relevant for people running ProFASi simulations
          with tens of thousands of MPI processes on large supercomputers.
          */
        void init_options(int argc, char *argv[]);

        //! List of all instructions in command line and settings files
        /**
         This returns the combined list of InstructionStrings, possibly
         containing both command line arguments and settings file
         instructions.
          */
        inline std::list<InstructionString> & get_options() { return cmds; }

        //! Analyze the command line to find given options.
        /**
         This function creates a list of options specified on the command line
          */
        void analyze(int argc,char *argv[]);

        //! Number of left over arguments on the command line
        /**
        Are there any arguments passed on the command line that
        were not recognized as any of the known options or any
        values processed along with those options ? Such command
        line arguments are stored separately as "spare arguments".
        This function returns the number of such arguments.
        */
        int n_spare_args() const {return remn.size();}

        //! i'th left over argument
        std::string spare_args(int i);

        //! List available options
        /**
        Lists them along with the expected number of arguments and
        the corresponding short help texts.
        */
        void write_available();

        //! Disable an option even if it is recognized
        /**
        This is useful if you want to inherit from a class which has a
        ProgArgs object, and want to edit the options supported. Adding
        new options is never a problem. But some of the options of the
        base class might not make sense for the derived class. They can
        then be "disabled". A disabled option can be enabled later. If
        such an option was meant to take arguments, it still will. But
        the values of the arguments will never be used. The function
        option_given always returns false for a disabled option.
        */
        void disable(std::string opnname);

        //! Enable an option explicitly, after having disabled it.
        void enable(std::string opnname);

        //! Set rank, if you need to parse secondary settings
        inline void set_rank(int i) { rank=i; }
        //! Assign a different settings file name
        inline void settings_file_name(std::string flnm) { settingsfile=flnm; }
        //! Choose whether to use settings file
        /**
          * One can choose whether the application uses the settings file
          * parsing features provided here. Calling this function with
          * argument 0 turns off settings file parsing altogether. If it
          * is 1, the program will always try to find a settings file, and
          * use it. A value 2 here means the use of a settings file is
          * conditional. By default no settings file will be looked for.
          * But if the user explicitly provides a settings file with the
          * "-st" option and that file exists, it will be parsed.
          */
        inline void settings_use_type(int i) {
            settings_use=i;
            if (i==0) disable("settings_file"); else enable("settings_file");
        }

        inline void clear_cache() { cmds.clear(); }

        //! Convert command line arguments into InstructionStrings
        /**
         Convert whatever is passed on the command line into a lsit of
         InstructionStrings and add to the queue. This does not preserve
         the order in which the commands are passed. This is help with
         the flexibility that command line arguments can be passed in
         any order. Therefore, if you want commands to be executed in
         one particular order, use the settings file.

         Also, by choosing when to call this function, the main program
         can set the priority of command line arguments relative to the
         settings file. Since (normally) one would execute the commands
         in the queue in sequence, the ones that come later will have the
         "last word". If this function is called before "get_settings()"
         the commands in the settings file will overwrite whatever was
         set with the command line. If called after get_settings(), the
         command line arguments will have priority. Normally, you should
         not have to call this function directly at all. Just be happy
         with the sequence "opts.init_opts(argc,argv);
         cmdlist=opts.get_options();", unless you really can't do without
         fine grained control over the order of command line and settings
         file instructions.
          */
        void add_cmd_line_instructions();

        //! Get instructions from the settings files
        int get_settings();
        int get_primary_settings();
        int get_secondary_settings();

        //! Is option optname given on the command line ?
        bool option_given(std::string optname);
        //! Is an option with this name currently available ?
        bool option_available(std::string optname);

        //! Value passed for option "optname" on the command line
        std::string option(std::string optname);
        //! i'th value passed for option "optname" on the command line
        std::string option_arr(std::string, int);

        //! This is to query the state of the switch
        /**
          * If the user specified the switch through the command line or
          * the settings file, the value specified will be returned.
          * If the user didn't specify anything, the default value for
          * the respective switch will be returned.
          */
        bool state_of_switch(std::string swnm);
        //! This is to query if a switch was specified by the user
        /**
          * This just returns whether the user specified a value. Not
          * what value.
          */
        bool switch_given(std::string swnm);

    private:
        void set_availability(std::string opnm, bool avl);
        std::deque<ProgOpt> opts;
        std::deque<ProgSwitch> swtchs;
        std::deque<std::string> remn;
        std::map<std::string,bool> available;
        std::list<InstructionString> cmds,cmdlni;
        std::string settingsfile;
        int rank,settings_use;
    };
}

#endif
