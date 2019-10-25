//Contents of this file are required only for documentation
/**
* \mainpage ProFASi: The Protein Folding and Aggregation Simulator


\section intro_sec Introduction
<div style="text-align:justify">
PROFASI (PROtein Folding and Aggregation SImulator [\ref I]) is a C++ program  package for Monte Carlo simulations of protein folding and aggregation. It provides an implementation of an all-atom protein model with fixed bond lengths and bond angles, an implicit water simplified force field, and a set of tools to perform Monte Carlo simulations with the model.

This documentation refers to PROFASI version 1.5, publicly available from 1 March
2011. The current version number is 1.4.9. The release version of PROFASI 1.5
will be distributed through an open access git repository. This is pre-release
beta code.

The documentation is automatically generated using the <a href="http://www.doxygen.org">doxygen</a> program. If doxygen is installed on your system, you can generate this documentation locally on your computer by running "make docs" in the main PROFASI directory. The documentation on <a href="http://cbbp.thep.lu.se/activities/profasi/"> the PROFASI homepage</a>
is by now out of date, and does not refer to this version.
</div>

<div style="text-align:justify">
The model implemented in PROFASI successfully describes the folding and thermodynamic behaviour of a number of peptides of both \f$ \alpha \f$-helical and \f$ \beta \f$-sheet secondary structure with about 20 residues [\ref IIa, \ref IIb, \ref III]. At least on one instance, it also describes the folding of a 49-residue protein with both \f$ \alpha \f$-helical and  \f$ \beta \f$-sheet secondary structure elements and a complex topology [\ref IV], using unbiased replica exchange Monte Carlo simulations starting from random initial conformations.

PROFASI has been used in a number of studies of amyloid aggregation, with up to 30 peptide chains in full atomic detail [\ref V, \ref VI, \ref VII]. Other interesting applications include  studies of mechanical and thermal unfolding of globular proteins [\ref VIII, \ref IX, \ref X] as well as studies of small semiconductor-binding peptides [\ref XI].

For a detailed description of the model and the interaction potential, please refer to [\ref IIa]. A gallery of some interesting simulations done with this program can be found <a href="http://www.fz-juelich.de/jsc/slbio/PROFASI_Gallery/index.html">here</a>. If any publication should result using this program package, please cite it through [\ref I] and the web address of the program: "http://cbbp.thep.lu.se/activities/profasi/".

PROFASI is freely available under a licence similar to the GNU General Public Licence, but is restricted to academic users.
</div>

\li We welcome your involvement in this project. If you use PROFASI and find
errors in the code, causing it to crash where it should not, give wrong answers,
hang etc., we would appreciate if you send a note to "profasi at thep.lu.se". Even
better, if you fix the problem yourself, it would be a nice thing to do to share
it with other people using PROFASI for their research. Please send a patch file
containing your fix/enhancement to the same email address. We will review the
patch and do the necessary corrections to ensure compatibility with other changes
in the code, and then make the patch available in the PROFASI updates page, with
due credit to you. A patch will be considered for inclusion if and only if it does
not break the build. \n

\li For a quick introduction to get started with PROFASI take a look at \ref
tutorial_1 . Whether or not you have any experience with the program, this
tutorial should give you a feel for it quickly.  \n

\li If you have used an older version of PROFASI for some time, you can get an
idea of the main changes in this version by reading \ref new_features_15 \n

\li The \ref toc page will help you navigate through the documentation.

\section references References
\anchor I \e [I] <b>PROFASI: A Monte Carlo simulation package for protein folding and aggregation</b>, A. Irb&auml;ck and S. Mohanty, <em>J. Comput. Chem.</em> <b>27</b>, 1548-1555 (2006)

\anchor IIa \e [IIa] <b>An effective all-atom potential for proteins</b>, A. Irb&auml;ck, S. Mitternacht and S. Mohanty, <em>PMC Biophysics</em> <b>2</b>, 2 (2009)

\anchor IIb \e [IIb] <b>Folding thermodynamics of peptides</b>, A. Irb&auml;ck and S. Mohanty, <em>Biophys. J.</em> <b>88</b>, 1560-1569 (2005)

\anchor III \e [III] <b>Folding of proteins with diverse folds</b>, S. Mohanty and U.H.E. Hansmann, <em>Biophys. J.</em> <b>91</b>, 3573-3578 (2006)

\anchor IV \e [IV] <b>Simulation of Top7-CFr: A transient helix extension guides folding</b>, S. Mohanty, J.H. Meinke, O. Zimmermann and U.H.E. Hansmann, <em>Proc. Natl. Acad. Sci. USA</em> <b>105</b>, 8004-8007 (2008)

\anchor V \e  [V] <b>Oligomerization of amyloid \f$ \mathbf{A\beta_{16-22}} \f$ peptides using hydrogen bonds and hydrophobicity forces</b>, G. Favrin, A. Irb&auml;ck and S. Mohanty, <em>Biophys. J.</em> <b>87</b>, 3657-3664 (2004)

\anchor VI \e [VI] <b>Structural reorganisation and potential toxicity of oligomeric species formed during the</b>, M. Cheon, I. Chang, S. Mohanty, L.M. Luheshi, C. M. Dobson, M. Vendruscolo and G. Favrin, <em>PLoS Comput. Biol.</em> <b>3</b>, e173 (2007)

\anchor VII \e [VII] <b>Spontaneous \f$ \beta \f$-barrel formation: An all-atom Monte Carlo study of \f$ \mathbf{A\beta_{16-22}} \f$ oligomerization</b>, A. Irb&auml;ck and S. Mitternacht, <em>Proteins</em> <b>71</b>, 207-214 (2008)

\anchor VIII \e [VIII] <b>Dissecting the mechanical unfolding of ubiquitin</b>, A. Irb&auml;ck, S. Mitternacht and S. Mohanty, <em>Proc. Natl. Acad. Sci. USA</em> <b>102</b>, 13427-13432 (2005)

\anchor IX \e [IX] <b>Thermal versus mechanical unfolding of ubiquitin</b>, A. Irb&auml;ck and S. Mitternacht, <em/>Proteins</em> <b>65</b>, 759-766 (2006)

\anchor X \e [X] <b>Changing the mechanical unfolding pathway of FnIII-10 by tuning the pulling strength</b>, S. Mitternacht, S. Luccioli, A. Torcini, A. Imparato and A. Irb&auml;ck, submitted (2008)

\anchor XI \e [XI] <b>Differences in solution behavior among four semiconductor-binding peptides</b>, S. Mitternacht, S. Schnabel, M. Bachmann, W. Janke and A. Irb&auml;ck, <em>J. Phys. Chem.</em> <b>B 111</b>, 4355-4360 (2007)

*/

