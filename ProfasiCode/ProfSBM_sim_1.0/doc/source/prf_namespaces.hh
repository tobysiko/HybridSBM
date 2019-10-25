#ifndef PRF_NAMESPACES_HH
#define PRF_NAMESPACES_HH
//Contents of this file are required only for documentation

//! The namespace containing the key classes of PROFASI
/**
* The prf namespace contains classes representing the main concepts behind 
* PROFASI, for instance: the building blocks like Atom, AtomCoordinates, 
* AminoAcid, Protein, Population; the energy terms like Bias or 
* Hydrophobicity; the conformational updates like BGS or Rotation; the Monte 
* Carlo methods, like SimTemp (simulated tempering)... These classes together
* (roughly) define this model. In principle, energy terms are a different 
* kind of classes than updates, and they could be organized in separate 
* namespaces inside the prf namespace. This may be done in the future.   
*/
namespace prf {}
//! The namespace containing auxiliary useful classes used inside PROFASI
/**
 * Classes in this namespace are of more generic nature and usefulness. They
 * have little to do with the protein folding model implemented here, but 
 * they have been found as very useful tools. Two examples would be the 
 * Vector3, to represent vectors in 3-space, and the Histogram class for a 
 * variety of possible applications in connection with gathering statistics.
 */
namespace prf_utils {}
//! Application side classes expected to be used frequently
/**
 * This class presently contains only the ObsHandler and InterfaceBase. These
 * will, in all likelihood, appear frequently in various types of applications.
 * They do not define what PROFASI is, but rather constiture tools for 
 * conveniently using the model in an organized manner. 
 */
namespace prf_app {}
#endif
