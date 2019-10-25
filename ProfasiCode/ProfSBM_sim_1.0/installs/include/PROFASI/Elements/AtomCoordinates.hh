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

#ifndef AtomCoordinates_HH
#define AtomCoordinates_HH
#include "../Aux/Vector3.hh"
#include <valarray>
#include <cstdio>

namespace prf
{
//! A class to store the coordinates of an atom.
    /**
     * \ingroup building_blocks
     * AtomCoordinates class is meant to work in two ways. On the one hand it is
     * the x,y and z coordinates of an atom. So, if you subtract two AtomCoordinate
     * objects you get a Vector3. But there are some advantages in not
     * using the Vector3 directly for the coordinates of an atom.
     * In order to impose the periodic boundary condition, one has to check
     * the extent of each molecule of the system along all three axes separately.
     * This is cumbersome if each component of the position vectors is only
     * accessible through the three vector. As far as periodicity is concerned,
     * the different spatial directions have nothing to do with eachother.
     *
     * Besides, in this implementation, the AtomCoordinates class will actually
     * store coordinates of all atoms in static valarrays. The unique_id of
     * each atom marks the element in the static arrays reserved for that atom.
     * While this is by no means crucial, it has been done for two reasons...
     * <ul>
     * <li> For certain frequent operations, such as translation or rotation of a
     * block of atoms, one is concerned with an operation in space. A certain set
     * of vectors in space are translated or rotated. It does not matter that
     * those vectors belong to atoms and so on. Conformational updates are often
     * such operations on space. So, the space of all coordinates is a natural
     * and important concept. Having a class storing all coordinates at one place
     * simplifies implementation of conformational updates.</li>
     * <li> An atom only needs to store one integer id for its coordinates. So,
     * any number of copies of the atom can be made without duplicating the spatial
     * coordinates. When an AtomCoordinates object is copied, for instance when it
     * is passed as an argument to a function, only the integer label is copied.
     * Only at the final step, when it is really necessary, is it converted to
     * a Vector3 object. Number of copying operations is reduced. </li>
     * </ul>
     * AtomCoordinates also saves a backup set of coordinates for each atom.
     * This helps quickly restore the system when an update is rejected.
     */

    class AtomCoordinates
    {

    public:
        explicit AtomCoordinates(int i);
        //!<create coordinates for atom with unique_id i
        ~AtomCoordinates();
        AtomCoordinates(const AtomCoordinates &);
        //! An important initialization function
        /**
         * The function Initialize(tot_atoms) allocates space for the
         * coordinates of tot_atoms Atoms. It should be called once, at
         * the beginning when the total number of atoms has been determined.
         */
        static void Initialize(int tot_atoms);
        //! Update backup coordinates
        /** Copy a block of coordinates "forward", so that backups are
         * synchronized with the current values */
        static void update(int strt, int iend);
        //! Revert to backup coordinates
        /** Copy a block of coordinates "backward", so that current values are
         * abandoned and replaced with backup values */
        static void revert(int strt, int iend);
        //! Center of mass of a block of coordinates, expensive!
        static Vector3 CenterOfMass(int strt, int iend);
        
        static Vector3 CenterOfMass3(int i1, int i2, int i3);
        //! Move center of mass to origin
        static void center();
        //! Radius of gyration of a block, again, use carefully!
        static double rgyr(int istrt, int iend);
        //! Rotate a block of atoms
        /**
         * Rotate a block of atoms between strt and iend by angle tht about the
         * axis passing through the atomCoordinates ax2 and ax1.
         */
        static void BlockRotate(double tht,int strt, int iend,
                                AtomCoordinates ax1, AtomCoordinates ax2);
        //! Rotate block about given axis and origin
        static void BlockRotate(double tht,int strt, int iend,
                                Vector3 org,
                                Vector3 axs);
        //! translate a block of coordinates by given amounts.
        static void BlockTranslate(const Vector3 & amnt,
                                   int strt, int iend);
        static void BoundingBox(int i1, int i2,
                                       double &xmn, double &xmx,
                                       double &ymn, double &ymx,
                                       double &zmn, double &zmx);
        //!< find bounding box for atoms between i1 and i2
        static void EnforceBC(int i1, int i2);
        //!< enforce periodic boundary conditions on atoms between i1 and i2
        static void write_pdb_box_info(FILE *fp);
        //!< Periodic box information to a pdb file
        static void CheckDiff();
        //!< check that the differences between current and backup coordinates
        //!< are not too big. Mostly a sanity check function
        //! Switch to make the back up coordinates the active ones
        /**
         * This does not dump the current coordinates, but causes energy
         * calculations or for that matter any external object to access the
         * backup values as the coordinates of different atoms.
         */
        static inline void switch_to_backups() {uc=&btc;}

        //! Switch back to current coordinates
        /**
         * This reverses what is done in switch_to_backups
         */
        static inline void switch_to_regular() {uc=&atc;}

        //! Set size of periodic box
        static inline void SetBox(double bx) {BOXL=bx;HALFBX=0.5*BOXL;}

        //! access box size
        static inline double boxL() {return BOXL;}

        static inline double halfBL() {return HALFBX;}

        static inline int numberOfAtoms() {return ntotatoms;}

        static double total_error();
        //!< difference between present coordinates and backups
        //@{
        /**
         * @name Accessing coordinate data
         * access or modify properties of an individual AtomCoordinates object
         */
        inline double x() const {return (*uc)[locn];}

        inline double y() const {return (*uc)[locn+1];}

        inline double z() const {return (*uc)[locn+2];}

        inline void x(double val) {(*uc)[locn]=val;}

        inline void y(double val) {(*uc)[locn+1]=val;}

        inline void z(double val) {(*uc)[locn+2]=val;}

        //! return the unique_id corresponding to the AtomCoordinates object
        inline int Id() const {return locn/3;}

        inline void Id(int i) {locn=3*i;}

        //@}
        //! Access coordinates without an AtomCoordinate object
        //@{
        /**
         * @name Static coordinate access without AtomCoordinates object
         * It is also possible to access the x,y, z coordinates of the i'th
         * atom without ever referring to the atom or an AtomCoordinates
         * object. This is done using static functions using the unique_id
         * of that atom as the argument.
         */
        static inline double xval(int i) {return (*uc)[3*i];}

        static inline double yval(int i) {return (*uc)[3*i+1];}

        static inline double zval(int i) {return (*uc)[3*i+2];}

        //! difference between current and backup coordinates for atom i
        static inline Vector3 diff_from_backup(int i) {
            i*=3;
            return Vector3(atc[i]-btc[i],atc[i+1]-btc[i+1],atc[i+2]-btc[i+2]);
        }

        //@}
        //@{
        /**
         * @name Three-vector like operations
         * One can manipulate AtomCoordinates just like three vectors.
         * Algebraic operations are supported among them which yield answers
         * which would be expected had they really been three vectors. Even
         * algebraic operations between an ordinary three vector and an
         * AtomCoordinates object should work seamlessly.
         */
        //! extract three vector out of the AtomCoordinates object
        inline Vector3 value() const { return Vector3((*uc)[locn],(*uc)[locn+1],(*uc)[locn+2]); }

        //! assign x,y and z coordinates
        inline void value(double gx, double gy, double gz)
        {(*uc)[locn]=gx;(*uc)[locn+1]=gy;(*uc)[locn+2]=gz;}

        //! statically get the Vector3 representing the coordinates of atom i
        static inline Vector3 vec(int i) {i*=3;return Vector3((*uc)[i],(*uc)[i+1],(*uc)[i+2]);}

        //! statically assign x,y and z values to atom i
        static inline void vec(int i, double gx, double gy, double gz)
        {i*=3;(*uc)[i]=gx;(*uc)[i+1]=gy;(*uc)[i+2]=gz;}

        //! Copy x,y and z coordinates from another AtomCoordinates object
        AtomCoordinates & operator= (const AtomCoordinates &);
        //! Assign a Vector3 to an AtomCoordinates
        inline AtomCoordinates & operator= (const Vector3 &);

        inline Vector3 operator+ (const AtomCoordinates &) const;
        inline Vector3 operator+ (const Vector3 &) const;
        inline AtomCoordinates & operator+= (const AtomCoordinates &);
        inline AtomCoordinates & operator+= (const Vector3 &);
        friend Vector3 operator+ (const Vector3 &,
                                  const AtomCoordinates &);
        inline Vector3 operator- (const AtomCoordinates &) const;
        inline Vector3 operator- (const Vector3 &)
        const;
        inline AtomCoordinates & operator-= (const AtomCoordinates &);
        inline AtomCoordinates & operator-= (const Vector3 &);
        friend Vector3  operator- (const Vector3 &,
                                   const AtomCoordinates &);
        //! Scalar product of two AtomCoordinates
        inline double dot(const AtomCoordinates &) const;
        //! Scalar product of an AtomCoordinates object with a Vector3
        inline double dot(const Vector3 &) const;
        //! Vector (cross) product of two AtomCoordinates
        inline Vector3 operator*(const AtomCoordinates &) const;
        //! Vector product of an AtomCoordinates object with a Vector3
        inline Vector3 operator*(const Vector3 &) const;
        //@}
        //@{
        /**
         * @name Distance measure related functions
         */
        //!Periodic separation between values x1 and x2
        static double per_sep(double x1, double x2);
        //!Positive definite periodic separation
        static double pos_per_sep(double x1, double x2);
        //! components of the periodic dist. between atoms indexed i1 and i2
        static inline double xdist(int i1, int i2)
        {return per_sep((*uc)[3*i1], (*uc)[3*i2]);}

        static inline double ydist(int i1, int i2)
        {return per_sep((*uc)[3*i1+1], (*uc)[3*i2+1]);}

        static inline double zdist(int i1, int i2)
        {return per_sep((*uc)[3*i1+2], (*uc)[3*i2+2]);}

        //!non-periodic versions of the separations along x,y and z axes
        static inline double npxdist(int i1, int i2)
        {return(*uc)[3*i1]- (*uc)[3*i2];}

        static inline double npydist(int i1, int i2)
        {return(*uc)[3*i1+1]- (*uc)[3*i2+1];}

        static inline double npzdist(int i1, int i2)
        {return(*uc)[3*i1+2]- (*uc)[3*i2+2];}

        //!periodic separation between atom i and a given point along an axis
        static inline double sepx(int i, double xpt)
        {return per_sep((*uc)[3*i],xpt);}

        static inline double sepy(int i, double ypt)
        {return per_sep((*uc)[3*i+1],ypt);}

        static inline double sepz(int i, double zpt)
        {return per_sep((*uc)[3*i+2],zpt);}

        static inline Vector3 sepv(int i, Vector3 v) {
            i*=3;
            return Vector3(per_sep((*uc)[i],v.x()),
                           per_sep((*uc)[i+1],v.y()),
                           per_sep((*uc)[i+2],v.z()));
        }

        //!positive definite periodic separation between atom i and a point
        static inline double psepx(int i, double xpt)
        {return pos_per_sep((*uc)[3*i],xpt);}

        static inline double psepy(int i, double ypt)
        {return pos_per_sep((*uc)[3*i+1],ypt);}

        static inline double psepz(int i, double zpt)
        {return pos_per_sep((*uc)[3*i+2],zpt);}


        //!bare cartessian distance measure between atoms i1 and i2
        static inline double s2(int i1,int i2) {
            i1*=3;
            i2*=3;
            return sqr((*uc)[i1] - (*uc)[i2])+
                   sqr((*uc)[i1+1] - (*uc)[i2+1])+
                   sqr((*uc)[i1+2] - (*uc)[i2+2]);
        }

        static inline double s(int i1,int i2) {return sqrt(s2(i1,i2));}

        //!distance measure respecting periodic box
        static inline double dist2(int i1, int i2);
        static inline double dist(int i1, int i2) {
            return sqrt(dist2(i1,i2));
        }

        //! Periodic distance measure from a given point
        static inline double seps2(int i,
                                   double xpt,
                                   double ypt,
                                   double zpt);
        //! Distance^2 from a point using positive definite periodic measures
        static inline double pseps2(int i,
                                    double xpt,
                                    double ypt,
                                    double zpt);
        static inline double seps(int i,
                                  double xpt,
                                  double ypt,
                                  double zpt)
        {return sqrt(seps2(i,xpt,ypt,zpt));}

        static inline double pseps(int i,
                                   double xpt,
                                   double ypt,
                                   double zpt)
        {return sqrt(pseps2(i,xpt,ypt,zpt));}

        //! Non-periodic vector difference between i1 and i2
        static inline Vector3 diff(int i1, int i2) {
            i1*=3;
            i2*=3;
            return Vector3((*uc)[i1]-(*uc)[i2],
                           (*uc)[i1+1]-(*uc)[i2+1],
                           (*uc)[i1+2]-(*uc)[i2+2]);
        }

        //! Periodic vector difference between i1 and i2
        static inline Vector3 sep(int i1, int i2) {
            i1*=3;
            i2*=3;
            return Vector3(per_sep((*uc)[i1],(*uc)[i2]),
                           per_sep((*uc)[i1+1],(*uc)[i2+1]),
                           per_sep((*uc)[i1+2],(*uc)[i2+2]));
        }

        //! Positive-definite periodic Vector difference between i1 and i2
        static inline Vector3 psep(int i1, int i2) {
            i1*=3;
            i2*=3;
            return Vector3(pos_per_sep((*uc)[i1],(*uc)[i2]),
                           pos_per_sep((*uc)[i1+1],(*uc)[i2+1]),
                           pos_per_sep((*uc)[i1+2],(*uc)[i2+2]));
        }

    private:
        static std::valarray<double> atc,btc;
        static double BOXL,HALFBX;
        static std::valarray<double> *uc;
        static int ntotatoms;
        static inline double sqr(double x) {return x*x;}

    private:
        int locn;
    };
}

#include "AtomCoordinates.i.cc"
#endif
