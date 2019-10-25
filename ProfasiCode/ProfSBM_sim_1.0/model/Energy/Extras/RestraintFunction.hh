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

#ifndef RESTRAINTFUNCTION_HH
#define RESTRAINTFUNCTION_HH
#include "../../Aux/prf_xml.hh"
#include <deque>
#include <string>
#include <vector>

namespace prf {
    //! Simple harmonic potential well for arbitrary variable
    /**
      The base class restraint function implements the interface used by
      all other restraint types. But it also provides a useful restraint
      function itself: a harmonic oscillator potential (quadratic). There
      are two adjustable parameters: mean and overall weight.
      */
    class RestraintFunction {
    public:
        //! Default constructor
        RestraintFunction();
        virtual ~RestraintFunction();
        //! Set parameters form an XML node
        /**
          The node must be called \c parameters . It must have a field
          called \c mean and another called \c weight. These fields are
          inheritted by all other restraint functions, so they don't need
          their own.
          */
        
        virtual int set_pars(prf_xml::XML_Node *pars);
        //! Evaluate the function for a given value of the coordinate
        virtual double operator()(double x);
        //! Estimate a maximum based on some "large" scale
        virtual double estimate_max(double scale_large);
        //! Estimate a minimum
        /**
          The default value 0 is convenient for most cases. When the variable
          is at its mean location, the function should be 0, otherwise it
          should be some positive value.
          */
        virtual double estimate_min();
        virtual void setPerRestraintEnergy(double e);
        inline double getPerRestraintEnergy(){return perRestraintEnergy;}
    protected:
        double offset, weight;
        double perRestraintEnergy;
    };

    //! Generic power law potential
    /**
      Just like the quadratic or harmonic restraint term
      (\ref prf::RestraintFunction ), but with a general power law shape.
      There is one additional parameter called \c exponent .
      */
    class PowerLawRestraint : public RestraintFunction {
    public:
        PowerLawRestraint();
        ~PowerLawRestraint();
        //! Set parameters from an XML node
        /**
          Parameters are mean, weight and exponent.
          */

        int set_pars(prf_xml::XML_Node *pars);
        double operator()(double x);
    private:
        double exponent;
    };

    //! Power law with a small distance cut
    /**
      If the distance is greater than what would be the mean for the
      power law, it is equivalent to the power law. But if the atoms are
      closer, the potential does not grow but stays at 0.
      */
    class FlattenedPL : public RestraintFunction {
    public:
        FlattenedPL();
        ~FlattenedPL();
        int set_pars(prf_xml::XML_Node *pars);
        double operator()(double x);
    private:
        double exponent;
    };

    //! Restraint of the form of a Gaussian
    /**
      This is an inverted bell shaped curve centred at a given position.
      There is one parameter in addition to those in RestraintFunction, namely
      width.
      */
    class GaussianRestraint : public RestraintFunction {
    public:
        GaussianRestraint();
        ~GaussianRestraint();
        int set_pars(prf_xml::XML_Node *pars);
        double operator()(double x);
        double estimate_max(double scale_large);
    private:
        double width;
    };

    //! Circular normal potential well
    /**
      This is an adaptation of the circular normal or von Mises distribution
      for use as a potential well in angle space. The von Mises distribution
      is given by\n
      \f$
      f(\theta)= \frac{\exp(\kappa \cos(\theta-\mu))}{2\pi I_0 (\kappa)}
      \f$
      where \f$I_0\f$ is the modified Bessel function. The parameter \c kappa
      roughly corresponds to the inverse of the square of the width of the
      distribution. The angle in question does not appear in the denominator.
      It is therefore convenient to modify the function to :\n
       \f$
      f(\theta)=N \left( 1 - \frac{\exp(\kappa \cos(\theta-\mu))}{\exp(\kappa)} \right)
      \f$
      This way of writing the function makes sure that the value is 0 when the
      the angle takes its mean value, and N when it is not. The user can choose the
      location (mu), width ( roughly inverse square-root of kappa ) and the
      depth (N) of the potential well.
      */
    class CircularNormal : public RestraintFunction {
    public:
        CircularNormal();
        ~CircularNormal();
        //! Set parameters for circular normal potential
        /**
          The XML node passed should be a \c parameters node with up to 3
          child fields called \c kappa , \c mean and \c weight in
          arbitrary order. Example:
          \verbatim
          <parameters>
              <mean>1.1342</mean>
              <kappa>8</kappa>
              <weight>10</weight>
          <parameters>
          \endverbatim
          */
        int set_pars(prf_xml::XML_Node *pars);
        double operator()(double x);
        double estimate_max(double scale_large);
    private:
        double kappa;
    };
    
    //! Lennard-Jones native contact potential well
    /**
      Go-like native contact potential..
      */
    class LennardJonesNativeRestraint : public RestraintFunction { //NEW! -TS
    public:
    	LennardJonesNativeRestraint();
        ~LennardJonesNativeRestraint();
        //! Set parameters for LJ Go-like potential
        /**
          The XML node passed should be a \c parameters node with up to 2
          child fields called \c sigma ,  and \c epsilon in
          arbitrary order. Example:
          \verbatim
          <parameters>
              <sigma>1.1342</mean>
              <epsilon>8</kappa>
          <parameters>
          \endverbatim
          */
        int set_pars(prf_xml::XML_Node *pars);
        double operator()(double x);
        double estimate_min();
        double estimate_max(double scale_large);
        void setPerRestraintEnergy(double e);
    private:
        double minimum, epsilon;
    };
    
    class GaussianNativeRestraint : public RestraintFunction { //NEW! -TS
        public:
        	GaussianNativeRestraint();
            ~GaussianNativeRestraint();
            //! Set parameters for LJ Go-like potential
            /**
              The XML node passed should be a \c parameters node with up to 2
              child fields called \c sigma ,  and \c epsilon in
              arbitrary order. Example:
              \verbatim
              <parameters>
                  <sigma>1.1342</mean>
                  <epsilon>8</kappa>
              <parameters>
              \endverbatim
              */
            int set_pars(prf_xml::XML_Node *pars);
            double operator()(double x);
            double estimate_min();
            double estimate_max(double scale_large);
            void setPerRestraintEnergy(double e);
        private:
            double radius, steepness, minimum, width, depth;
    };
    
    class DualGaussianNativeRestraint : public RestraintFunction { //NEW! -TS
        public:
        	DualGaussianNativeRestraint();
        	~DualGaussianNativeRestraint();
            //! Set parameters for LJ Go-like potential
            /**
              The XML node passed should be a \c parameters node with up to 2
              child fields called \c sigma ,  and \c epsilon in
              arbitrary order. Example:
              \verbatim
              <parameters>
                  <sigma>1.1342</mean>
                  <epsilon>8</kappa>
              <parameters>
              \endverbatim
              */
            int set_pars(prf_xml::XML_Node *pars);
            double operator()(double x);
            double estimate_min();
            double estimate_max(double scale_large);
            void setPerRestraintEnergy(double e);
        private:
            double radius, steepness, minimum, width, depth, minimum2, width2, depth2;
    };

    class FixedDepthDualGaussianNativeRestraint : public RestraintFunction { //NEW! -TS
            public:
            	FixedDepthDualGaussianNativeRestraint();
            	~FixedDepthDualGaussianNativeRestraint();
                //! Set parameters for LJ Go-like potential
                /**
                  The XML node passed should be a \c parameters node with up to 2
                  child fields called \c sigma ,  and \c epsilon in
                  arbitrary order. Example:
                  \verbatim
                  <parameters>
                      <sigma>1.1342</mean>
                      <epsilon>8</kappa>
                  <parameters>
                  \endverbatim
                  */
                int set_pars(prf_xml::XML_Node *pars);
                double operator()(double x);
                double estimate_min();
                double estimate_max(double scale_large);
                void setPerRestraintEnergy(double e);
            private:
                double radius, steepness,minimum,width, depth,  minimum2, width2;
        };
    
    class FixedDepthMultiGaussianNativeRestraint : public RestraintFunction { //NEW! -TS
                public:
                	FixedDepthMultiGaussianNativeRestraint();
                	~FixedDepthMultiGaussianNativeRestraint();
                    //! Set parameters for LJ Go-like potential
                    /**
                      The XML node passed should be a \c parameters node with up to 2
                      child fields called \c sigma ,  and \c epsilon in
                      arbitrary order. Example:
                      \verbatim
                      <parameters>
                          <sigma>1.1342</mean>
                          <epsilon>8</kappa>
                      <parameters>
                      \endverbatim
                      */
                    int set_pars(prf_xml::XML_Node *pars);
                    double operator()(double x);
                    double estimate_min();
                    double estimate_max(double scale_large);
                    void setPerRestraintEnergy(double e);
                private:
                	std::vector<double> minima;
                    double radius, steepness, width, depth;
                    
                    std::vector<double> convertStringVectortoDoubleVector(const std::vector<std::string>& stringVector);
                    //double elemtodouble (const std::string& val);
                    std::vector<std::string> split(const std::string &s, char delim);
                    std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
                  };
    class FixedDepthMultiGaussianSmoothNativeRestraint : public RestraintFunction { //NEW! -TS
                public:
                	FixedDepthMultiGaussianSmoothNativeRestraint();
                	~FixedDepthMultiGaussianSmoothNativeRestraint();
                    //! Set parameters for LJ Go-like potential
                    /**
                      The XML node passed should be a \c parameters node with up to 2
                      child fields called \c sigma ,  and \c epsilon in
                      arbitrary order. Example:
                      \verbatim
                      <parameters>
                          <sigma>1.1342</mean>
                          <epsilon>8</kappa>
                      <parameters>
                      \endverbatim
                      */
                    int set_pars(prf_xml::XML_Node *pars);
                    double operator()(double x);
                    double estimate_min();
                    double estimate_max(double scale_large);
                    void setPerRestraintEnergy(double e);
                private:
                	
                    double low,high,radius, steepness, width, depth, lowcut, highcut,Rij,repulsion,val,wells,Gij;
                    
                  };

}

#endif // RESTRAINTFUNCTION_HH
