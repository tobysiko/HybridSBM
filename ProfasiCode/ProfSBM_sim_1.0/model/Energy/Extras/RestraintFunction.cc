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

#include "RestraintFunction.hh"
#include "../../Aux/Constants.hh"
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>

#include <algorithm>

#include <stdio.h>
//#include <transform>

using namespace prf;


RestraintFunction::RestraintFunction() : offset(0), weight(1) {
	perRestraintEnergy = 0;
}

RestraintFunction::~RestraintFunction() {}

int RestraintFunction::set_pars(prf_xml::XML_Node *pars)
{
    if (pars==NULL or pars->name()!="parameters") return 0;
    if (pars->child("mean")!=NULL)
        offset=strtod(pars->child("mean")->value().c_str(),NULL);
    else if (pars->child("offset")!=NULL)
        offset=strtod(pars->child("offset")->value().c_str(),NULL);
    if (pars->child("weight")!=NULL)
        weight=strtod(pars->child("weight")->value().c_str(),NULL);
    else if (pars->child("scale")!=NULL)
        weight=strtod(pars->child("scale")->value().c_str(),NULL);

    return 1;
}

double RestraintFunction::operator()(double x)
{
    return weight*(x-offset)*(x-offset);
}

double RestraintFunction::estimate_max(double scale_large)
{
    double offbkp=offset;
    offset=0;
    double lrg=operator()(scale_large);
    offset=offbkp;
    return lrg;
}

double RestraintFunction::estimate_min()
{
    return 0;
}

void RestraintFunction::setPerRestraintEnergy(double e){
	perRestraintEnergy = 0;
}


PowerLawRestraint::PowerLawRestraint() : RestraintFunction(), exponent(2) {}

PowerLawRestraint::~PowerLawRestraint() {}

int PowerLawRestraint::set_pars(prf_xml::XML_Node *pars)
{
    if (RestraintFunction::set_pars(pars)==0) return 0;
    if (pars->child("exponent")!=NULL)
        exponent=strtod(pars->child("exponent")->value().c_str(),NULL);
    return 1;
}


double PowerLawRestraint::operator()(double x)
{
    return weight*pow(fabs(x-offset),exponent);
}

FlattenedPL::FlattenedPL() : RestraintFunction(), exponent(2) {}

FlattenedPL::~FlattenedPL() {}

int FlattenedPL::set_pars(prf_xml::XML_Node *pars)
{
    if (RestraintFunction::set_pars(pars)==0) return 0;
    if (pars->child("exponent")!=NULL)
        exponent=strtod(pars->child("exponent")->value().c_str(),NULL);
    return 1;
}

double FlattenedPL::operator()(double x)
{
    return (x>offset)?weight*pow(x-offset,exponent):0;
}

GaussianRestraint::GaussianRestraint() : RestraintFunction(), width(1) {}

GaussianRestraint::~GaussianRestraint() {}

int GaussianRestraint::set_pars(prf_xml::XML_Node *pars)
{
    if (RestraintFunction::set_pars(pars)==0) return 0;
    if (pars->child("width")!=NULL)
        width=strtod(pars->child("width")->value().c_str(),NULL);
    return 1;
}

double GaussianRestraint::operator()(double x)
{
    return weight*(1-exp(-0.5*(x-offset)*(x-offset)/width/width));
}

double GaussianRestraint::estimate_max(double scale_large)
{
    return weight;
}

CircularNormal::CircularNormal() : RestraintFunction(), kappa(8) {}

CircularNormal::~CircularNormal() {}

int CircularNormal::set_pars(prf_xml::XML_Node *pars)
{
    if (RestraintFunction::set_pars(pars)==0) return 0;
    if (pars->child("kappa")!=NULL)
        kappa=strtod(pars->child("kappa")->value().c_str(),NULL);
    return 1;
}

double CircularNormal::operator ()(double x)
{
    return weight*(1-exp(kappa*(cos(x-offset)-1)));
}

double CircularNormal::estimate_max(double scale_large)
{
    return weight*(1-exp(-2*kappa));
}






LennardJonesNativeRestraint::LennardJonesNativeRestraint() : RestraintFunction(), minimum(1), epsilon(1) {
	if (perRestraintEnergy != 0) epsilon=perRestraintEnergy;
}  //NEW! -TS

LennardJonesNativeRestraint::~LennardJonesNativeRestraint() {}

int LennardJonesNativeRestraint::set_pars(prf_xml::XML_Node *pars)
{
    if (RestraintFunction::set_pars(pars)==0) return 0;
    if (pars->child("minimum")!=NULL)
    	minimum=strtod(pars->child("minimum")->value().c_str(),NULL);
    if (pars->child("epsilon")!=NULL)
            epsilon=strtod(pars->child("epsilon")->value().c_str(),NULL);
    return 1;
}

double LennardJonesNativeRestraint::operator()(double x)
{
	double ratio = minimum / x;
	double value = epsilon * ( pow(ratio, 12.0) - 2 * pow(ratio, 6.0 ) );
	//prf::cout <<"x:"<<x<<" min: "<<minimum<<" LJ:"<< value<<" \n";
	
    return value;
}

double LennardJonesNativeRestraint::estimate_max(double scale_large)
{
	prf::cout << "LJ estimate max: large scale = " << scale_large << "\n";
    return 0;
}
double LennardJonesNativeRestraint::estimate_min()
{
	prf::cout << "LJ estimate min \n";
    return epsilon * -1.0;
}

void LennardJonesNativeRestraint::setPerRestraintEnergy(double e){
	setPerRestraintEnergy(e);
	epsilon = e;
}



GaussianNativeRestraint::GaussianNativeRestraint() : RestraintFunction(), radius(1), steepness(1), minimum(1), width(1), depth(1) {
	if (perRestraintEnergy != 0) depth=perRestraintEnergy;
}  //NEW! -TS

GaussianNativeRestraint::~GaussianNativeRestraint() {}

int GaussianNativeRestraint::set_pars(prf_xml::XML_Node *pars)
{
    if (RestraintFunction::set_pars(pars)==0) return 0;
    
    if (pars->child("radius")!=NULL)
            radius=strtod(pars->child("radius")->value().c_str(),NULL);
    
    if (pars->child("steepness")!=NULL)
            steepness=strtod(pars->child("steepness")->value().c_str(),NULL);
        
    
    if (pars->child("minimum")!=NULL)
    	minimum=strtod(pars->child("minimum")->value().c_str(),NULL);
    
    if (pars->child("width")!=NULL)
        width=strtod(pars->child("width")->value().c_str(),NULL);
    
    if (pars->child("depth")!=NULL)
    	depth=strtod(pars->child("depth")->value().c_str(),NULL);
    
    return 1;
}

double GaussianNativeRestraint::operator()(double x)
{
	// with steepness=1 the function is 1 when x=radius
	double Rij;
	if (radius == 0.0)
		Rij = steepness * pow(radius/x, 12.0); // LJ repulsion
	else
		Rij = 0.0;
	double Gij = - exp( - pow(x-minimum, 2.0) / (2 * width * width) );
    
	//double result = Rij + depth * Gij + Rij * Gij;
	//prf::cout << "E:" << result << " x:" << x << " min:" << minimum << " depth:" << depth << " width:" << width << "\n";
	//prf::cout << "PerRestraintEnergy:" << perRestraintEnergy << "; depth:" << depth << "\n";
	return Rij + depth * Gij + Rij * Gij;//result;
}

double GaussianNativeRestraint::estimate_max(double scale_large)
{
	prf::cout << "LJ estimate max: large scale = " << scale_large << "\n";
    return 0;
}
double GaussianNativeRestraint::estimate_min()
{
	prf::cout << "LJ estimate min \n";
    return depth * -1.0;
}

void GaussianNativeRestraint::setPerRestraintEnergy(double e){
	perRestraintEnergy = e;
	depth = e;
}



DualGaussianNativeRestraint::DualGaussianNativeRestraint() : RestraintFunction(),  radius(1), steepness(1), minimum(1), width(1), depth(1), minimum2(1), width2(1), depth2(1) {}  //NEW! -TS

DualGaussianNativeRestraint::~DualGaussianNativeRestraint() {}

int DualGaussianNativeRestraint::set_pars(prf_xml::XML_Node *pars)
{
    if (RestraintFunction::set_pars(pars)==0) return 0;

    if (pars->child("radius")!=NULL)
            radius=strtod(pars->child("radius")->value().c_str(),NULL);

    if (pars->child("steepness")!=NULL)
            steepness=strtod(pars->child("steepness")->value().c_str(),NULL);


    if (pars->child("minimum")!=NULL)
    	minimum=strtod(pars->child("minimum")->value().c_str(),NULL);

    if (pars->child("width")!=NULL)
        width=strtod(pars->child("width")->value().c_str(),NULL);



    if (pars->child("depth")!=NULL)
    	depth=strtod(pars->child("depth")->value().c_str(),NULL);

    if (pars->child("minimum2")!=NULL)
    	minimum2=strtod(pars->child("minimum2")->value().c_str(),NULL);

    if (pars->child("width2")!=NULL)
            width2=strtod(pars->child("width2")->value().c_str(),NULL);

    if (pars->child("depth2")!=NULL)
        depth2=strtod(pars->child("depth2")->value().c_str(),NULL);
    return 1;
}

double DualGaussianNativeRestraint::operator()(double x)
{
	// with steepness=1 the function is 1 when x=radius
	double Rij;
	if (radius == 0.0)
		Rij = 0.0;

	else
		Rij = steepness * pow(radius/x, 12.0); // LJ repulsion


	double Gij;
	if (depth == 0.0 || width == 0.0)
		Gij = 0.0;
	else
		Gij = - exp( - pow(x-minimum, 2.0) / (2 * width * width) );

	double Gij2;
	if (depth2 == 0.0 || width2 == 0.0)
		Gij2 = 0.0;
	else
		Gij2 = - exp( - pow(x-minimum2, 2.0) / (2 * width2 * width2) );

	//double result = Rij + depth * Gij + Rij * Gij    + depth2 * Gij2 + Rij * Gij2;
	//prf::cout << result <<"\n";
	return Rij + depth * Gij + Rij * Gij    + depth2 * Gij2 + Rij * Gij2;//result ;
}
double DualGaussianNativeRestraint::estimate_max(double scale_large)
{
	prf::cout << "LJ estimate max: large scale = " << scale_large << "\n";
    return 0;
}
double DualGaussianNativeRestraint::estimate_min()
{
	prf::cout << "LJ estimate min \n";
    return (depth+depth2) * -1.0;
}
void DualGaussianNativeRestraint::setPerRestraintEnergy(double e){
	perRestraintEnergy = e;

}




FixedDepthDualGaussianNativeRestraint::FixedDepthDualGaussianNativeRestraint() : RestraintFunction(),  radius(1), steepness(1), minimum(1), width(1), depth(1), minimum2(1), width2(1) {
	if (perRestraintEnergy != 0) depth=perRestraintEnergy;
}  //NEW! -TS

FixedDepthDualGaussianNativeRestraint::~FixedDepthDualGaussianNativeRestraint() {}

int FixedDepthDualGaussianNativeRestraint::set_pars(prf_xml::XML_Node *pars)
{
    if (RestraintFunction::set_pars(pars)==0) return 0;

    if (pars->child("radius")!=NULL)
            radius=strtod(pars->child("radius")->value().c_str(),NULL);

    if (pars->child("steepness")!=NULL)
            steepness=strtod(pars->child("steepness")->value().c_str(),NULL);


    if (pars->child("minimum")!=NULL)
    	minimum=strtod(pars->child("minimum")->value().c_str(),NULL);

    if (pars->child("width")!=NULL)
        width=strtod(pars->child("width")->value().c_str(),NULL);



    if (pars->child("depth")!=NULL){
    	depth=strtod(pars->child("depth")->value().c_str(),NULL);
    	setPerRestraintEnergy(depth);
    }

    if (pars->child("minimum2")!=NULL)
    	minimum2=strtod(pars->child("minimum2")->value().c_str(),NULL);

    if (pars->child("width2")!=NULL)
            width2=strtod(pars->child("width2")->value().c_str(),NULL);

    return 1;
}

double FixedDepthDualGaussianNativeRestraint::operator()(double x)
{
	// with steepness=1 the function is 1 when x=radius
	double Rij;
	if (radius == 0.0)
		Rij = 0.0;

	else
		Rij = steepness * pow(radius/x, 12.0); // LJ repulsion

	double Gij;
	double Gij2;
	if (depth==0.0)
		return 0.0;
	else
		if (width == 0.0)
			Gij = 0.0;
		else
			Gij = - exp( - pow(x-minimum, 2.0) / (2 * width * width) );

		if (width2 == 0.0)
			Gij2 = 0.0;
		else
			Gij2 = - exp( - pow(x-minimum2, 2.0) / (2 * width2 * width2) );

	return depth * ( (1 + (1.0/depth)*Rij) * (1 + Gij) * (1 + Gij2) -1) ;
}
double FixedDepthDualGaussianNativeRestraint::estimate_max(double scale_large)
{
	prf::cout << "LJ estimate max: large scale = " << scale_large << "\n";
    return 0;
}
double FixedDepthDualGaussianNativeRestraint::estimate_min()
{
	prf::cout << "LJ estimate min \n";
    return depth * -1.0;
}
void FixedDepthDualGaussianNativeRestraint::setPerRestraintEnergy(double e){
	perRestraintEnergy = e;
	depth = e;
}

FixedDepthMultiGaussianNativeRestraint::FixedDepthMultiGaussianNativeRestraint() : RestraintFunction(), minima(1), radius(1), steepness(1), width(1), depth(1) {
	if (perRestraintEnergy != 0) depth=perRestraintEnergy;
}  //NEW! -TS

FixedDepthMultiGaussianNativeRestraint::~FixedDepthMultiGaussianNativeRestraint() {}

int FixedDepthMultiGaussianNativeRestraint::set_pars(prf_xml::XML_Node *pars)
{
    if (RestraintFunction::set_pars(pars)==0) return 0;

    if (pars->child("minima")!=NULL){

    	std::vector<std::string> minimaStrings = split(pars->child("minima")->value().c_str(), ',');
    	minima = convertStringVectortoDoubleVector(minimaStrings);
    	//prf::cout << "Contact minima found:\n";
    	//for (int i=0; i<minima.size(); i++ ){
    	//	prf::cout << minima[i]<<" , ";
    	//}
    	//prf::cout << "\n";
    }

    if (pars->child("radius")!=NULL){
    	radius=strtod(pars->child("radius")->value().c_str(),NULL);
    	//prf::cout << "Rep radius: "<< radius<< "\n";
    }

    if (pars->child("steepness")!=NULL){
            steepness=strtod(pars->child("steepness")->value().c_str(),NULL);
            //prf::cout << "Rep steepness: "<< steepness << "\n";
    }
    if (pars->child("width")!=NULL){
        width=strtod(pars->child("width")->value().c_str(),NULL);
        //prf::cout << "Width: "<< width << "\n";
    }

    if (pars->child("depth")!=NULL){
    	depth=strtod(pars->child("depth")->value().c_str(),NULL);
    	setPerRestraintEnergy(depth);
    	//prf::cout << "Depth: "<< depth << "\n";
    }


    return 1;
}

double FixedDepthMultiGaussianNativeRestraint::operator()(double x)
{
	// with steepness=1 the function is 1 when x=radius
	double Rij;
	double repulsion;
	
	if (radius==0.0){
		Rij=0.0;
		repulsion=1.0;
	}else{
		Rij = steepness * pow(radius/x, 12.0); 
		repulsion = (1.0 + (Rij/depth));
	}
	
	if (width==0 && depth==0){
		//prf::cout << "Non-native contact: x = " << x << ", radius = " << radius << ", Rij = " << Rij << "\n";
		return Rij;
	}
	
	double wells = 1.0;
	
	double Gij;
	for (int i=0; i< (int)minima.size(); i++ ){
		Gij = - exp( - pow(x-minima[i], 2.0) / (2.0 * width * width) );
		wells = wells * (1.0 + Gij);
	}
	
	double val =  (depth *  repulsion * wells ) - depth;
	
	if (val < -depth)
		prf::cout << "ERROR: FixedDepthMultiGaussianNativeRestraint invalid: " << val << "\n";
	return val;
}
double FixedDepthMultiGaussianNativeRestraint::estimate_max(double scale_large)
{
	prf::cout << "LJ estimate max: large scale = " << scale_large << "\n";
    return 0;
}
double FixedDepthMultiGaussianNativeRestraint::estimate_min()
{
	//prf::cout << "LJ estimate min \n";
    return depth * -1.0;
}
void FixedDepthMultiGaussianNativeRestraint::setPerRestraintEnergy(double e){
	perRestraintEnergy = e;
	depth = e;
}


std::vector<std::string> &FixedDepthMultiGaussianNativeRestraint::split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> FixedDepthMultiGaussianNativeRestraint::split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

double elemtodouble ( const std::string& val){
	//return std::stod(val);
	double t;
	std::stringstream convert( val  );
	convert>>t;
	return t;
	}

std::vector<double> FixedDepthMultiGaussianNativeRestraint::convertStringVectortoDoubleVector(const std::vector<std::string>& stringVector){
	std::vector<double> doubleVector(stringVector.size());
	
	transform(stringVector.begin(), stringVector.end(), doubleVector.begin(), elemtodouble );
	return doubleVector;
}





FixedDepthMultiGaussianSmoothNativeRestraint::FixedDepthMultiGaussianSmoothNativeRestraint() : RestraintFunction(), low(0), high(0), radius(0), steepness(0), width(0), depth(0), lowcut(0), highcut(0), wells(0), Rij(0), repulsion(0), val(0), Gij(0) {
	if (perRestraintEnergy != 0) depth=perRestraintEnergy;

}  //NEW! -TS

FixedDepthMultiGaussianSmoothNativeRestraint::~FixedDepthMultiGaussianSmoothNativeRestraint() {}

int FixedDepthMultiGaussianSmoothNativeRestraint::set_pars(prf_xml::XML_Node *pars)
{
    //if (RestraintFunction::set_pars(pars)==0) return 0;
    
    if (pars->child("low")!=NULL){
    	//low distance end of well
    	low = strtod(pars->child("low")->value().c_str(),NULL);
    }
    if (pars->child("high")!=NULL){
    	//high distance end of well
    	high = strtod(pars->child("high")->value().c_str(),NULL);
    }
    
    if (pars->child("radius")!=NULL){
    	radius=strtod(pars->child("radius")->value().c_str(),NULL);
    	//prf::cout << "Rep radius: "<< radius<< "\n";
    }
    
    if (pars->child("steepness")!=NULL){
            steepness=strtod(pars->child("steepness")->value().c_str(),NULL);
            //prf::cout << "Rep steepness: "<< steepness << "\n";
    }
    if (pars->child("width")!=NULL){
        width=strtod(pars->child("width")->value().c_str(),NULL);
        //prf::cout << "Width: "<< width << "\n";
        lowcut = low - (5.0 * width);
        highcut = high + (5.0 * width);
        //prf::cout << "Width: "<< width << "lowcut=" << lowcut << "; highcut=" << highcut <<"\n";
    }
    
    if (pars->child("depth")!=NULL){
    	depth=strtod(pars->child("depth")->value().c_str(),NULL);
    	setPerRestraintEnergy(depth);
    	//prf::cout << "Depth: "<< depth << "\n";
    }


    return 1;
}

double FixedDepthMultiGaussianSmoothNativeRestraint::operator()(double x)
{
	if (x>=highcut || x<=lowcut){
		val = -0.0;
	} else if (x<=low){
		val =  ( depth * (1.0 - exp( - pow(x-low , 2.0) / (2.0 * width * width) ) ) ) - depth;
	} else if (x >= high){
		val =  ( depth * (1.0 - exp( - pow(x-high, 2.0) / (2.0 * width * width) ) ) ) - depth;
	} else
		val = -depth;
	
	if (val < -depth)
		prf::cout << "ERROR: FixedDepthMultiGaussianSmoothNativeRestraint invalid: " << val << "\n";

	//prf::cout << "FixedDepthMultiGaussianSmoothNativeRestraint: " << val << " " << depth << " "<<x<< " "<<low<<" " <<high<<" "<<lowcut<<" "<<highcut<<"\n";
	//for (int i=0; i<minima.size(); i++ ){
	//	prf::cout << x-minima[i]<<",";
	//}
	//prf::cout << "\n";
	return val ;
}
double FixedDepthMultiGaussianSmoothNativeRestraint::estimate_max(double scale_large)
{
	//prf::cout << "FixedDepthMultiGaussianSmoothNativeRestraint estimate max: large scale = " << scale_large << "\n";
    return 0;
}
double FixedDepthMultiGaussianSmoothNativeRestraint::estimate_min()
{
	//prf::cout << "FixedDepthMultiGaussianSmoothNativeRestraint estimate min \n";
    return depth * -1.0;
}
void FixedDepthMultiGaussianSmoothNativeRestraint::setPerRestraintEnergy(double e){
	perRestraintEnergy = e;
	depth = e;
}
