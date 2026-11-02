
/*--------------------------------------------------------------------------------------------
This the header for the parameter declaration function. It produces a .prm file that containts 
the simulation parameters. This document is read at run time and hence eliminate the need for
recompilation. 

The function containts two structures: Refractors and Parameters. 

Refractors containts the information about the position and radiu of each refractor as well as 
there permitivity. 

Parameters containts the rest of the simulaion prameters: Mesh and Geometry, Incoming wave, 
number of refractors ... 
 

--------------------------------------------------------------------------------------------*/ 

#pragma once 
#include<vector>
#include<string>

struct Refractor {
	double cx , cy, radius, eps ;
}; 

struct Parameters {

	// Mesh & Geometry 
	double x_min , y_min, x_max , y_max; // dimentions of the wave guide
	double h_max;  		      // maximum size of cell during generation
	double W ; 
	unsigned int n_refinements; 	      // number of global refinements
	std::vector<Refractor> refractors ; // list containts the information about the refractors
	
	// Physics 
	
	double c , omega; // speed of ligh and wave angular velocity 
	unsigned int n ;  // mode number
	
	// output 
	
	std::string output_file_name; 
	
	// post proscesing 
	
	unsigned int n_sub;
	unsigned int n_points;
	unsigned int m; 
	std::vector<double> x_values; 
};

Parameters read_parameters(const std::string &parameter_file); 

