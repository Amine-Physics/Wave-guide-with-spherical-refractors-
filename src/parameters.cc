
#include "parameters.h" 
#include<deal.II/base/parameter_handler.h> 
#include<deal.II/base/exceptions.h>
#include<deal.II/base/utilities.h> 
#include<deal.II/base/point.h>
#include<vector>
#include<string> 
#include<cmath>
#include<iostream> 
#include<fstream> 




namespace helpers {
	using namespace dealii; 
	void declare(ParameterHandler &param) {
	
		// Mesh and Geometry	
		param.enter_subsection("Mesh & Geometry");
			param.declare_entry("Number Of Refinements", "2", Patterns::Integer(0), "Number of Refinement Steps"); 
			param.declare_entry("Wave Guide Corners", "-1.0 , 0.0, 1.0, 0.5", Patterns::List(Patterns::Double() , 4 , 4), 
					"x_min, y_min, x_max, y_max"); 
			param.declare_entry("Max Cell Size", "0.05", Patterns::Double(0), " Maxium cell size during mesh generation"); 
		param.leave_subsection(); 
		
		//Refractors
		
		param.enter_subsection("Refractors"); 
			param.declare_entry("Centers x", "0.5" , Patterns::List(Patterns::Double()), " a list of the x coordinates of the refractors centers"); 
			param.declare_entry("Centers y",  "0.25" , Patterns::List(Patterns::Double()), " a list of the y coordinates of the refractors centers"); 
			param.declare_entry("Radii" , "0.1", Patterns::List(Patterns::Double()), "a list of the radii of the refractors"); 
			param.declare_entry("Permitivities", "1.5",  Patterns::List(Patterns::Double()) , "a list of the relative permitivties of the refractors"); 
		param.leave_subsection();
		
		// physics 
		param.enter_subsection("Physics"); 
			param.declare_entry("Omega", "3e9", Patterns::Double(0), "waves angular velocity"); 
			param.declare_entry("C", "3e8", Patterns::Double(0), "Light's velocity in empty space"); 
			param.declare_entry("Mode Number", "1",  Patterns::Integer(0) , "Incoming wave mode number");
		param.leave_subsection(); 
		
		// output 
		param.enter_subsection("OutPut"); 
			param.declare_entry("OutPut File Name", "solution", Patterns::Anything(), "Output file name without entension"); 
		param.leave_subsection();
		
		// PostProcessing 
		
		param.enter_subsection("PostProcessing");
			param.declare_entry("Number Of Subdivisions","400",  Patterns::Integer(0) , "Number of transversal subdivisions  needed in defining vertical integration points"); 
			param.declare_entry("Mode", "2" , Patterns::Integer(0) , " the mode number for which we cant to compute the mode amplitude in the solution."); 					
			param.declare_entry("Computation Positions", "-0.95 , -0.25 , 0.25 ,  0.95", Patterns::List(Patterns::Double()) , "the list of the x values "
												"where we want to the modal amplitude and intensity");
			param.declare_entry("Number Of Points", "1000", Patterns::Integer(0) , "Number of longititunal points for graphs computation");											    
		param.leave_subsection(); 
		
			
	}
	
	std::vector<double> get_list (ParameterHandler &param, const std::string &key) {
		std::vector<double> list = Utilities::string_to_double(Utilities::split_string_list(param.get(key))) ; 
		return list; 
	
	}

} 

using namespace helpers ; 
Parameters read_parameters( const std::string & parameter_file) { 
	ParameterHandler param ;
	declare(param); 
	//std::ofstream out(parameter_file);
	//param.print_parameters(out, ParameterHandler::ShortPRM);
	param.parse_input(parameter_file); 
	
	Parameters p; 
	
	// filling the parameters sructure 
	
	param.enter_subsection("Mesh & Geometry"); 
	p.n_refinements = param.get_integer("Number Of Refinements"); 
	
	const auto corners = get_list(param, "Wave Guide Corners");
	p.x_min = corners[0]; 
	p.y_min = corners[1]; 
	p.x_max = corners[2]; 
	p.y_max = corners[3]; 
	
	p.h_max = param.get_double("Max Cell Size"); 
	p.W = p.y_max - p.y_min ;
	param.leave_subsection(); 
	
	param.enter_subsection("Physics"); 
	p.c = param.get_double("C"); 
	p.omega = param.get_double("Omega"); 
	p.n = param.get_integer("Mode Number"); 
	param.leave_subsection(); 
	
	param.enter_subsection("Refractors"); 
	
	const auto centers_x = get_list(param, "Centers x"); 
	const auto centers_y = get_list(param, "Centers y"); 
	const auto radii = get_list(param, "Radii");
	const auto eps_r = get_list(param , "Permitivities"); 
	param.leave_subsection(); 
	AssertThrow(centers_x.size()  == centers_y.size() && centers_x.size() == radii.size() && centers_x.size() == eps_r.size(), 
		ExcMessage("Refractor lists should be the same size")); 
		
	std::vector<Point<2>> centers(centers_x.size()) ; 
	for (unsigned int i = 0 ; i < centers.size() ; ++i) {
		centers[i] = Point<2>( centers_x[i] , centers_y[i]) ; 
	}
	 
	std::vector<std::pair<int,int>> overlaping_refractors;
	for (unsigned int i = 0 ; i < centers.size() ; i++){
		for (unsigned int j = 0 ; j < centers.size() ; j++){
			if (j > i && (centers[i] - centers[j]).norm() < radii[i] + radii[j]) {
				overlaping_refractors.emplace_back(i,j);
			}
		}
	}
	
	AssertThrow(overlaping_refractors.size() == 0 , ExcMessage("The refractors domain ovelape !")); 
	if(overlaping_refractors.size() != 0) {
		std::cout<< "The overlaping refractors are" << std::endl ; 
		for (unsigned int i = 0 ; i < overlaping_refractors.size() ; i++){
			std::cout<< i << "| " << overlaping_refractors[i].first << " | " << overlaping_refractors[i].second << std::endl;
		}
	}
	
	 for (unsigned int i = 0 ; i < centers.size();i++){
	 	p.refractors.push_back({centers_x[i], centers_y[i], radii[i], eps_r[i]}); 
         		
         }
	for (const auto &d : p.refractors){
    		AssertThrow(d.cx - d.radius > p.x_min && d.cx + d.radius < p.x_max &&
                  		  d.cy - d.radius > p.y_min && d.cy + d.radius < p.y_max,
                		   ExcMessage("A refractor does not fit inside the waveguide")); }
                
        
         param.enter_subsection("OutPut"); 
         p.output_file_name = param.get("OutPut File Name");
         param.leave_subsection();
         
         param.enter_subsection("PostProcessing"); 
         p.n_sub = param.get_integer("Number Of Subdivisions");
         p.m = param.get_integer("Mode"); 
         p.x_values = get_list(param , "Computation Positions");
         p.n_points = param.get_integer("Number Of Points");
         param.leave_subsection();
         return p;
	
}


























