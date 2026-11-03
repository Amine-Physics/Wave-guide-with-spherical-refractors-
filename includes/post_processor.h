/*---------------------------------------------------------------------------------------------------
This the post processor class declaration. This class takes a solution vector, a dof_handler, and a
parameter type object and provide information about the the amplitude of the refection and transmission 
coefficients, the amplitude of the wave at a given cross section caracterised with an x value. 
---------------------------------------------------------------------------------------------------*/

# pragma once 
#include<deal.II/lac/vector.h> 
#include<deal.II/dofs/dof_handler.h> 
#include<deal.II/numerics/fe_field_function.h>
#include"parameters.h"
#include<vector> 
#include<complex>


using namespace dealii; 

// declariing ScatteringResults structure 


struct ScatteringResults{
	std::complex<double> r,t ; //complex valued reflection and transmission coefficients 
	double R,T , balance ; // real values scattering and transmission coefficients and energy balance i.e R = |r|² , T = |t|² and balance = R + T

} ; 

struct AmplitudeSample{
	double position ; // x
	std::complex<double> value; // complex model amplitude a_m(x)

};

struct IntensitySample{
	double position ; // x 
	double value; // transverse energy 
};

struct Integration {
	std::vector<double> mapped_quad_points, mapped_weights; 

};

struct ModeCoefficients {
	unsigned int mode ; 
	double beta ; 
	std::complex<double>  A_left ; 
	std::complex<double>  B_left ;   
	std::complex<double>  A_right ; 
	std::complex<double>  B_right ;   
};

template<int dim> 
class PostProcessor {
	public: 
	        PostProcessor(const dealii::DoFHandler<dim> &dofh , const dealii::Vector<double>  &sol , const Parameters &prm, const Integration &inte) ;  // class constructor 

	        std::vector<AmplitudeSample> modal_amplitude_graph() const; 
	        std::vector<AmplitudeSample> modal_amplitude_point_list() const; 
	        std::vector<IntensitySample> intensity_graph() const ; 
	        std::vector<IntensitySample> intensity_point_list() const ;
	       
	        std::vector<ModeCoefficients> power_spectrum() const ;  
	        
	        ScatteringResults scattering_results()const ; 
	        
	private:
	         const dealii::DoFHandler<dim> &dof_handler; 
	         const dealii::Vector<double> &solution; 
	         const Parameters  &param;
	         const Integration integration_params;
	         Functions::FEFieldFunction < dim , DoFHandler<dim> , Vector<double>> field; 
	         
	         std::complex<double>  modal_amplitude(const double x_value, const unsigned int m) const ; // helper function that takes a an x value in the wave guide and comptes the amplitude of the m mode in the transver wave 
	         double intensity(double x_value) const; // helper function that taks an x value within the wave guide and compute the L 2 norm of the transverse wave
	         
	         std::vector<AmplitudeSample> modal_amplitude_list(const std::vector<double> &x_values)const ; // helper function that takes a list of points and computes the modal amplitude at each entry 
	         std::vector<IntensitySample> intensity_list(const std::vector<double> &x_values)const ; // helper function that takes a list of points and computes the intensity at each entry
	         
	         ModeCoefficients transmission_reflection( const unsigned int mode) const ; //  helper function that computes the reflection and transmission coefficients for a given mode 
	         
};

void get_integration_points_and_weights(const Parameters &p , Integration &inte);

void write_intensity_csv(const std::vector<IntensitySample> &intensitues , const std::string &output_file_name); 
void write_amplitude_csv(const Parameters &p, const std::vector<AmplitudeSample> &amplitudes , const std::string &output_file_name); 
void write_scattering_results_csv(const ScatteringResults &res , const std::string &output_file_name); 
































