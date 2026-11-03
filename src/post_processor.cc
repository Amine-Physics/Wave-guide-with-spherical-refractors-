/*-----------------------------------------------------------------------------------------------------
This is the amplimentation of the post processor class 
------------------------------------------------------------------------------------------------------*/

#include<deal.II/numerics/vector_tools.h>
#include<deal.II/base/quadrature_lib.h>  
#include<deal.II/base/numbers.h>
#include<cmath> 
#include<fstream> 
#include<iostream> 
#include"post_processor.h"

using namespace dealii;

/*-- implementing the helper function get_integration_points_and_weights */ 
void get_integration_points_and_weights(const Parameters &p , Integration &inte) {
	
	const double W = p.y_max - p.y_min ; // wave guide width 
	const double h = W / p.n_sub ; // size of subdivision 
	
	QGauss<1> gauss(4); // gauss quad rule  
	
	std::vector<double> quad_points , weights ; 
	
	// loop over subdivision 
	for(unsigned int s = 0 ; s < p.n_sub ; ++s){ 
		// loop over quadrature points within a subdivision 
		for(unsigned int q = 0 ; q < gauss.size() ; ++q){
			double y = p.y_min + (s + gauss.point(q)[0]) * h ; 
			quad_points.push_back(y) ; 
			weights.push_back(gauss.weight(q) * h); 
		}// loop over quadrature points within a subdivision
		
	
	}// loop over subdivision 
	
	inte.mapped_quad_points = quad_points; 
	inte.mapped_weights = weights;

} // implementing the helper function get_integration_points_and_weights 




/*-- Implementing the class PostProcessor --*/ 

// constructor 
template<int dim> 
PostProcessor<dim>::PostProcessor(const DoFHandler<dim> &dofh , const Vector<double> &sol , const Parameters &prm, const Integration &inte): dof_handler(dofh) , 
													    solution(sol), 
													    param(prm),
													    integration_params(inte),
													    field(dofh , sol){} 

/*-- Implementing Helper Functions --*/

// modal_amplitude
template<int dim> 

std::complex<double> PostProcessor<dim>::modal_amplitude(const double x_value, const unsigned int m) const {
	
	AssertThrow( x_value< param.x_max && x_value > param.x_min , ExcMessage(" The x coordinate x sould lay stickly inside the wave guide compuation domaine")); 
	
	std::vector<Point<dim>> points; 
	 
	 for(unsigned int k = 0 ; k < integration_params.mapped_quad_points.size() ; ++k){
	 	points.emplace_back(x_value , integration_params.mapped_quad_points[k]);
	 }
	 
	 std::vector<Vector<double>> values(points.size(), Vector<double>(2)); 
	 field.vector_value_list(points, values); 
	 
	 std::complex<double> sum = 0.0;
	 for(unsigned int k = 0 ; k < points.size() ; ++k){
	 	double phi = std::sqrt(2.0/param.W) * std::sin((points[k](1) - param.y_min) / param.W * numbers::PI * m) ; 
	 	sum += std::complex<double>(values[k](0) , values[k](1)) * phi * integration_params.mapped_weights[k] ; 
	 }
	return sum ;
	 
} 

// modal_amplitude_list 
template<int dim> 
std::vector<AmplitudeSample> PostProcessor<dim>::modal_amplitude_list(const std::vector<double> &x_values) const {
	
	std::vector<AmplitudeSample> output; 
	
	for(unsigned int k = 0; k < x_values.size() ; ++k){
		output.push_back({x_values[k] , PostProcessor<dim>::modal_amplitude(x_values[k], param.m)});
	}
	return output;
	
}// modal_amplitude_list 

// intensity 

template<int dim> 
double PostProcessor<dim>::intensity(const double x_value) const {
	AssertThrow( x_value< param.x_max && x_value > param.x_min , ExcMessage(" The x coordinate x sould lay stickly inside the wave guide compuation domaine"));
	
	std::vector<Point<dim>> points;
	for(unsigned int k = 0 ; k < integration_params.mapped_quad_points.size() ; ++k){
	 	points.push_back({x_value , integration_params.mapped_quad_points[k]});
	 }
	 
	 std::vector<Vector<double>> values(points.size(), Vector<double>(2)); 
	 field.vector_value_list(points, values); 
	 
	 double sum = 0; 
	 
	 for(unsigned int k = 0 ; k < points.size() ; ++k){
	 	sum += (values[k](0) * values[k](0) + values[k](1) * values[k](1)) * integration_params.mapped_weights[k];
	 }
	 
	 return sum;
} // intensity

// intensity list
template<int dim> 
std::vector<IntensitySample> PostProcessor<dim>::intensity_list(const std::vector<double> &x_values) const {
	std::vector<IntensitySample> output; 
	
	for(unsigned int k = 0; k < x_values.size() ; ++k){
		output.push_back({x_values[k] , PostProcessor<dim>::intensity(x_values[k])});
	}
	return output;
	

} // intensity list

// implementing helper function transmission_reflection 

template<int dim> 
ModeCoefficients PostProcessor<dim>::transmission_reflection( const unsigned int mode) const {
	// compute beta
	const double k_0 = param.omega / param.c ; // k_0 
	const double cut_off = mode * numbers::PI / param.W ; // mode cut_off 
	const double beta_2 = k_0 * k_0 - cut_off * cut_off ; // beta_square
	
	AssertThrow(beta_2 > 0 , ExcMessage("The mode: " + std::to_string(mode) + " is evanescant!")); 
	double beta = std::sqrt(beta_2); // beta 
	
	//find the first and last planes incluing refractors
	double ref_min = 0.5 * (param.x_min + param.x_max);
	double ref_max = ref_min;
	if (!param.refractors.empty()) {
		ref_min = param.x_max;
		ref_max = param.x_min;
		for (const auto &r : param.refractors) {
			ref_min = std::min(ref_min, r.cx - r.radius);
			ref_max = std::max(ref_max, r.cx + r.radius);
		}
	}
	
	double gap_left = ref_min - param.x_min ; // gap between the input port and the first refractor 
	double gap_right = param.x_max - ref_max ; // gap between the output port and the second refractor 
	
	double d = numbers::PI/(4.0 * beta); // distance seperating the two planes where the coefficient will be calculated 
	
	AssertThrow(gap_left > d , ExcMessage("Not enough room left for two planes at the input")); 
	AssertThrow(gap_right > d , ExcMessage("Not enough room left for two planes at the output")); 
	
	double x_left_1 = param.x_min + 0.5 * (gap_left - d) ; // first plane at the left  
	double x_left_2 = x_left_1 + d ;  // second plane at the left 
	
	double x_right_1 = param.x_max - 0.5 * (gap_right - d) ; // first plane at the right  
	double x_right_2 = x_right_1 - d ; // second plane at the right 
	
	const std::complex<double> I (0.0 , 1.0) ; // the complex number i 
	
	// computing transmission and reflection at the input port  
	std::complex<double> a_left_1 = modal_amplitude(x_left_1 , mode); 
	std::complex<double> a_left_2 = modal_amplitude(x_left_2 , mode); 
	
	std::complex<double> em1 = std::exp(-I* beta * x_left_1), ep1 = std::exp(I * beta * x_left_1) ; 
	std::complex<double> em2 = std::exp(-I * beta * x_left_2), ep2 = std:: exp(I * beta * x_left_2); 
	
	std::complex<double> det = em1 * ep2 - em2 * ep1 ; 
	std::complex<double> A_left = (a_left_1 * ep2 - a_left_2 * ep1) / det ; 
	std::complex<double> B_left = (a_left_2 * em1 - a_left_1 * em2) / det ; 
	
	
	// computing transmisson and reflection and the output port 
	
	std::complex<double> a_right_1 = modal_amplitude(x_right_1 , mode); 
	std::complex<double> a_right_2 = modal_amplitude(x_right_2 , mode); 
	
	std::complex<double> emr1 = std::exp(-I* beta * x_right_1), epr1 = std::exp(I * beta * x_right_1) ; 
	std::complex<double> emr2 = std::exp(-I * beta * x_right_2), epr2 = std::exp(I * beta * x_right_2); 
	
	std::complex<double> detr = emr1 * epr2 - emr2 * epr1 ; 
	std::complex<double> A_right = (a_right_1 * epr2 - a_right_2 * epr1) / detr ; 
	std::complex<double> B_right = (a_right_2 * emr1 - a_right_1 * emr2) / detr ;
	
	
	
	ModeCoefficients result{mode , beta, A_left, B_left, A_right, B_right}; 
	return result ; 
	  
} // transmission reflection 



/*-- Implementing main funcitons --*/ 

template<int dim> 
std::vector<AmplitudeSample> PostProcessor<dim>::modal_amplitude_point_list() const{
	return PostProcessor<dim>::modal_amplitude_list(param.x_values);
}

template<int dim> 
std::vector<AmplitudeSample> PostProcessor<dim>::modal_amplitude_graph()const {
	const double L = param.x_max - param.x_min ; 
	std::vector<double> x_values; 
	for(unsigned int k = 0 ; k < param.n_points ; ++k){
		x_values.push_back(param.x_min + (k + 0.5) * L / param.n_points); 
	}
	
	return PostProcessor<dim>::modal_amplitude_list(x_values);
}

template<int dim> 
std::vector<IntensitySample> PostProcessor<dim>::intensity_point_list() const {
	return PostProcessor<dim>::intensity_list(param.x_values);
}

template<int dim> 
std::vector<IntensitySample> PostProcessor<dim>::intensity_graph()const {
	const double L = param.x_max - param.x_min ; 
	std::vector<double> x_values; 
	for(unsigned int k = 0 ; k < param.n_points ; ++k){
		x_values.push_back(param.x_min + (k + 0.5) * L / param.n_points); 
	}
	
	return PostProcessor<dim>::intensity_list(x_values);
}


template<int dim>
ScatteringResults PostProcessor<dim>::scattering_results() const {

	ModeCoefficients mode_coeff = transmission_reflection(param.n); 

	ScatteringResults res;
	res.r       = mode_coeff.B_left / mode_coeff.A_left;
	res.t       = mode_coeff.A_right / mode_coeff.A_left;
	res.R       = std::norm(res.r);
	res.T       = std::norm(res.t);
	res.balance = res.R + res.T;
	
	std::cout<<"Incoming wave amplitudes: A_left " << std::abs(mode_coeff.A_left)  << ". Reflected wave amplitude = " << std::abs(mode_coeff.B_left) << std::endl;
	return res;
}


template<int dim> 
std::vector<PowerSample> PostProcessor<dim>:: power_spectrum() const {
	
	const double k_0 = param.omega / param.c ; 
	const double cut_off = numbers::PI / param.W ; 
	const double ratio = k_0 / cut_off ; 
	
	AssertThrow(ratio > 1 , ExcMessage("The wave number is below the cut off"));

	const unsigned int max_mode = static_cast<unsigned int>(std::ceil(ratio)) - 1; 	
	
	ModeCoefficients incoming_wave = PostProcessor<dim>::transmission_reflection(param.n); 
	double incoming_power = incoming_wave.beta * std::norm(incoming_wave.A_left);
	
	std::vector<PowerSample> power_spec ; 
	
	
	
	for (unsigned int m = 1 ; m <= max_mode ; ++m) {
		ModeCoefficients mode_decom = PostProcessor<dim>::transmission_reflection(m); 
		double R_power = mode_decom.beta * std::norm(mode_decom.B_left) / incoming_power ; 
		double T_power = mode_decom.beta * std::norm(mode_decom.A_right) / incoming_power ; 
		double port_ref_in = ((mode_decom.mode == param.n)? 0.0 : 1.0) *  mode_decom.beta * std::norm(mode_decom.A_left) / incoming_power ; 
		double port_ref_out = mode_decom.beta * std::norm(mode_decom.B_right) / incoming_power; 
		power_spec.push_back({mode_decom.mode , mode_decom.beta , R_power , T_power , port_ref_in , port_ref_out}); 
	
	}
	
	return power_spec ; 
	
	

}



/*---------------------------------------- Pringting Function --------------------------------------------------*/

// printintg intensity computatoin results 
void write_intensity_csv(const std::vector<IntensitySample> &intensities , const std::string &output_file_name)  {

	std::ofstream output(output_file_name); 
	AssertThrow(output.is_open(), ExcMessage("Cannot open " + output_file_name)); 
	
	output << std::setprecision(16); 
	
	//write the columns headers 
	
	output << "x , Intensities" << std::endl ; 
	
	for (unsigned int k = 0 ; k < intensities.size() ; ++k){
		output << intensities[k].position << "," << intensities[k].value << std::endl; 
	}
	output.close();

}

// printing amplitude results 

void write_amplitude_csv(const Parameters &p ,const std::vector<AmplitudeSample> &amplitudes , const std::string &output_file_name)  {
	std::ofstream output(output_file_name); 
	AssertThrow(output.is_open(), ExcMessage("Cannot open " + output_file_name)); 
	output << std::setprecision(16); 
	
	// write colums headers 
	
	output << "mode , x , Re , Im , ABS" << std::endl; 
	
	for(unsigned int k = 0 ; k < amplitudes.size() ; ++k){
	output << p.m << "," << amplitudes[k].position << "," << amplitudes[k].value.real() << "," 
		<< amplitudes[k].value.imag() << "," << std::sqrt(std::norm(amplitudes[k].value)) << std::endl ;
	
	} 
	
	output.close() ;


}

// printing Scattering Results 

void write_scattering_results_csv (const ScatteringResults &res , const std::string &output_file_name)  {
	
	std::ofstream output(output_file_name); 
	AssertThrow(output.is_open(), ExcMessage("Cannot open " + output_file_name)); 
	output << std::setprecision(16); 
	
	output << "Re(r) , Im(r) , Re(t) , Im(t) , R , T , balance" << std::endl ;
	output << res.r.real() << "," << res.r.imag() << ","  
	       << res.t.real() << "," << res.t.imag() << "," 
	       << res.R << "," << res.T << "," 
	       << res.balance << std::endl ;
	output.close();

} 

// printing power spectrum 

void write_power_spectrum_csv(const std::vector<PowerSample> &power_spec, const std::string &output_file_name) {

	std::ofstream output(output_file_name); 
	AssertThrow(output.is_open() , ExcMessage("Cannot open" + output_file_name));
	output << std::setprecision(16); 
	
	output << "mode , beta , R , T , port_ref_in , port_ref_out" << std::endl;
	
	for(unsigned int m = 0 ; m < power_spec.size() ; ++m) {
		output<< power_spec[m].mode << "," << power_spec[m].beta << "," << power_spec[m].R << "," 
		<< power_spec[m].T << "," << power_spec[m].port_ref_in << "," << power_spec[m].port_ref_out << std::endl;
	}
	output.close();	

}















template class PostProcessor<2> ;

























