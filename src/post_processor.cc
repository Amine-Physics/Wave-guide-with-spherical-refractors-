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

	const double W       = param.y_max - param.y_min;
	const double k0      = param.omega / param.c;
	const double cutoff  = param.n * numbers::PI / W;
	const double beta2   = k0 * k0 - cutoff * cutoff;
	AssertThrow(beta2 > 0, ExcMessage("The mode is evanescent"));
	const double beta = std::sqrt(beta2);

	// x-extent of the refractors (guide midpoint if there are none)
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

	// two planes left of the refractors, a quarter guided wavelength apart
	// (best conditioning for separating the e^{-i beta x} and e^{+i beta x} parts)
	const double d   = numbers::PI / (4.0 * beta);
	const double gap = ref_min - param.x_min;
	AssertThrow(gap > d, ExcMessage("Not enough room left of the refractors for the two planes"));
	const double x1      = param.x_min + 0.5 * (gap - d);
	const double x2      = x1 + d;
	const double x_right = 0.25 * (3*ref_max + param.x_max);

	const std::complex<double> I(0.0, 1.0);
	const std::complex<double> a1    = modal_amplitude(x1, param.n);
	const std::complex<double> a2    = modal_amplitude(x2, param.n);
	const std::complex<double> a_out = modal_amplitude(x_right, param.n);

	// solve a_j = A e^{-i beta x_j} + B e^{+i beta x_j}, j = 1,2
	const std::complex<double> em1 = std::exp(-I * beta * x1), ep1 = std::exp(I * beta * x1);
	const std::complex<double> em2 = std::exp(-I * beta * x2), ep2 = std::exp(I * beta * x2);
	const std::complex<double> det = em1 * ep2 - ep1 * em2;
	const std::complex<double> A   = (a1 * ep2 - a2 * ep1) / det;   // incident
	const std::complex<double> B   = (em1 * a2 - em2 * a1) / det;   // reflected

	ScatteringResults res;
	res.r       = B / A;
	res.t       = a_out * std::exp(I * beta * x_right) / A;
	res.R       = std::norm(res.r);
	res.T       = std::norm(res.t);
	res.balance = res.R + res.T;
	std::cout<<"Reflacted Wave Amplitudes: A = " << std::abs(A)  << ", B = " << std::abs(B) << std::endl;
	std::cout<< "Error( |A| - sqrt(2/W)) = " << std::abs(A) - std::sqrt(2/W) << std::endl;
	return res;
}

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

// printin amplitude results 

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

// print ScatteringResults 

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















template class PostProcessor<2> ;

























