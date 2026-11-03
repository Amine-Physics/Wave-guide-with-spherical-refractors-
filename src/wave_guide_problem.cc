/*This is the source file for the WaveGuideProblem class it implements the different functions*/

#include "wave_guide_problem.h" 

using namespace dealii;


template<int dim> 
WaveGuideProblem<dim>::WaveGuideProblem(Parameters &param_h) : param(param_h), 
						dof_handler(triangulation), 
						fe(FE_Q<dim>(1) ^2 ) {} 
						

template<int dim> 
void WaveGuideProblem<dim>::make_grid(){
	
	std::cout <<"Generating the Grid ..." << std::endl; 
	Timer timer; 
	
	build_mesh(param , "grid.msh"); // generating the mesh using Gmsh
	
	GridIn<dim> grid_in; 
	grid_in.attach_triangulation(triangulation); 
	std::ifstream input_grid_file ("grid.msh"); 
	Assert(dim == 2, ExcNotImplemented());
	grid_in.read_msh(input_grid_file); 
	
	//set the manifold id 
	for( auto &cell : triangulation.active_cell_iterators()) {
		for(unsigned int f : cell->face_indices()){
			if (!cell->face(f)->at_boundary() && cell->material_id() != cell->neighbor(f)->material_id()){ // a face is at the boundary between the 
											//background and a refractor if the material id changes when we cross the face
				unsigned int refractor_id = std::max(cell->material_id() , cell->neighbor(f)->material_id()) - 2; 
				cell->face(f)->set_manifold_id(refractor_id); 
				cell->set_manifold_id(refractor_id); 
				cell->neighbor(f)->set_manifold_id(refractor_id);
				
			}
		}
	}
	
	for(unsigned int ref = 0 ; ref < param.refractors.size() ; ++ref){
		const Point<dim> ref_center (param.refractors[ref].cx , param.refractors[ref].cy);
		triangulation.set_manifold(ref , SphericalManifold<dim>(ref_center));
	}
	
	triangulation.refine_global(param.n_refinements);
	
	timer.stop(); 
	std::cout<<"done (" << timer.cpu_time() << "s )" << std::endl; 
	std::cout<<"Number of active cells :" << triangulation.n_active_cells()<< std::endl;
	
	// reading boundary ids //used during debugging  
	//for (const auto id : triangulation.get_boundary_ids()){
 		//std::cout << "boundary id " << static_cast<int>(id) << std::endl;}
 		
	//outputing grid in svg format 
	GridOut grid_out ; 
	std::string grid_svg = "grid.svg"; 
	std::ofstream grid_output(grid_svg); 
	grid_out.write_svg(triangulation, grid_output);
} // make_grid 


template<int dim> 

void WaveGuideProblem<dim>::setup_system() {
	
	std::cout<<"Setting up the system ... " << std::endl ; 
	
	Timer timer ; 
	
	dof_handler.distribute_dofs(fe); // distribute dofs
	
	DynamicSparsityPattern dsp(dof_handler.n_dofs() , dof_handler.n_dofs()); 
	DoFTools::make_sparsity_pattern(dof_handler , dsp);
	sparsity_pattern.copy_from(dsp);
	
	
	system_matrix.reinit(sparsity_pattern); // initilize system_matrix
	
	solution.reinit(dof_handler.n_dofs()); 
	system_rhs.reinit(dof_handler.n_dofs());
	
	timer.stop(); 
	
	std::cout<<"done ("<< timer.cpu_time() << "s )" << std::endl; 
	std::cout<<"Number of degrees of freedom:" << dof_handler.n_dofs() << std::endl;

}


template<int dim> 

void WaveGuideProblem<dim>::assemble_system() {

	std::cout<<"Assembling the system ..." << std::endl; 
	Timer timer ; 
	
	// setting up quadrature rule for cells and faces
	QGauss<dim> quadrature_rule(fe.degree + 1) ; 
	QGauss<dim - 1> face_quadrature_rule(fe.degree + 1); 
	
	// setting up fe_values for cells and faces
	FEValues<dim> fe_values(fe, quadrature_rule , update_values |update_gradients |update_JxW_values ) ; 
	FEFaceValues<dim> fe_face_values(fe,face_quadrature_rule, update_values|update_JxW_values|update_quadrature_points); 
	
	// usefull variables and data points 
	
	const double GW_width = param.y_max - param.y_min ; // wave guide width
	const double beta = std::sqrt(std::pow(param.omega / param.c , 2) - std::pow(param.n * numbers::PI / GW_width , 2)) ; // incoming wave propagation constant 
	double eps_r = 0; 
	const unsigned int n_q_points = quadrature_rule.size(), 
		         n_face_q_points = face_quadrature_rule.size(), 
		         dofs_per_cell = fe.n_dofs_per_cell(); 
		         
	FullMatrix<double> cell_matrix(dofs_per_cell , dofs_per_cell); // local matrix
	Vector<double> cell_rhs(dofs_per_cell); // local rhs
	std::vector<types::global_dof_index> local_dofs_indices (dofs_per_cell); // connecting local to global dofs indices 
	
	// main loop over cells 
	
	for (const auto & cell : dof_handler.active_cell_iterators()){
		// initilizing local data structures 
		
		fe_values.reinit(cell); 
		cell_matrix = 0; 
		cell_rhs = 0; 
		eps_r = 0;
		//determinig relative permitivity index 
		if(cell->material_id() !=1){
			eps_r = param.refractors[cell->material_id() - 2].eps ; 
		}
		else{
			eps_r = 1; 
						
		}
					
		// first loop over dofs 
		for(unsigned int i = 0 ; i < dofs_per_cell ; ++i) {
			// second loop over dofs
			for(unsigned int j = 0 ; j < dofs_per_cell ; ++j) {
				
				//cheking both dofs coorespond the the same compoenent (real or imaginary part 
				
				if(fe.system_to_component_index(i).first == fe.system_to_component_index(j).first){
					//loop over quadrature points 
					
					for(unsigned int q = 0; q < n_q_points ; q++){
					
						cell_matrix(i,j) +=(( - std::pow(param.omega , 2) * fe_values.shape_value(i,q) * fe_values.shape_value(j,q) ) + 
							   (std::pow(param.c , 2) * fe_values.shape_grad(i,q) * (1.0/eps_r) * fe_values.shape_grad(j,q))) * 
							   fe_values.JxW(q); 
					
					}// loop over quadrature points 
					
				} //cheking both dofs coorespond the the same compoenent (real or imaginary part 
			} // second loop over dofs
			
		} // first loop over dofs 
		
		// loop over cell faces for including port boundary conditions 
		for(const auto &face_no: cell->face_indices()){
			
			// cheking that the face is at the boudary at the boundary 
			if(cell->face(face_no)->at_boundary() && cell->face(face_no)->boundary_id() !=3){
				fe_face_values.reinit(cell,face_no);
				
				//loop over first dofs index 
				for(unsigned int i =0 ; i < dofs_per_cell ; ++i) {
					// cheking for input or output to update cell_rhs
					if (cell->face(face_no)->boundary_id() == 1) {
						// loop over quadrature points 
						for(unsigned int q_face = 0 ; q_face < n_face_q_points ; ++q_face){
							Point<dim> quad_point = fe_face_values.get_quadrature_points()[q_face];
						   	cell_rhs(i) += -2 * param.c *param.c *beta* std::sqrt(2.0 / param.W) * ((fe.system_to_component_index(i).first == 0)? -std::cos(beta * quad_point[0]) : 
						   	std::sin(beta * quad_point[0])) * std::sin(numbers::PI * (quad_point[1] - param.y_min)/ GW_width) *
						   	fe_face_values.shape_value(i,q_face) * fe_face_values.JxW(q_face);
						 }
					}
					// loop over second dofs index
					for (unsigned int j = 0; j < dofs_per_cell ; ++j){
					
						//cheking that dofs belong to deferent components 
						if((fe.system_to_component_index(i).first != fe.system_to_component_index(j).first) && 
						   (fe.has_support_on_face(i,face_no))&&
						   (fe.has_support_on_face(j,face_no))){
						   	
						   	// loop over quadrature points 
						   	for(unsigned int q_face = 0 ; q_face < n_face_q_points; ++q_face){
						   		cell_matrix(i,j) += ((fe.system_to_component_index(i).first == 0) ? -1 : 1)  * beta * param.c * param.c * 
						   				fe_face_values.shape_value(i,q_face) * fe_face_values.shape_value(j,q_face) * fe_face_values.JxW(q_face);
						   		
						   	} // loop over face quadrature points
						   
						   }//cheking that dofs belong to deferent components
						
						
					}// loop over second dofs index
				}// loop over first dofs index 
			
			}// cheking that the face is at the boudary at the boundary
				
		
		} // loop over cell faces
		
		// felling system_matrix and system_rhs 
		
		cell->get_dof_indices(local_dofs_indices); 
		
		// loop over first dof index 
		for (unsigned int i = 0; i < dofs_per_cell ; ++i){
			system_rhs(local_dofs_indices[i]) += cell_rhs(i);
			//loop over second dof index 
			for (unsigned int j = 0 ; j < dofs_per_cell ; ++j){
				system_matrix.add(local_dofs_indices[i], local_dofs_indices[j] , cell_matrix(i,j));
			}
		}
		
	} // main loop over cells 
	
	std::map<types::global_dof_index , double> boundary_values ; 
	VectorTools::interpolate_boundary_values(dof_handler, 3 /*boundary_id*/, DirchletBoundaryCondition<dim>(), boundary_values); 
	MatrixTools::apply_boundary_values(boundary_values, system_matrix, solution, system_rhs); 
		
	timer.stop(); 
	std::cout<<"done(" << timer.cpu_time() <<"s)" << std::endl; 
}

//solve
template<int dim> 
void WaveGuideProblem<dim>:: solve(){
	std::cout<<"Solving the system .." << std::endl; 
	Timer timer; 
		
	SparseDirectUMFPACK A_direct;  
	solution = system_rhs ; 
	A_direct.solve(system_matrix,solution); 
		
	timer.stop(); 
	std::cout<<"done ("<< timer.cpu_time() << "s )" << std::endl ; 
	
}// solve
	
//output results 
	
template<int dim> 
void WaveGuideProblem<dim>::output_results() const {
	std::cout<<"Generating output files ..." << std::endl; 
		
	Timer timer ; 
	
		
	DataOut<dim> data_out ; 
	data_out.attach_dof_handler(dof_handler); 
 		
	std::string output_file_name_hdf5 = param.output_file_name + "_" +  std::to_string(param.refractors.size()) +  ".h5" ; 
	std::string output_file_name_vtu = param.output_file_name +"_"+  std::to_string(param.refractors.size()) + ".vtu" ;
	std::ofstream output (output_file_name_vtu); 
	
	DataOutBase::DataOutFilterFlags flags(true , true);
	DataOutBase::DataOutFilter filter(flags);	
	std::vector<std::string> solution_names ; 
	solution_names.emplace_back("Re"); 
	solution_names.emplace_back("IM");
		
	data_out.add_data_vector(solution , solution_names); 
	
	Vector<float> material(triangulation.n_active_cells());
	unsigned int k = 0;
	for (const auto &cell : triangulation.active_cell_iterators())
  		material(k++) = cell->material_id();

	data_out.add_data_vector(material, "material_id");
		
	data_out.build_patches(); 
	
	data_out.write_filtered_data(filter);
	data_out.write_hdf5_parallel(filter , output_file_name_hdf5 , MPI_COMM_WORLD);
	data_out.write_vtu(output);
	timer.stop(); 
		
	std::cout<<"done (" << timer.cpu_time() << " s)" << std::endl; 
	std::cout<< "The solution has been writen to: " << output_file_name_hdf5 << " as an h5 file "<< std::endl 
	<< "and to: " << output_file_name_vtu << " as a vtu file " << std::endl;
	
}

template<int dim> 
void WaveGuideProblem<dim>:: output_postprocessing_results() const {
	
	std::cout<<"Postprocessing ..." << std::endl ;
	Timer timer ; 
	
	Integration inte ; 
	get_integration_points_and_weights(param , inte); 
	
	PostProcessor<dim> post_processor(dof_handler, solution , param, inte); 
	
	std::vector<AmplitudeSample> point_amplitudes = post_processor.modal_amplitude_point_list() ; 
	std::vector<AmplitudeSample> graph_amplitudes = post_processor.modal_amplitude_graph(); 
	
	std::vector<IntensitySample> point_intensities = post_processor.intensity_point_list();
	std::vector<IntensitySample> graph_intensities = post_processor.intensity_graph(); 
	
	ScatteringResults res = post_processor.scattering_results(); 
	
	std::vector<PowerSample> power_spec= post_processor.power_spectrum();
	
	std::string point_amplitude_file = "point_amplitude.csv" ; 
	std::string graph_amplitude_file ="graph_amplitude.csv" ; 
	std::string point_intensity_file ="point_intensity.csv" ; 
	std::string graph_intensity_file ="graph_intensity.csv" ; 
	std::string scattering_results_file = "scattering_results.csv" ; 
	std::string power_spectrum_file = "power_spectrum.csv" ; 
	
	write_intensity_csv(point_intensities , point_intensity_file);
	write_intensity_csv(graph_intensities , graph_intensity_file);
	
	write_amplitude_csv(param, point_amplitudes , point_amplitude_file);
	write_amplitude_csv(param, graph_amplitudes , graph_amplitude_file);
	
	write_scattering_results_csv(res , scattering_results_file);
	
	write_power_spectrum_csv(power_spec , power_spectrum_file) ;
	
	timer.stop(); 
	std::cout<< "done (" << timer.cpu_time() << "s )" << std::endl; 
	
	std::cout<<"----------Scattering Results-----------" << std::endl; 
	std::cout<< "   R      |      T   |   B      " << std::endl;
	std::cout<< " " << res.R << " | "  << res.T << " | " << res.balance << std::endl;
}


	


template<int dim> 
void WaveGuideProblem<dim>::run(){
	make_grid();
	setup_system();
	assemble_system();
	std::cout<< "System right hand side norm:" << system_rhs.l2_norm() << std::endl;
	solve();
	output_results();
	output_postprocessing_results();
}


 template class WaveGuideProblem<2>;
