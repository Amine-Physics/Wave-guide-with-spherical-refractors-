/*This is the main header file for the FEM solver here we declare the class WaveGuideProblem with takes
 in a parameter file and produces the solution for the the propagation of the EM wave. */

#pragma once 

#include<deal.II/base/function.h> 
#include<deal.II/base/quadrature_lib.h> 
#include<deal.II/base/timer.h> 
#include<deal.II/base/utilities.h>
#include<deal.II/base/numbers.h> 
#include<deal.II/base/hdf5.h>

#include<deal.II/lac/dynamic_sparsity_pattern.h> 
#include<deal.II/lac/sparse_matrix.h> 
#include<deal.II/lac/full_matrix.h> 
#include<deal.II/lac/vector.h> 
#include<deal.II/lac/sparse_direct.h> 

#include<deal.II/dofs/dof_handler.h> 
#include<deal.II/dofs/dof_tools.h> 

#include<deal.II/fe/fe_system.h> 
#include<deal.II/fe/fe_q.h> 
#include<deal.II/fe/fe_values.h> 

#include<deal.II/grid/grid_in.h> 
#include<deal.II/grid/grid_out.h> 
#include<deal.II/grid/manifold_lib.h>
#include<deal.II/grid/tria.h> 

#include<deal.II/numerics/data_out.h> 
#include<deal.II/numerics/vector_tools.h> 
#include<deal.II/numerics/vector_tools.h> 
#include<deal.II/numerics/matrix_tools.h>

#include<iostream> 
#include<fstream> 
#include<cmath> 

#include "parameters.h" 
#include "mesh_builder.h"
#include "dirchlet_boundary_condition.h"
#include "post_processor.h"


using namespace dealii;

template<int dim> 
class WaveGuideProblem {
	public: WaveGuideProblem(Parameters &p); 
	        void run(); 
	
	private: void make_grid (); 
	         void setup_system(); 
	         void assemble_system(); 
	         void solve (); 
	         void output_results() const ; 
	         void output_postprocessing_results() const ;
	         
	         Parameters &param; 
	         Triangulation<dim> triangulation ; 
	         DoFHandler<dim> dof_handler ; 
	         const FESystem<dim>  fe ; 
	         
	         SparsityPattern sparsity_pattern ; 
	         SparseMatrix<double> system_matrix; 
	         Vector<double> solution, system_rhs; 
		
}; 