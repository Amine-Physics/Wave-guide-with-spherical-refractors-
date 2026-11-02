#include"wave_guide_problem.h"

int main (int argc , char* argv[]){

	using namespace dealii;
	Utilities::MPI::MPI_InitFinalize mpi_initilization(argc , argv);
	const int dim = 2; 
	
	Parameters p = read_parameters("wave_guide_problem.prm"); 
	WaveGuideProblem<dim> wave_guide_problem_2_D(p); 
	wave_guide_problem_2_D.run();

}

