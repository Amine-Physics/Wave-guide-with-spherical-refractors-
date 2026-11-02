/*----------------------------------------------------------------------------------------------------

This the declaration of the Dirclet Boundary Condition class. It provides information about the value 
of the wave at the wave guide walls. 
----------------------------------------------------------------------------------------------------*/


#pragma once 
#include<deal.II/base/function.h>
#include "parameters.h"
#include<deal.II/lac/vector.h> 
#include<vector> 
#include<deal.II/base/point.h>
using namespace dealii;
template<int dim> 

class DirchletBoundaryCondition : public dealii::Function<dim> {
	public: DirchletBoundaryCondition () ; 
		virtual void vector_value( const Point<dim> &p , Vector<double> &values) const override; 
		virtual void vector_value_list( const std::vector<Point<dim>> &points , std::vector<Vector<double>> &values_list) const override ;  

};

// implementation 

template<int dim> 
DirchletBoundaryCondition<dim>::DirchletBoundaryCondition(): dealii::Function <dim> (2) {}

template<int dim> 
void DirchletBoundaryCondition<dim>::vector_value( const Point<dim> &p, Vector<double> &values) const {
	(void) p ;
	AssertDimension(values.size() , 2 ); 
	values[0] = 0 ; 
	values[1] = 0 ; 
}

template<int dim> 
void DirchletBoundaryCondition<dim>:: vector_value_list (const std::vector<Point<dim>> &points , std::vector<Vector<double>> &values_list) const {
	AssertDimension(points.size() , values_list.size()) ; 
	for (unsigned int i = 0 ; i < points.size() ; ++i) {
		DirchletBoundaryCondition<dim>::vector_value(points[i] , values_list[i]);
	}

} 