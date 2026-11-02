#include "parameters.h"
#include "mesh_builder.h" 
#include <iostream> 
#include <fstream> 
#include<gmsh.h> 

int main (){
	Parameters p = read_parameters("mesh_test_param_1_ref.prm"); 
	build_mesh(p , "test.mesh"); 
	
	  // reopen and inspect
  	gmsh::initialize();
  	gmsh::option::setNumber("General.Terminal", 1);
  	gmsh::open("test.msh");

  	// physical groups: expect surfaces 1 (background) and 2 (refractor),
  	// curves 1 (input), 2 (output), 3 (walls)
  	std::vector<std::pair<int, int>> groups;
  	gmsh::model::getPhysicalGroups(groups);
  	for (const auto &g : groups)
    	{
      	std::vector<int> ents;
      	gmsh::model::getEntitiesForPhysicalGroup(g.first, g.second, ents);
      	std::cout << "dim " << g.first << ", group " << g.second << ": "
                	<< ents.size() << " entities" << std::endl;
    	}

  	// 2D element types: expect only type 3 (4-node quadrangle)
  	std::vector<int> types;
  	std::vector<std::vector<std::size_t>> tags, nodes;
  	gmsh::model::mesh::getElements(types, tags, nodes, 2);
  	for (std::size_t i = 0; i < types.size(); ++i)
   	 std::cout << "2D element type " << types[i] << ": " << tags[i].size()
              << " elements"<<std::endl;

 	 gmsh::finalize();
}

