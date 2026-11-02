
#include "mesh_builder.h"
#include <gmsh.h> 
#include<set>
#include<vector>
#include<utility>
#include<iostream>
#include<fstream>

void build_mesh (const Parameters &p, const std::string &file_name) {
	namespace occ = gmsh::model::occ ; 
	
	gmsh::initialize(); 
	gmsh::option::setNumber("General.Terminal" , 0); 
	gmsh::model::add("waveguide"); 
	
	const int rect = occ::addRectangle(p.x_min, p.y_min,0, 
				     p.x_max - p.x_min , p.y_max - p.y_min); 
	
	std::vector<std::pair<int,int>> disks ; 
	
	for(const auto &r : p.refractors){
		disks.emplace_back(2, occ::addDisk(r.cx , r.cy, 0 , r.radius, r.radius)); 
		
	}
	
	std::vector<std::pair<int,int>> out; 
	std::vector<std::vector<std::pair<int,int>>> map; 
	occ::fragment({{2,rect}} , disks, out,map);
	occ::synchronize();
	
	// material ids : tag 1 = background, 2+i = refractor i 
	// map[0]holds all pices of the rectangle, map[i+1] the pieces of the disk i
	
	std::set<int> disk_tags; 
	for(std::size_t i = 0 ; i < disks.size() ; ++i) {
		std::vector<int> t ; 
		for (auto &e : map[i+1]) {
			t.push_back(e.second) ;
			disk_tags.insert(e.second); 
		}
		gmsh::model::addPhysicalGroup(2, t, 2+i); 
	}
	
	std::vector<int> outer; 
	for(auto &e : map[0]){
		if(!disk_tags.count(e.second)) outer.push_back(e.second);
	}
	
	gmsh::model::addPhysicalGroup(2, outer, 1 ); 
	
	
	//boundary ids found by bounding box: 1 = input port,2 = output port, 3 = walls 
	
	const double eps = 1e-6; 
	auto curves = [&] (double x0, double y0, double x1, double y1){
		std::vector<std::pair<int,int>> d ; 
		gmsh::model::getEntitiesInBoundingBox(x0- eps, y0-eps,-eps, x1 + eps , y1 + eps, eps, d,1); 
		std::vector<int> t; 
		for(auto &e : d) {
			t.push_back(e.second); 
			
		}
		return t;};
	
	
	//std::cout << "left "   << curves(p.x_min,p.y_min,p.x_min,p.y_max).size() used during debugging
          //<< " right  " << curves(p.x_max,p.y_min,p.x_max,p.y_max).size()
          //<< " bottom " << curves(p.x_min,p.y_min,p.x_max,p.y_min).size()
          //<< " top "    << curves(p.x_min,p.y_max,p.x_max,p.y_max).size()
          // << std::endl;
	gmsh::model::addPhysicalGroup(1, curves(p.x_min,p.y_min,p.x_min,p.y_max), 1);
  	gmsh::model::addPhysicalGroup(1, curves(p.x_max,p.y_min,p.x_max,p.y_max), 2);
  	auto walls = curves(p.x_min, p.y_min, p.x_max, p.y_min);        // bottom
	const auto top = curves(p.x_min, p.y_max, p.x_max, p.y_max);    // top
	walls.insert(walls.end(), top.begin(), top.end());
	gmsh::model::addPhysicalGroup(1, walls, 3);

	gmsh::option::setNumber("Mesh.MeshSizeMax", p.h_max);
	gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);
	gmsh::option::setNumber("Mesh.Algorithm", 6);                // Frontal-Delaunay
	gmsh::option::setNumber("Mesh.RecombineAll", 1);
	gmsh::option::setNumber("Mesh.RecombinationAlgorithm", 1);   // blossom
	gmsh::option::setNumber("Mesh.SubdivisionAlgorithm", 1);     // split into all quads
	gmsh::model::mesh::generate(2);
  	gmsh::write(file_name);
  	gmsh::finalize();
}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

