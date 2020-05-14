

#include <iomanip>
#include <string>
#include <memory>
#include <fstream>							 
#include <ctime>
#include <stdio.h>
#include <inttypes.h>
#include <cstddef>

#include "System.h"
#include "SystemBuilder.h"
#include "Storage.h"


#include "pugixml/include/pugixml.hpp"




std::shared_ptr<System> createSystem(const char* schemeFile, std::shared_ptr<SystemBuilder> builder)	{
	
	pugi::xml_document doc;
	pugi::xml_parse_result parseResult = doc.load_file(schemeFile);

	if (!parseResult) {
		std::cout << "parse error in createNodeSystem: " << parseResult.description() << std::endl;
		return nullptr;
	}
	pugi::xml_node root = doc.child("data");
	pugi::xml_node nodes = root.child("nodes"); 
	pugi::xml_node edges = root.child("edgeinfos");
	pugi::xml_node elems = root.child("elems");
	pugi::xml_node elem2edges = root.child("elem2edges");
	pugi::xml_node edge2elems = root.child("edge2elems");
	
	pugi::xml_node props = root.child("settings");

	//read in parameters from settings
	if (auto p = props.child("Tau")){
		builder->defaultTau = (p.text().as_double());
		std::cout<<"Setting tau: "<< builder->defaultTau << std::endl;
	} 

	if (auto p = props.child("KBT")){
		builder->defaultKBT = (p.text().as_double());
		std::cout<<"Setting kbt: "<< builder->defaultKBT << std::endl;
	}

	if (auto p = props.child("Linear_Const")){
		builder->defaultLinear_Const = (p.text().as_double());
		std::cout<<"Setting linear const: "<< builder->defaultLinear_Const << std::endl;
	}

	if (auto p = props.child("Area_Const")) {
		builder->defaultArea_Const = (p.text().as_double());
		std::cout<<"Setting area const: "<< builder->defaultArea_Const << std::endl;
	}

	if (auto p = props.child("Bend_Const")) {
		builder->defaultBending_Const = (p.text().as_double());
		std::cout<<"Setting bending const: "<< builder->defaultBending_Const << std::endl;
	}

	if (auto p = props.child("LJ_Eps")){ 
		builder->defaultLJ_Eps = (p.text().as_double());
		std::cout<<"Setting lj eps: "<< builder->defaultLJ_Eps << std::endl;
	}

	if (auto p = props.child("LJ_Rmin")){
		builder->defaultLJ_Rmin = (p.text().as_double());
		std::cout<<"Setting lj rmin: "<< builder->defaultLJ_Rmin << std::endl;
	}

	if (auto p = props.child("LJ_Rmax")){
		builder->defaultLJ_Rmax = (p.text().as_double());
		std::cout<<"Setting lj rmax: "<< builder->defaultLJ_Rmax << std::endl;
	}

	if (auto p = props.child("LJ_Const")){
		builder->defaultLJ_Const = (p.text().as_double());
		std::cout<<"Setting lj const: "<< builder->defaultLJ_Const << std::endl;
	}

	if (auto p = props.child("LJ_X")){
		builder->defaultLJ_X = (p.text().as_double());
		std::cout<<"Setting lj x: "<< builder->defaultLJ_X << std::endl;
	}

	if (auto p = props.child("LJ_Y")){
		builder->defaultLJ_Y = (p.text().as_double());
		std::cout<<"Setting lj y: "<< builder->defaultLJ_Y << std::endl;
	}

	if (auto p = props.child("LJ_Z")){
		builder->defaultLJ_Z = (p.text().as_double());
		std::cout<<"Setting lj z: "<< builder->defaultLJ_Z << std::endl;
	}

	//add nodes
	double x,y,z;
	for (auto node = nodes.child("node"); node; node = node.next_sibling("node")) {

		const char* text = node.text().as_string();

		//std::cout<<"mass: " << mass << std::endl;
		if (3 != sscanf(text, "%lf %lf %lf", &x, &y, &z)) {
			std::cout << "parse node error\n";
			return 0;
		}
		builder->addNode(x,y,z);
	}
	
	//add edges
	unsigned from, to;	
	for (auto edge = edges.child("edgeinfo"); edge; edge = edge.next_sibling("edgeinfo")) {
		if (2 != sscanf(edge.text().as_string(""), "%u %u" , &from, &to)) {
			std::cout << "parse link error\n";
			return 0;
		}
		//std::cout << "putting spring between: " << from << ' ' <<to<<  std::endl;
		builder->addEdge(from-1, to-1); //adds edges into saved vectors
	}
	
	//add triangles
	unsigned E1, E2, E3;	
	for (auto elem = elems.child("elem"); elem; elem = elem.next_sibling("elem")) {
		if (3 != sscanf(elem.text().as_string(""), "%u %u %u" , &E1, &E2, &E3)) {
			std::cout << "parse elem error\n";
			return 0;
		}
		//std::cout << "putting spring between: " << from << ' ' <<to<<  std::endl;
		builder->addElement(E1-1, E2-1, E3-1); //adds edges into saved vectors
	}

	//add indices of edges in each element
	for (auto elem2edge = elem2edges.child("elem2edge"); elem2edge; elem2edge = elem2edge.next_sibling("elem2edge")) {
		if (3 != sscanf(elem2edge.text().as_string(""), "%u %u %u" , &E1, &E2, &E3)) {
			std::cout << "parse elem2edge error\n";
			return 0;
		}
		//std::cout << "putting spring between: " << from << ' ' <<to<<  std::endl;
		builder->addElement2Edge(E1-1, E2-1, E3-1); //adds edges into saved vectors
	}
	
	//add indices of elems in each edge
	for (auto edge2elem = edge2elems.child("edge2elem"); edge2elem; edge2elem = edge2elem.next_sibling("edge2elem")) {
		if (2 != sscanf(edge2elem.text().as_string(""), "%u %u" , &E1, &E2)) {
			std::cout << "parse elem2edge error\n";
			return 0;
		}
		//std::cout << "putting spring between: " << from << ' ' <<to<<  std::endl;
		builder->addEdge2Elem(E1-1, E2-1); //adds edges into saved vectors
	}

	//set system variables
	std::shared_ptr<System> virus_system = builder->createSystem();
    return virus_system;
}


std::string generateOutputFileName(std::string inputFileName)
{
	time_t now;																   
	const int MAX_DATE = 64;
	char theDate[MAX_DATE];
	theDate[0] = '\0';

	now = time(nullptr);

	if (now != -1) {
		strftime(theDate, MAX_DATE, "_%Y.%m.%d_%H-%M-%S", gmtime(&now));
		return inputFileName + theDate;
	}
	return "";
}

void run(int argc, char** argv) {

	
	time_t t0,t1;
	t0 = time(0);

	double forceStep = 0.0;
	double timestep = 0.001;
	unsigned solve_time = 10000;
	bool time_found = false;

	for (int i = -1; i < argc-1; i++) {

		std::string arg = argv[i];
		unsigned pos = arg.find('=');

		std::string key = arg.substr(0, pos);
		std::string val = arg.substr(pos + 1);
		
		std::cout<<"argc: "<< argc <<std::endl;
		std::cout<<"arg: "<< arg <<std::endl;
		std::cout<<"pos: "<< pos <<std::endl;
		std::cout<<"key: "<< key <<std::endl;
		std::cout<<"val: "<< val <<std::endl;

		if (key == "-dt") {
			time_found = true;
			timestep = std::atof(val.c_str());
			std::cout<<"setting timestep: "<< timestep << std::endl;
			continue;
		}
		if (key == "-solve_time") {

			solve_time = std::atof(val.c_str());
			std::cout<<"setting solve time: "<< solve_time << std::endl;
			continue;
		}
	}

	if (time_found) {
		auto builder = std::make_shared<SystemBuilder>(timestep, solve_time);

		//sets all parameters and edges on device side
		
		std::shared_ptr<System> system = createSystem(argv[argc-1], builder);
		
		auto storage = std::make_shared<Storage>(system);
		
		system->assignStorage(storage);
		
		//std::cout << "post fdiagram in main" << std::endl;


		std::cout << "solving system in main" << std::endl;
		system->solveSystem();
	}
	t1 = time(0);  //current time at the end of solving the system.
	int total,hours,min,sec;
	total = difftime(t1,t0);
	hours = total / 3600;
	min = (total % 3600) / 60;
	sec = (total % 3600) % 60;
	// текущее время расчета
	
	std::cout << "Total time hh: " << hours << " mm:" << min << " ss:" << sec <<"\n";
}


int main(int argc, char** argv)
{
	std::cout << argc;
	std::cout << std::endl;

	run(argc - 2, argv + 2);
	
	return 0;
}
