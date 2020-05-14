//#include <curand.h>
#include <set>
#include <list>					 
#include <vector>
#include <memory>
#include "System.h"
#include "SystemBuilder.h"
#include "SystemStructures.h"
# define M_PI 3.14159265358979323846  /* pi */

SystemBuilder::SystemBuilder(double timestep_, unsigned solve_time_){
	dt = timestep_;
	solve_time = solve_time_;
};

void SystemBuilder::addNode(double x,double y, double z ) {
	hostSetInfoVecs.hostNodeLocX.push_back(x);
	hostSetInfoVecs.hostNodeLocY.push_back(y);
	hostSetInfoVecs.hostNodeLocZ.push_back(z);

	hostSetInfoVecs.hostIsNodeFixed.push_back(false);


	//std::cout<<"adding node: "<< x << " " << y << " "<< z <<  std::endl;
}

void SystemBuilder::addEdge(unsigned idL, unsigned idR ) {
	hostSetInfoVecs.hostEdges2Nodes_1.push_back(idL);
	hostSetInfoVecs.hostEdges2Nodes_2.push_back(idR);


	double xL = hostSetInfoVecs.hostNodeLocX[idL];
	double yL = hostSetInfoVecs.hostNodeLocY[idL];
	double zL = hostSetInfoVecs.hostNodeLocZ[idL];
	double xR = hostSetInfoVecs.hostNodeLocX[idR];
	double yR = hostSetInfoVecs.hostNodeLocY[idR];
	double zR = hostSetInfoVecs.hostNodeLocZ[idR];
	double dist = std::sqrt( (xL-xR)*(xL-xR) + (yL-yR)*(yL-yR) + (zL-zR)*(zL-zR));
	hostSetInfoVecs.hostEdge_initial_length.push_back(dist);
	//std::cout<<"adding edge: "<< idL << " " << idR << std::endl;
}

void SystemBuilder::addElement(unsigned idA, unsigned idB, unsigned idC ) {
	hostSetInfoVecs.hostTriangles2Nodes_1.push_back(idA);
	hostSetInfoVecs.hostTriangles2Nodes_2.push_back(idB);	
	hostSetInfoVecs.hostTriangles2Nodes_3.push_back(idC);
	
}

void SystemBuilder::addElement2Edge(unsigned idA, unsigned idB, unsigned idC ) {
	hostSetInfoVecs.hostTriangles2Edges_1.push_back(idA);
	hostSetInfoVecs.hostTriangles2Edges_2.push_back(idB);	
	hostSetInfoVecs.hostTriangles2Edges_3.push_back(idC);
	
}

void SystemBuilder::addEdge2Elem(unsigned idA, unsigned idB ) {
	hostSetInfoVecs.hostEdges2Triangles_1.push_back(idA);
	hostSetInfoVecs.hostEdges2Triangles_2.push_back(idB);	
}


//adds all constraints to the nodesystem model so that it can use the constraints.
std::shared_ptr<System> SystemBuilder::createSystem() {


	//now all the edges and variables are set. 
	//so set the system and return a pointer.
	std::shared_ptr<System> host_ptr_System = std::make_shared<System>();
	//hand in individually set parameters here:


	//set individual parameters
	host_ptr_System->generalParams.dt = dt;
	host_ptr_System->generalParams.solve_time = solve_time;//unsigned
	host_ptr_System->generalParams.tau = defaultTau; 
	host_ptr_System->generalParams.kT = defaultKBT;

	host_ptr_System->linearSpringInfoVecs.spring_constant = defaultLinear_Const;
	host_ptr_System->areaTriangleInfoVecs.spring_constant = defaultArea_Const;
	host_ptr_System->bendingTriangleInfoVecs.spring_constant = defaultBending_Const;

	host_ptr_System->ljInfoVecs.epsilon = defaultLJ_Eps;
	host_ptr_System->ljInfoVecs.Rmin = defaultLJ_Rmin;
	host_ptr_System->ljInfoVecs.Rcutoff = defaultLJ_Rmax;
	host_ptr_System->ljInfoVecs.spring_constant = defaultLJ_Const;
	host_ptr_System->ljInfoVecs.LJ_PosX = defaultLJ_X;
	host_ptr_System->ljInfoVecs.LJ_PosY = defaultLJ_Y;
	host_ptr_System->ljInfoVecs.LJ_PosZ = defaultLJ_Z;


	//set vectors and allocate memory here:
	host_ptr_System->initializeSystem(
		hostSetInfoVecs
	);
	
	
	return host_ptr_System;							 

}

