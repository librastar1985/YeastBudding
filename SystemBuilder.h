
/*
 * SystemBuilder.h
 *
 *  Created on: 25 авг. 2014 г.
 *      Author: yan
 */

#ifndef SystemBuilder_H_
#define SystemBuilder_H_


#include "SystemStructures.h" 

//struct HostSetInfoVecs;
class System;


class SystemBuilder {
public:
	//set by constructor using command line
	double dt;
	int solve_time;

	//set by xml input file. 
	double defaultTau = 1.0; 
	double defaultKBT = 1.0;
	double defaultLinear_Const = 9.0;
	double defaultArea_Const = 10.0;
	double defaultBending_Const = 4.0 ;
	double defaultLJ_Eps = 0.1;
	double defaultLJ_Rmin = 2.0;
	double defaultLJ_Rmax = 2.0*1.4;
	double defaultLJ_Const = 1.0;
	double defaultLJ_X = 0.0;
	double defaultLJ_Y = 0.0;
	double defaultLJ_Z = -0.1;

	double defaultEdgeEq = 1.0;//0.1503;
	double defaultAreaEq = 0.0;
	double defaultAngleEq = 0.0;

	HostSetInfoVecs hostSetInfoVecs;


public:

	SystemBuilder(double timestep, int solve_time);
	//collection of set up host vectors in SystemStructures.h
	
	
	void addNode(double x, double y, double z);

	void addNndata(double x1,double x2, double x3, double x4,double x5, double x6, double x7,double x8, double x9, double x10,double x11, double x12 );

	void addEdge(int idL, int idR );

    void addEdge(int idL, int idR, double edge_initial_length);

 	void addElement(int idA, int idB, int idC );
 	
	void addElement2Edge(int idA, int idB, int idC );
 
	void addEdge2Elem(int idA, int idB );

	void fixNodes(int id);
	
	void addCapsidNode(double x, double y, double z);
	//void setSystemForParallelComputation();
	std::shared_ptr<System> createSystem();


};

#endif /* SystemBuilder_H_ */
