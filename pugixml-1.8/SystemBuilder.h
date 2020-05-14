
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
	unsigned solve_time;

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

	HostSetInfoVecs hostSetInfoVecs;


public:

	SystemBuilder(double timestep, unsigned solve_time);
	//collection of set up host vectors in SystemStructures.h
	
	
	void addNode(double x, double y, double z);

	void addEdge(unsigned idL, unsigned idR );

 	void addElement(unsigned idA, unsigned idB, unsigned idC );
 	
	void addElement2Edge(unsigned idA, unsigned idB, unsigned idC );
 
	void addEdge2Elem(unsigned idA, unsigned idB );
	//void setSystemForParallelComputation();
	std::shared_ptr<System> createSystem();


};

#endif /* SystemBuilder_H_ */
