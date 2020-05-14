

#ifndef SYSTEM_H_
#define SYSTEM_H_

#pragma once

#include <memory>
#include <math.h>  
#include <thrust/extrema.h>
#include <thrust/sort.h>
#include <thrust/copy.h>
#include <thrust/for_each.h>
#include <thrust/transform.h> 
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/remove.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/copy.h>
#include <thrust/execution_policy.h> 
#include <thrust/pair.h>
#include <stdint.h>

//Data Structure for node location. velocity and force
struct CoordInfoVecs {



	thrust::device_vector<bool> isNodeFixed;
	//GLOBAL COORDS
	// X,Y,Z, location, velocity and force of all nodes
	thrust::device_vector<double> prevNodeLocX;
	thrust::device_vector<double> prevNodeLocY;
	thrust::device_vector<double> prevNodeLocZ;	

	thrust::device_vector<double> prevNodeForceX;
	thrust::device_vector<double> prevNodeForceY;
	thrust::device_vector<double> prevNodeForceZ;

 	thrust::device_vector<double> nodeLocX;
	thrust::device_vector<double> nodeLocY;
	thrust::device_vector<double> nodeLocZ;
	
	thrust::device_vector<double> nodeForceX;
	thrust::device_vector<double> nodeForceY;
	thrust::device_vector<double> nodeForceZ;
	
	//LOCAL COORDS
	//indices of each triangle
	unsigned num_triangles;
	thrust::device_vector<unsigned> triangles2Nodes_1;
	thrust::device_vector<unsigned> triangles2Nodes_2;
	thrust::device_vector<unsigned> triangles2Nodes_3;

	
	//indices of each edge
	unsigned num_edges;
	thrust::device_vector<unsigned> edges2Nodes_1;
	thrust::device_vector<unsigned> edges2Nodes_2;

	//indices of 2 triangle on each edge
	thrust::device_vector<unsigned> edges2Triangles_1;
	thrust::device_vector<unsigned> edges2Triangles_2;

	//indices of edges on each triangle.
	thrust::device_vector<unsigned> triangles2Edges_1;
	thrust::device_vector<unsigned> triangles2Edges_2;
	thrust::device_vector<unsigned> triangles2Edges_3;
	

};
 


//struct used for linking of nodes in network 
struct AuxVecs {
	// bucket key means which bucket ID does a certain point fit into
	thrust::device_vector<unsigned> bucketKeys;	//bucket id
	// bucket value means global rank of a certain point
	thrust::device_vector<unsigned> bucketValues;//node id
	// bucket key expanded means what are the bucket IDs are the neighbors of a certain point
	thrust::device_vector<unsigned> bucketKeysExpanded;
	// bucket value expanded means each point ( represented by its global rank) will have multiple copies
	thrust::device_vector<unsigned> bucketValuesIncludingNeighbor;
	
	// begin position of a keys in bucketKeysExpanded and bucketValuesIncludingNeighbor
	//entry keyBegin[bucketKey] returns start of indices to link
	thrust::device_vector<unsigned> keyBegin;
	// end position of a keys in bucketKeysExpanded and bucketValuesIncludingNeighbor
	thrust::device_vector<unsigned> keyEnd;
	
	unsigned endIndexBucketKeys; 
};



struct DomainParams {
	double minX;
	double maxX;
	double minY;
	double maxY;
	double minZ;
	double maxZ;
	double originMinX;
	double originMaxX;
	double originMinY;
	double originMaxY;
	double originMinZ;
	double originMaxZ;
	double gridSpacing = 0.5;
	unsigned XBucketCount;
	unsigned YBucketCount;
	unsigned ZBucketCount;
	unsigned totalBucketCount;
};

struct LJInfoVecs{
	double LJ_PosX;
	double LJ_PosY;
	double LJ_PosZ;
	double Rmin;
	double Rcutoff;

	double epsilon;
	double spring_constant;
	
	thrust::device_vector<unsigned> node_id_close;
	double lj_energy;
	double forceX;
	double forceY;
	double forceZ;

};

struct AreaTriangleInfoVecs {

	unsigned factor = 3;//used for reduction
	double initial_area = 0.433;
	double spring_constant;

	double area_triangle_energy;

	thrust::device_vector<unsigned> tempNodeIdUnreduced;
	thrust::device_vector<double> tempNodeForceXUnreduced;
	thrust::device_vector<double> tempNodeForceYUnreduced;
	thrust::device_vector<double> tempNodeForceZUnreduced;
	
	thrust::device_vector<unsigned> tempNodeIdReduced;
	thrust::device_vector<double> tempNodeForceXReduced;
	thrust::device_vector<double> tempNodeForceYReduced; 
	thrust::device_vector<double> tempNodeForceZReduced;

};

struct BendingTriangleInfoVecs {
	unsigned numBendingSprings=0;

	unsigned factor = 4;//used for reduction
	double spring_constant = 4.0;
	double initial_angle = 0.0;//radians

	double bending_triangle_energy;

	thrust::device_vector<unsigned> tempNodeIdUnreduced;
	thrust::device_vector<double> tempNodeForceXUnreduced;
	thrust::device_vector<double> tempNodeForceYUnreduced;
	thrust::device_vector<double> tempNodeForceZUnreduced;
	
	thrust::device_vector<unsigned> tempNodeIdReduced;
	thrust::device_vector<double> tempNodeForceXReduced;
	thrust::device_vector<double> tempNodeForceYReduced;
	thrust::device_vector<double> tempNodeForceZReduced;
};
struct LinearSpringInfoVecs {
	
	unsigned factor = 2;//used for reduction
	double spring_constant;
	
	double linear_spring_energy;

	thrust::device_vector<double> edge_initial_length;

	thrust::device_vector<unsigned> tempNodeIdUnreduced;
	thrust::device_vector<double> tempNodeForceXUnreduced;
	thrust::device_vector<double> tempNodeForceYUnreduced;
	thrust::device_vector<double> tempNodeForceZUnreduced;
	
	thrust::device_vector<unsigned> tempNodeIdReduced;
	thrust::device_vector<double> tempNodeForceXReduced;
	thrust::device_vector<double> tempNodeForceYReduced;
	thrust::device_vector<double> tempNodeForceZReduced;
};

struct GeneralParams{
	double kT;
	double tau;
	double solve_time;
	unsigned iteration = 0;
	unsigned maxNodeCount;
	//parameters for advancing timestep and determining equilibrium
	
	double dt;
	double nodeMass = 1.0;

};


class Storage;
struct HostSetInfoVecs;

class System {
public:
	GeneralParams generalParams;
	DomainParams domainParams;
	AuxVecs auxVecs;
	CoordInfoVecs coordInfoVecs;

	LinearSpringInfoVecs linearSpringInfoVecs;
	BendingTriangleInfoVecs bendingTriangleInfoVecs;
	AreaTriangleInfoVecs areaTriangleInfoVecs;
	LJInfoVecs ljInfoVecs;

	std::shared_ptr<Storage> storage;
	
public:

	System();

	void PrintForce();

	void initializeSystem(HostSetInfoVecs& hostSetInfoVecs);
	
	void assignStorage(std::shared_ptr<Storage> _storage);

	void solveSystem();

	void setExtras();
};


#endif /*POLYMERSYSTEM_H_*/