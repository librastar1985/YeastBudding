

#ifndef SYSTEM_H_
#define SYSTEM_H_

#pragma once

//#include <gsl/gsl_matrix.h>
#include <fstream>
//#include <gsl/gsl_vector.h>
//#include <gsl/gsl_math.h>
//#include <gsl/gsl_multimin.h>
//#include <gsl/gsl_errno.h>

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

//
struct CapsidInfoVecs {
	thrust::device_vector<int> id_bucket;	//bucket id
	// bucket value means global rank of a certain point
	thrust::device_vector<int> id_value;//node id

 	thrust::device_vector<double> nodeLocX;
	thrust::device_vector<double> nodeLocY;
	thrust::device_vector<double> nodeLocZ;


	//used for capside linear springs binding to membrane
	thrust::device_vector<int> tempMembraneId;
	thrust::device_vector<int> tempCapsideId;

	thrust::device_vector<double> tempLengthsPairs;
	thrust::device_vector<double> tempNodeForceX;
	thrust::device_vector<double> tempNodeForceY;
	thrust::device_vector<double> tempNodeForceZ;

	//used for capside repulsion
	int factor = 10;//number of default nodes repulsion can interact with.
	thrust::device_vector<int> tempNodeIdUnreduced;
	thrust::device_vector<double> tempNodeForceXUnreduced;
	thrust::device_vector<double> tempNodeForceYUnreduced;
	thrust::device_vector<double> tempNodeForceZUnreduced;

	thrust::device_vector<int> tempNodeIdReduced;
	thrust::device_vector<double> tempNodeForceXReduced;
	thrust::device_vector<double> tempNodeForceYReduced;
	thrust::device_vector<double> tempNodeForceZReduced;

	double spring_constant = 10.0;
	double length_zero = 0.97;
	double length_cutoff = 1.63;
	int maxNodeCount;
	double viscosity = 1.0;
	double forceX = 0.0;
	double forceY = 0.0;
	double forceZ = 0.0;

	int num_connections=0;

};

//Data Structure for node location. velocity and force
struct CoordInfoVecs {

	thrust::device_vector<double> scaling_per_edge;

	thrust::device_vector<int> nodes2Triangles_1;
	thrust::device_vector<int> nodes2Triangles_2;
	thrust::device_vector<int> nodes2Triangles_3;
	thrust::device_vector<int> nodes2Triangles_4;
	thrust::device_vector<int> nodes2Triangles_5;
	thrust::device_vector<int> nodes2Triangles_6;
	thrust::device_vector<int> nodes2Triangles_7;
	thrust::device_vector<int> nodes2Triangles_8;
	thrust::device_vector<int> nodes2Triangles_9;
//	thrust::device_vector<int> nodes2Triangles_10;
//	thrust::device_vector<int> nodes2Triangles_11;
//	thrust::device_vector<int> nodes2Triangles_12;

	thrust::device_vector<double> SurfaceNormalX;
	thrust::device_vector<double> SurfaceNormalY;
	thrust::device_vector<double> SurfaceNormalZ;

	thrust::device_vector<int> nndata1;
	thrust::device_vector<int> nndata2;
	thrust::device_vector<int> nndata3;
	thrust::device_vector<int> nndata4;
	thrust::device_vector<int> nndata5;
	thrust::device_vector<int> nndata6;
	thrust::device_vector<int> nndata7;
	thrust::device_vector<int> nndata8;
	thrust::device_vector<int> nndata9;
//	thrust::device_vector<int> nndata10;
//	thrust::device_vector<int> nndata11;
//	thrust::device_vector<int> nndata12;

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

	thrust::device_vector<double> tempNodeForceXReduced;
	thrust::device_vector<double> tempNodeForceYReduced;
	thrust::device_vector<double> tempNodeForceZReduced;
	thrust::device_vector<double> tempNodeForceXUnreduced;
	thrust::device_vector<double> tempNodeForceYUnreduced;
	thrust::device_vector<double> tempNodeForceZUnreduced;

	//LOCAL COORDS
	//indices of each triangle
	int num_triangles;
	thrust::device_vector<int> triangles2Nodes_1;
	thrust::device_vector<int> triangles2Nodes_2;
	thrust::device_vector<int> triangles2Nodes_3;


	//indices of each edge
	int num_edges;
	thrust::device_vector<int> edges2Nodes_1;
	thrust::device_vector<int> edges2Nodes_2;

	//indices of 2 triangle on each edge
	thrust::device_vector<int> edges2Triangles_1;
	thrust::device_vector<int> edges2Triangles_2;

	//indices of edges on each triangle.
	thrust::device_vector<int> triangles2Edges_1;
	thrust::device_vector<int> triangles2Edges_2;
	thrust::device_vector<int> triangles2Edges_3;

};



//struct used for linking of nodes in network
struct AuxVecs {
	// bucket key means which bucket ID does a certain point fit into
	thrust::device_vector<int> id_bucket;	//bucket id
	// bucket value means global rank of a certain point
	thrust::device_vector<int> id_value;//node id
	// bucket key expanded means what are the bucket IDs are the neighbors of a certain point
	thrust::device_vector<int> id_bucket_expanded;
	// bucket value expanded means each point ( represented by its global rank) will have multiple copies
	thrust::device_vector<int> id_value_expanded;

	// begin position of a keys in id_bucket_expanded and id_value_expanded
	//entry keyBegin[bucketKey] returns start of indices to link
	thrust::device_vector<int> keyBegin;
	// end position of a keys in id_bucket_expanded and id_value_expanded
	thrust::device_vector<int> keyEnd;

	int endIndexid_bucket;
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
	double gridSpacing = 1.5;//bucket scheme search distance, must be larger than cutoff for capsid
	int XBucketCount;
	int YBucketCount;
	int ZBucketCount;
	int totalBucketCount = 0;
};

struct LJInfoVecs{
	double LJ_PosX;
	double LJ_PosY;
	double LJ_PosZ;
	thrust::device_vector<double> LJ_PosX_all;
	thrust::device_vector<double> LJ_PosY_all;
	thrust::device_vector<double> LJ_PosZ_all;



	double Rmin_M=0.97;//1.0;
	double Rcutoff_M=0.97;//1.0;
	double Rmin_LJ;
	double Rcutoff_LJ;

	double epsilon_M=1.0;
	double epsilon_M_att1;
	double epsilon_M_att2;
	double epsilon_M_rep1;
	double epsilon_M_rep2;
	double epsilon_LJ;
	double epsilon_LJ_rep1;
	double epsilon_LJ_rep2;
	double spring_constant;

	thrust::device_vector<int> node_id_close;
	double lj_energy_M;
	double lj_energy_LJ;
	double forceX;
	double forceY;
	double forceZ;
	thrust::device_vector<double> forceX_all;
	thrust::device_vector<double> forceY_all;
	thrust::device_vector<double> forceZ_all;

};

struct AreaTriangleInfoVecs {
	double dummy;
	int factor = 3;//used for reduction
	double initial_area = 0.0048013;
	double spring_constant;
	double spring_constant_weak;

	double area_triangle_energy;

	thrust::device_vector<int> tempNodeIdUnreduced;
	thrust::device_vector<double> tempNodeForceXUnreduced;
	thrust::device_vector<double> tempNodeForceYUnreduced;
	thrust::device_vector<double> tempNodeForceZUnreduced;

	thrust::device_vector<int> tempNodeIdReduced;
	thrust::device_vector<double> tempNodeForceXReduced;
	thrust::device_vector<double> tempNodeForceYReduced;
	thrust::device_vector<double> tempNodeForceZReduced;

};

struct BendingTriangleInfoVecs {
	int numBendingSprings=0;

	int factor = 4;//used for reduction
	double spring_constant;
	double spring_constant_weak;
	double spring_constant_raft;
	double spring_constant_coat;
	double initial_angle = 0.0;//radians
	double initial_angle_raft;
	double initial_angle_coat;
	double initial_angle_bud;

	double bending_triangle_energy;

	thrust::device_vector<int> tempNodeIdUnreduced;
	thrust::device_vector<double> tempNodeForceXUnreduced;
	thrust::device_vector<double> tempNodeForceYUnreduced;
	thrust::device_vector<double> tempNodeForceZUnreduced;

	thrust::device_vector<int> tempNodeIdReduced;
	thrust::device_vector<double> tempNodeForceXReduced;
	thrust::device_vector<double> tempNodeForceYReduced;
	thrust::device_vector<double> tempNodeForceZReduced;
};
struct LinearSpringInfoVecs {

	int factor = 2;//used for reduction
	double spring_constant;
	double spring_constant_weak;
	double spring_constant_att1;
	double spring_constant_att2;
	double spring_constant_rep1; //This is the "D" in Morse potential
	double spring_constant_rep2; //This is the "a" in Morse potential

	double linear_spring_energy;
	double memrepulsion_energy;
	double scalar_edge_length;
	
	thrust::device_vector<double> edge_initial_length;

	thrust::device_vector<int> tempNodeIdUnreduced;
	thrust::device_vector<double> tempNodeForceXUnreduced;
	thrust::device_vector<double> tempNodeForceYUnreduced;
	thrust::device_vector<double> tempNodeForceZUnreduced;

	thrust::device_vector<int> tempNodeIdReduced;
	thrust::device_vector<double> tempNodeForceXReduced;
	thrust::device_vector<double> tempNodeForceYReduced;
	thrust::device_vector<double> tempNodeForceZReduced;
};

struct GeneralParams{
	double kT;
	double kT_growth;
	double tau;
	int solve_time=100;
	double Rmin = 1.0; //Current mesh minimum edge length, subject to change.
	double Rmin_growth;
	double abs_Rmin;
	int iteration = 0;
	int maxNodeCount;
	int maxNodeCountLJ;
	double length_scale;
	//parameters for advancing timestep and determining equilibrium

	double dt;
	double nodeMass = 1.0;
	int edges_in_upperhem_list_length;
	double growth_energy_scaling;
	thrust::device_vector<int> edge_to_ljparticle;
	thrust::device_vector<int> nodes_in_upperhem;
	thrust::device_vector<int> edges_in_upperhem;
	thrust::device_vector<int> edges_in_upperhem_list;
	thrust::device_vector<int> triangles_in_upperhem;
	thrust::device_vector<int> boundaries_in_upperhem;
	thrust::device_vector<double> angle_per_edge;
	double centerX = 0.0;
	double centerY = 0.0;
	double centerZ = 0.0;

	double current_total_volume;
	double true_current_total_volume;
	double eq_total_volume;
	double volume_spring_constant;
	double volume_energy;
	double eq_total_boundary_length;
	double line_tension_energy;
	double line_tension_constant;
	double safeguardthreshold;
	

	//int num_of_triangles;
	//int num_of_edges;
	int true_num_edges;

	double insertion_energy_cost;
	double strain_threshold;

	int SCALE_TYPE;
	double scaling_pow;
	double gausssigma;
	double hilleqnconst;
	double hilleqnpow;

	thrust::device_vector<int> no_weakening;
	double septin_ring_z;
	double boundary_z;

};


class Storage;
class SystemBuilder;
struct HostSetInfoVecs;

class System {
public:
	std::weak_ptr<SystemBuilder> weak_bld_ptr;
	GeneralParams generalParams;
	DomainParams domainParams;
	AuxVecs auxVecs;
	CoordInfoVecs coordInfoVecs;

	CapsidInfoVecs capsidInfoVecs;
	LinearSpringInfoVecs linearSpringInfoVecs;
	BendingTriangleInfoVecs bendingTriangleInfoVecs;
	AreaTriangleInfoVecs areaTriangleInfoVecs;
	LJInfoVecs ljInfoVecs;

	std::shared_ptr<Storage> storage;

	//gsl_vector* df;
	//gsl_vector* locations;
    



public:

	System();

	void set_weak_builder(std::weak_ptr<SystemBuilder> _weak_bld_ptr);

	void PrintForce();

	void initializeSystem(HostSetInfoVecs& hostSetInfoVecs);

	void assignStorage(std::shared_ptr<Storage> _storage);
 
	void solveSystem();

	void setExtras();

	void setBucketScheme();

	//void Solve_Forces(const gsl_vector* temp_locations);
	void Solve_Forces();
	//double Solve_Energy(const gsl_vector* temp_locations);

	//void dev_to_gsl_loc_update(gsl_vector* temp_locations);
	//void gsl_to_dev_loc_update(const gsl_vector* temp_locations);
	

};


#endif /*POLYMERSYSTEM_H_*/
