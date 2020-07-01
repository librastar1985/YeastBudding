#include "System.h"
#include "SystemStructures.h" 
#include "AreaTriangles.h"
//#include "AreaTrianglesEnergy.h"
#include "BendingTriangles.h"
//#include "BendingTrianglesEnergy.h"
#include "MemRepulsionSprings_universal.h"
#include "MemRepulsionSprings_local.h"
#include "MemRepulsionEnergy.h"
#include "LinearSprings.h"
//#include "LinearSpringsEnergy.h"
//#include "LJSprings.h"
//#include "LJSprings_LJ.h"
#include "AreaCompBud.h"
#include "NodeAdvance.h"
#include "BucketScheme.h"
#include "Storage.h" 
#include "Edgeswap_test.h"
#include "SystemBuilder.h"
#include <vector>
#include "VolumeComp.h"
#include "VolumeSprings.h"
#include <bits/stdc++.h>
#include "LineTensionSprings.h"
//#include "Growth.h"
#include <math.h>
//#include "SurfaceNormal.h"
//#include "Nodes2Triangles.h"
#include "TurgorForce.h"

///////////////////////////////////////////////////////////////////
///////////////////////// WARNING ////////////////////////////////
//////////////////REMEMBER TO CHANGE THE /////////////////////////
/////////////EQUILIBRIUM LENGTH OF EACH TRIANGLE EDGE /////////////
//////////////// IN THE VECTOR INITIALIZATION ////////////////////
//////////////////SECTION TOWARD THE END OF THE CODE /////////////
////////////////////////////////////////////////////////////////////

 //somehow the gradient is not being set in my version

//bool IsPos (int i){return (i>=0);}
int count_bigger(const std::vector<int>& elems) {
    return std::count_if(elems.begin(), elems.end(), [](int c){return c >= 0;});
}

System::System() {};

void System::Solve_Forces(){

	thrust::fill(coordInfoVecs.nodeForceX.begin(), coordInfoVecs.nodeForceX.end(), 0.0);
	thrust::fill(coordInfoVecs.nodeForceY.begin(), coordInfoVecs.nodeForceY.end(), 0.0);
	thrust::fill(coordInfoVecs.nodeForceZ.begin(), coordInfoVecs.nodeForceZ.end(), 0.0);
	
	//setBucketScheme();
	ComputeLinearSprings(
		generalParams, 
		coordInfoVecs,
		linearSpringInfoVecs, 
		ljInfoVecs);
	
	// ComputeAreaTriangleSprings(
		
	// 	generalParams,
	// 	coordInfoVecs,
	// 	areaTriangleInfoVecs);

	ComputeTurgorSprings(
		generalParams,
		coordInfoVecs,
		areaTriangleInfoVecs
	);

	ComputeCosTriangleSprings(
		
		generalParams,
		coordInfoVecs,  
		bendingTriangleInfoVecs); 
	
	ComputeMemRepulsionSprings_universal(
		coordInfoVecs,
		linearSpringInfoVecs, 
		capsidInfoVecs,
		generalParams,
		auxVecs);

	// ComputeMemRepulsionSprings_local(
	// 	coordInfoVecs,
	// 	linearSpringInfoVecs, 
	// 	capsidInfoVecs,
	// 	generalParams,
	// 	auxVecs);

	// ComputeVolume(
	// 	generalParams,
	// 	coordInfoVecs,
	// 	linearSpringInfoVecs,
	// 	ljInfoVecs);


	/*ComputeVolumeSprings(
		coordInfoVecs,
		linearSpringInfoVecs, 
		capsidInfoVecs,
		generalParams,
		auxVecs);*/

	/* if (generalParams.true_current_total_volume/initial_volume >= 1.25){
	ComputeLineTensionSprings(
		generalParams,
		coordInfoVecs,
		linearSpringInfoVecs);
	} */
		
};


void System::solveSystem() {

	bool bud_tip_necessary = false;
	
	double current_center_x = 0.0;
	double current_center_y = 0.0;
	std::vector<int> septin_node = {76,	79,	111, 113,151,153,505,506,515,516,593,632};
	int GROWTH_COUNTER = 0;
	int min_num_edge_loop = 1;
	std::cout<<"min_num_edge_loop for edgeswap = "<<min_num_edge_loop<<std::endl;

	std::random_device rand_dev;
	std::mt19937 generator2(rand_dev());
	std::mt19937 generator_edgeswap(rand_dev());

	double MAX_VOLUME_RATIO = 2.0;
	double MAX_VOLUME_RATIO_stiffening = 1.5;
	// int MAX_GROWTH_NUMBER = 1;
	// std::cout<<"MAX_GROWTH_NUMBER (# of edge to expand) = "<<MAX_GROWTH_NUMBER<<std::endl;
	int GROWTH_FREQUENCY = 15;//100;
	std::cout<<"GROWTH_FREQ (how many times Max_Runtime has to be reached to perform growth"<<GROWTH_FREQUENCY<<std::endl;
	int MAX_GROWTH_TEST = 1;//100000;
	std::cout<<"MAX growth event (not attempt) per growth algorithm triggering = "<<MAX_GROWTH_TEST<<std::endl;
	double energy_gradient_threshold = 0.0;
	std::cout<<"ENERGY_GRADIENT_THRESHOLD = "<<energy_gradient_threshold<<std::endl;

	generalParams.kT_growth = 1.0;
	generalParams.SCALE_TYPE = 1;
	// 0:= Gaussian-like weakening
	// 1:= a1*(pow(x,b)) + a2*(1-pow(x,b)) type weakening
	// 2:= pure Gaussian weakening
	// 3:= isotropic
	// 4:= hill equation
	//Note that (3) is used in combination with sigma = INT_MAX;
	std::cout<<"SCALE TYPE = "<<generalParams.SCALE_TYPE<<std::endl;
	std::cout<<"0:= sigmoidal Gaussian-like weakening, 1:= a1*(pow(x,b)) + a2*(1-pow(x,b)) type weakening, 2:= pure Gaussian weakening, 3:= isotropic, 4:= hill equation"<<std::endl;
	generalParams.scaling_pow = 1.0;
	std::cout<<"scaling_pow (this is for SCALE_TYPE = 1 case) = "<<generalParams.scaling_pow<<std::endl;
	generalParams.gausssigma = 0.1;
	std::cout<<"gausssigma (this is for the SCALE_TYPE = 0 case) = "<<generalParams.gausssigma<<std::endl;
	//coordInfoVecs.scaling_per_edge.
	//generalParams.hilleqnconst = 0.9;
	//generalParams.hilleqnpow = 40.0;
	std::vector<int> nodes_in_growth;
	std::vector<int> triangles_in_growth;
	std::vector<int> edges_in_growth;
	double dtb; //dtb := distance to boundary
	double dtb_max; //dtb_max := the max distance used to calculate the distance ratio in the Hill equation.
	double sigma = 0.0;//INT_MAX; //if this is set to be INT_MAX then we assume isotropic weakening.
	double sigma_true = sqrt(0.5); //This is the variance used to calculate the scaling of the wall weakening.
	std::cout<<"initial sigma (for gradient distribution variance), based on initial distribution of Cdc42, if using true gaussian weakening = "<<sigma<<std::endl;
	std::cout<<"If sigma = INT_MAX, then we have isotropic weakening scenario"<<std::endl;
	std::cout<<"true sigma (for gaussian-related distribution variance) = "<<sigma_true<<std::endl;

	generalParams.insertion_energy_cost = -log(0.0025);
	std::cout<<"GROWTH: material insertion energy cost (dependent on local chemical concentration) = "<<generalParams.insertion_energy_cost<<std::endl;
	generalParams.strain_threshold = 0.05;// 0.05;//0.01;
	std::cout<<"GROWTH: critical strain threshold used for insertion probability calculation = "<<generalParams.strain_threshold<<std::endl;

	generalParams.growth_energy_scaling = 1.0;//0.01375;
	std::cout<<"GROWTH ENERGY SCALING FOR MATERIAL INSERTION PROBABILITY = "<<generalParams.growth_energy_scaling<<std::endl;
	generalParams.safeguardthreshold = 9;
	std::cout<<"NEIGHBOR SAFE GUARD THRESHOLD = "<<generalParams.safeguardthreshold<<std::endl;
	//safeguardthreshold is the maximum number of neighboring nodes a node can have.

	//////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////
	////////////////////////// PARAMETER SETTINGS ////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////

	generalParams.nodeMass = 400.0;//100.0;//
	std::cout<<"Damping coefficient for node movement = "<<generalParams.nodeMass<<std::endl;
	double Max_Runtime = generalParams.dt*50.0;
	double Max_RunStep = Max_Runtime/generalParams.dt;
	std::cout<<"Max runtime = "<<Max_Runtime<<std::endl;
	std::cout<<"Max runstep = "<<Max_RunStep<<std::endl;
	bool runSim = true;
	int num_edge_loop;
	double initial_kT;
	initial_kT = generalParams.kT;//This is for the acceptance of change after looping through every edge within proximity.
	double SAMPLE_SIZE = 0.05;
	std::cout<<"Sample ratio: "<<SAMPLE_SIZE<<std::endl;
	std::cout<<"If the Sample raio is 0, it means we have chosen a fixed number of attempt throughout the simulation"<<std::endl;
	//This determines the number of edges to test for bondflip remeshing

	auto edgeswap_ptr = std::make_shared<Edgeswap>(coordInfoVecs, generalParams);
	int RECORD_TIME = 1;//round(Max_RunStep/2);
	std::cout<<"Record frequency = "<<RECORD_TIME<<std::endl;
	//int GROWTH_TIME = 1;
	//std::cout<<"Growth frequency = "<<GROWTH_TIME<<std::endl;
	int translate_frequency = 200;
	//std::cout<<"translate + edgeswap frequency = "<<translate_frequency<<std::endl;
	//translate_frequency determines the frequency for the mesh to re-center and perform dynamical remeshing
	//int MAX_GROWTH_EVENT_NUMBER = 2000;//800;
	int NKBT = 2*8200;//10000;//7500; 
	std::cout<<"Number of edge-swap per kBT value (or total number of edge-swap event if kBT is fixed) = "<<NKBT<<std::endl;
	std::cout<<"Maximal iterative steps = "<<NKBT*Max_Runtime<<std::endl;
	int GROWTH_FREQUENCY_SCALE = 100;//4;
	std::cout<<"GROWTH FREQ SCALE: decides how many growth event must be checked before recording the result"<<std::endl;
	double min_kT = -0.1;//0.21;
	std::cout<<"min kT for simulation termination = "<<min_kT<<std::endl;
	int WHEN = 0;
	double old_total_energy = 0.0;
	double new_total_energy = 0.0;
	double energy_gradient = 0.0;
	double energy_rep = 0.0;
	int Num_of_step_run = 0;
	auto build_ptr = weak_bld_ptr.lock();//upgrade weak builder to access host variables.
	//std::cout<<"initial LJ-x : "<< ljInfoVecs.LJ_PosX <<std::endl;
	//std::cout<<"initial LJ-y : "<< ljInfoVecs.LJ_PosY <<std::endl;
	//std::cout<<"initial LJ-z : "<< ljInfoVecs.LJ_PosZ <<std::endl;
		

    
	double min_energy;
	generalParams.true_num_edges = 0;
	for (int i = 0; i < coordInfoVecs.num_edges; i++){
		if (coordInfoVecs.edges2Nodes_1[i] != INT_MAX && coordInfoVecs.edges2Nodes_2[i] != INT_MAX){
			generalParams.true_num_edges += 1;
		}
	}
	
	//double COMPRESS = 2.0227;
	// double COMPRESS2 = -2.0227;

	/////////////////////////////////////////////////////////////////
	/////////////////////// MEMBRANE RELATED ////////////////////////
	/////////////////////////////////////////////////////////////////
	
	std::vector<double> nodenormal_1(generalParams.maxNodeCount, 0.0);
	std::vector<double> nodenormal_2(generalParams.maxNodeCount, 0.0);
	std::vector<double> nodenormal_3(generalParams.maxNodeCount, 0.0);
	int reduce_counter = 0;

	double VOLUME_FACTOR = 1.0;//2.25;
	//VOLUME_FACTOR determines the target volume which equals to VOLUME_FACTOR*initial_volume.
	//double tip_depth = 0.5;
	//tip_depth is currently unused.

	double LINE_TENSION_THRESHOLD = 0.0;
	std::cout<<"LINE TENSION THRESHOLD for activation of line tension = "<<LINE_TENSION_THRESHOLD<<std::endl;
	double VOLUME_THRESHOLD = 0.0;
	std::cout<<"VOLUME THRESHOLD for activation of weakened membrane = "<<VOLUME_THRESHOLD<<std::endl;
	
	double weakened = 1.90;//6.0;
	//weakened determines the minimum height of the z-coordinate of the membrane node to be considered in the area of weakened mechanical properties.
	//double tip_base = 6.0;
	//tip_base currently unused.

	double EXPAN_THRESHOLD = 0.0;
	double EXPAN_THRESHOLD_weak = 0.0;//1.75;
	std::cout<<"EXPANSION THRESHOLD = "<<EXPAN_THRESHOLD<<std::endl;
	int RULES_OF_EXPAN = 1;	//EXPAN_THRESHOLD is the yielding ratio where a pair of triangles will perform expansion.
	
	std::cout<<"EXPANSION RULE = "<<RULES_OF_EXPAN<<std::endl;
	//EXPAN_THRESHOLD_weak is the secondary yielding ratio.
	//RULES_OF_EXPAN controls how the EXPAN_THRESHOLD is applied:
	// 1:= Both trianglular areas must exceed the threshold value.
	// 2:= If one trianglular area exceeds the treshold value while the other exceeds the secondary threshold value.
	// 3:= If the combined area of the two triangles exceed 2*EXPAN_THRESHOLD.
	// 4:= If a selected edges exceed the threshold value, split the triangles associated with the edge.

	for (int i = 0; i < generalParams.maxNodeCount; i++){
		generalParams.centerX += coordInfoVecs.nodeLocX[i];
		generalParams.centerY += coordInfoVecs.nodeLocY[i];
		generalParams.centerZ += coordInfoVecs.nodeLocZ[i];
	}
	generalParams.centerX = generalParams.centerX/generalParams.maxNodeCount;
	generalParams.centerY = generalParams.centerY/generalParams.maxNodeCount;
	generalParams.centerZ = generalParams.centerZ/generalParams.maxNodeCount;
	double displacementX, displacementY, displacementZ;
	double newcenterX, newcenterY, newcenterZ;
	//centerX, centerY, centerZ is determined as the referenced origin for recentering of the mesh.

	std::vector<int> VectorShuffleForGrowthLoop;
	std::vector<int> VectorShuffleForFilamentLoop;
	std::vector<int> VectorShuffleForEdgeswapLoop;

	double max_height = coordInfoVecs.nodeLocZ[35];
	double min_height = coordInfoVecs.nodeLocZ[38];
	int max_height_index = 35;
	/*double max_height = -10000.0;
	int max_height_index = -1;
	std::vector<int> Stiffness_gradient();
    for (int k = 0; k < generalParams.maxNodeCount; k++){
        if (coordInfoVecs. nodeLocZ[k] >= max_height){
			max_height = coordInfoVecs. nodeLocZ[k];
			max_height_index = k;
            }
	}*/
	//Max and min height of the membrane nodes, these have to be changed if the mesh used is changed.

	generalParams.Rmin = 0.301;//0.15;
	double distance_uniform_weak = 3.33*generalParams.Rmin;
	std::cout<<"distance from tip that undergoes uniform weakening = "<<distance_uniform_weak<<std::endl;
	//Equilibrium length of an edge of the triangle.
	//generalParams.Rmin_growth = 0.329;
	generalParams.abs_Rmin = generalParams.Rmin;//0.15;
	std::cout<<"abs_Rmin = "<<generalParams.abs_Rmin<<std::endl;
	//Equilibrium distance between membrane node for volume exclusion.
	areaTriangleInfoVecs.initial_area = 0.039;//0.03927344;//0.009817;
	std::cout<<"equilibrium triangular area = "<<areaTriangleInfoVecs.initial_area<<std::endl;
	//Equilibrium triangular area.
	ljInfoVecs.Rmin_M = 0.0;
	//Equilibrium distance between the nucleus particle and membrane.
	ljInfoVecs.Rcutoff_M = 0.0;
	//Maximal interaction range between the nucleus and membrane.
	ljInfoVecs.Rmin_LJ = 0.0;//3.0//1.0;
	//Equilibrium distance between nuclei.
	ljInfoVecs.Rcutoff_LJ = 0.0;//3.0;//1.0;
	//Maximal interaction range between the nuclei.
	ljInfoVecs.epsilon_M_att1 = 0.0;//6.0;//16.0;
	ljInfoVecs.epsilon_M_att2 = 0.0;//1.0;//1.0;
	std::cout<<"Morse_NM_D_att = "<<ljInfoVecs.epsilon_M_att1<<std::endl;
	std::cout<<"Morse_NM_a_att = "<<ljInfoVecs.epsilon_M_att2<<std::endl;
	//Coefficient for the attractive interaction between nuclei and membrane.
	ljInfoVecs.epsilon_M_rep1 = 0.0;//12.5;//16.0;
	ljInfoVecs.epsilon_M_rep2 = 0.0;//0.5;//1.0;
	std::cout<<"Morse_NM_D_rep = "<<ljInfoVecs.epsilon_M_rep1<<std::endl;
	std::cout<<"Morse_NM_a_rep = "<<ljInfoVecs.epsilon_M_rep2<<std::endl;
	//Coefficient for the repulsive interaction between nuclei and membrane.
	
	ljInfoVecs.epsilon_LJ_rep1 = 0.0;//10.0;//0.5;// 0.06;//7.5;
	ljInfoVecs.epsilon_LJ_rep2 = 0.0;//0.5;//1.0;//1.0;//1.0;
	std::cout<<"Morse_NN_D = "<<ljInfoVecs.epsilon_LJ_rep1<<std::endl;
	std::cout<<"Morse_NN_a = "<<ljInfoVecs.epsilon_LJ_rep2<<std::endl;
	//Coefficient of the interaction between nuclei.

	linearSpringInfoVecs.spring_constant_rep1 = 0.01;//0.023;
	linearSpringInfoVecs.spring_constant_rep2 = 9.0;//5.0;
	std::cout<<"Membrane volume exclusion Morse D = "<<linearSpringInfoVecs.spring_constant_rep1<<std::endl;
	std::cout<<"Membrane volume exclusion Morse a = "<<linearSpringInfoVecs.spring_constant_rep2<<std::endl;
	//The coefficient used for non-neighboring membrane node volume exclusion.
	//rep1 is the "D" and rep2 is the "alpha" in the standard form of Morse potential.

	generalParams.volume_spring_constant = 0.2;//(1.0/3.0)*areaTriangleInfoVecs.initial_area*1.0;
	std::cout<<"spring constant for surface normal expansion (pressure within the cell) = "<<generalParams.volume_spring_constant<<std::endl;
	generalParams.line_tension_constant = linearSpringInfoVecs.spring_constant;//75.0;//250.0;
	std::cout<<"spring constant for the septin ring = "<<generalParams.line_tension_constant<<std::endl;
	generalParams.length_scale = 1.0;//0.85;//0.1577;//1.0*generalParams.Rmin;// 0.8333;
	//std::cout<<"equilibrium length of each segment of the septin ring = "<<generalParams.length_scale<<std::endl;

	double scale_linear = linearSpringInfoVecs.spring_constant*1.0;//0.144;//0.25;//25.0/2.5;//75.0/15.0;
	double scale_bend = bendingTriangleInfoVecs.spring_constant*0.05;//0.052;//0.05;//10.0/1.0;//75.0/7.5;
	double scale_area = areaTriangleInfoVecs.spring_constant*0.0;//0.1;//0.25;//50.0/5.0;//75.0/15.0;
	std::cout<<"weakened region linear = "<<scale_linear<<std::endl;
	std::cout<<"weakened region bend = "<<scale_bend<<std::endl;
	std::cout<<"weakened region area = "<<scale_area<<std::endl;
	//linearSpringInfoVecs.spring_constant_weak = linearSpringInfoVecs.spring_constant/scale_linear;
	//bendingTriangleInfoVecs.spring_constant_weak = bendingTriangleInfoVecs.spring_constant/scale_bend;
	//areaTriangleInfoVecs.spring_constant_weak = areaTriangleInfoVecs.spring_constant/scale_area;
	
	// linearSpringInfoVecs.spring_constant_weak = scale_linear;
	// bendingTriangleInfoVecs.spring_constant_weak = scale_bend;
	// areaTriangleInfoVecs.spring_constant_weak = scale_area;
	//Scaling of the weakend mechanical properties.
	
	bendingTriangleInfoVecs.initial_angle = 0.0;//0.087249;//0.04335;
	//bendingTriangleInfoVecs.initial_angle_bud = 0.087249;
	bendingTriangleInfoVecs.initial_angle_raft = 0.0;//0.087249;//0.04335;
	bendingTriangleInfoVecs.initial_angle_coat = 0.0;//0.087249;//0.04335;
	std::cout<<"equilibrium bending angle of the membrane = "<<bendingTriangleInfoVecs.initial_angle<<std::endl;
	//raft and coat are current unused due to the assumption of uniform preferred curvature.
	
	bendingTriangleInfoVecs.spring_constant_raft = 0.0;//bendingTriangleInfoVecs.spring_constant;
	bendingTriangleInfoVecs.spring_constant_coat = 0.0;//bendingTriangleInfoVecs.spring_constant;
	bendingTriangleInfoVecs.spring_constant = bendingTriangleInfoVecs.spring_constant*(2.0/sqrt(3));
	bendingTriangleInfoVecs.spring_constant_raft = bendingTriangleInfoVecs.spring_constant_raft*(2.0/sqrt(3));
	bendingTriangleInfoVecs.spring_constant_coat = bendingTriangleInfoVecs.spring_constant_coat*(2.0/sqrt(3));
	std::cout<<"Effective bending coefficient is calculated by multiplying 2/sqrt(3)"<<std::endl;
	std::cout<<"effective bending coefficient of the membrane = "<<bendingTriangleInfoVecs.spring_constant<<std::endl;
	std::cout<<"effective bending coefficient of the membrane raft = "<<bendingTriangleInfoVecs.spring_constant_raft<<std::endl;
	std::cout<<"effective bending coefficient of the membrane coat = "<<bendingTriangleInfoVecs.spring_constant_coat<<std::endl;

	std::vector<int> pull_nodes_up;// = {35,    76,    79,   111,   113,   151,   153,   360,   361,   362,   363,   364,   365,   505,   506,   515,   516,   593,   632};//{35, 360,   361,   362,   363,   364,   365};
	std::vector<int> pull_nodes_down;// = {38,    86,    89,   121,   123,   144,   146,   378,   379,   380,   381,   382,   383,   535,   536,   545,   546,   602,   626};//{38, 378,   379,   380,   381,   382,   383};
	std::vector<int> push_nodes_down;
	std::vector<int> push_nodes_up;
	for (int i = 0; i < generalParams.maxNodeCount; i++){
		if (coordInfoVecs.nodeLocZ[i] >= 1.43026488631){
			pull_nodes_up.push_back(i);
		}
		if (coordInfoVecs.nodeLocZ[i] <= -1.43026488631){
			pull_nodes_down.push_back(i);
		}
	}

	/////////////////////////////////////////////////////////////////
	////////////////// END OF MEMBRANE RELATED //////////////////////
	/////////////////////////////////////////////////////////////////


	/////////////////////////////////////////////////////////////////
	//////////////////////// NULCEUS RELATED ////////////////////////
	/////////////////////////////////////////////////////////////////
	double beta1 = 0.0;
	double beta2 = 0.0;
	std::cout<<"manual push speed for the nucleus tip = "<<beta1<<std::endl;
	std::cout<<"manual push speed for the remainder of the nucleus = "<<beta2<<std::endl;
	//beta1 is the vertical speed (0, 0, beta1) applied to the nucleus tip.
	//beta2 is the vertical speed (0, 0, beta2) applied to the remainder of the nucleus.

	std::vector<double> V1 = {-0.0};/*, 0.0  ,  0.1966  ,  0.5547 ,  -0.4689 ,   0.2422 ,  -0.2229,
							   -0.4312 ,  -0.0185 ,   0.2887 ,   0.3187 ,   0.7140 ,  
								0.2231 ,  -0.1921 ,	  -0.5541 ,   -0.1542 ,   -0.1689 ,    0.4391 ,
							   -0.6661 ,  -0.6381 ,   0.6256 ,   0.0466 ,  -0.0610 ,   0.5134};
								*/
	std::vector<double> V2 = {0.0};/*, 0.0 ,  -0.4595 ,  -0.4129 ,   0.0954 ,   0.1764 ,   0.4186 ,
							  -0.5602 ,  -0.6082 ,  -0.5318 ,   0.3561 ,   0.0753 ,
							  -0.0917 ,  -0.2596 , 0.2871 ,  -0.3918 ,   0.5195 ,   0.5579 ,
							  -0.2805 ,   0.0133  , -0.0073 ,   0.7426 ,   0.0614 ,  -0.1506};
								*/
	std::vector<double> V3 = { 0.6390};/*, 0.0 ,  -0.5511 ,   0.0267 ,  -0.5240  , -0.4004 ,   0.2850 ,
							   0.2032 ,  -0.1771 ,   0.4048 ,   0.3461 ,  -0.2034 ,
							   0.5041 ,  -0.4535 ,	-0.1241 ,   0.5722 ,  -0.3748 ,  -0.1335 ,
							   -0.0851 ,   0.3213 ,   0.2389 ,   0.0044 ,  -0.7424 ,  -0.7450};
							   */
	//V1, V2, and V3 are the (x,y,z)-coordinate of the nucleus particles.

	for (int i = 0; i < V1.size(); i++){
		ljInfoVecs.LJ_PosX_all.push_back(V1[i]); 
		ljInfoVecs.LJ_PosY_all.push_back(V2[i]);
		ljInfoVecs.LJ_PosZ_all.push_back(V3[i]);
	}  
	
	double NUCLEUS_UPPERHEM_BASE = 0.5;
	double NUCLEUS_LOWERHEM_BASE = -0.6;
	//These values defines the z-coordinate requirement for nucleus particles to be considered tip-region or base-region. This is used to 
	// determine where to apply spring or constant force.

	//////////////////////////////////////////////////////////////////
	///////////////// END OF NUCLEUS RELATED /////////////////////////
	//////////////////////////////////////////////////////////////////

	/*std::vector<int> filament_base(generalParams.maxNodeCountLJ, -1); //= {0,1,2,3,4,5,6,7,8,9,10,11};//{35, 21, 38, etc if we need more points}
	double filament_strength = 0.0;
	double filament_strength_pull = 1.0*filament_strength;
	double filament_Rmin = ((max_height - min_height)/4.0);
	std::cout<<"filament_strength = "<<filament_strength<<std::endl;
	std::cout<<"filament_strength for vertical pull = "<<filament_strength_pull<<std::endl;
	std::cout<<"filament_Rmin = "<<filament_Rmin<<std::endl;
	
	//First, determine the initial membrane nodes having filament bridges
	//with the nuclei particles
	for (int i = 0; i < generalParams.maxNodeCountLJ; i++){
		if (i == 0){
			filament_base[i] = 35;
			continue;
		}
		for (int j = 0; j < generalParams.maxNodeCount; j++){
			double xsquared = (ljInfoVecs.LJ_PosX_all[i] - coordInfoVecs.nodeLocX[j])*
								(ljInfoVecs.LJ_PosX_all[i] - coordInfoVecs.nodeLocX[j]);
			double ysquared = (ljInfoVecs.LJ_PosY_all[i] - coordInfoVecs.nodeLocY[j])*
								(ljInfoVecs.LJ_PosY_all[i] - coordInfoVecs.nodeLocY[j]);
			double zsquared = (ljInfoVecs.LJ_PosZ_all[i] - coordInfoVecs.nodeLocZ[j])*
								(ljInfoVecs.LJ_PosZ_all[i] - coordInfoVecs.nodeLocZ[j]);
			double R = sqrt(xsquared + ysquared + zsquared);
			if (R < filament_Rmin*1.1 && j != 35){
				filament_base[i] = j;
				break;
			}
		}
	}*/
	
	//std::vector<double> filament_Rmin;
	//for (int i = 0; i < V3.size();i++){
	//	filament_Rmin.push_back(sqrt((V3[i] - coordInfoVecs.nodeLocZ[38])*(V3[i] - coordInfoVecs.nodeLocZ[38])));
	//}
	//double filament_Rmin = sqrt((V3.back() - coordInfoVecs.nodeLocZ[38])*(V3.back() - coordInfoVecs.nodeLocZ[38]));
	//This part calculates the filament connecting the minimum point (in terms of z-coordinate) to the base of the nuclei cluster.


	//////////////////////////////////////////////////////////////////
	/////////// IDENTIFYING REGIONS WITH DIFFERENT MECH PROP /////////
	//////////////////////////////////////////////////////////////////

	/*ljInfoVecs.forceX_all.reserve(ljInfoVecs.LJ_PosX_all.size());
	ljInfoVecs.forceY_all.reserve(ljInfoVecs.LJ_PosX_all.size());
	ljInfoVecs.forceZ_all.reserve(ljInfoVecs.LJ_PosX_all.size());

	generalParams.maxNodeCountLJ = ljInfoVecs.LJ_PosX_all.size();
	std::vector<int> nucleus_in_upperhem(generalParams.maxNodeCountLJ, -1);
	std::vector<int> nucleus_in_lowerhem(generalParams.maxNodeCountLJ, -1);
	for (int i = 0; i < generalParams.maxNodeCountLJ; i++){
		if (ljInfoVecs.LJ_PosZ_all[i] > NUCLEUS_UPPERHEM_BASE){
			nucleus_in_upperhem[i] = 1;
		}
		if (ljInfoVecs.LJ_PosZ_all[i] < NUCLEUS_LOWERHEM_BASE){
			nucleus_in_lowerhem[i] = 1;
		}
	}*/
	

	std::vector<int> out;
	//int ALPHA;

	std::vector<bool> boundary_edges;
	boundary_edges.reserve(coordInfoVecs.num_edges);
	for (int i = 0; i < coordInfoVecs.num_edges; i++){
		if (coordInfoVecs.edges2Triangles_1[i] == coordInfoVecs.edges2Triangles_2[i]){
			boundary_edges.push_back(true);
		}
		else {
			boundary_edges.push_back(false);
		}
	}

	std::vector<int> edgeIndices;
	edgeIndices.reserve(coordInfoVecs.num_edges);
	for (int i = 0; i < coordInfoVecs.num_edges; ++i){
		//edgeIndices.push_back(edge_to_ljparticle[i]);
		if (boundary_edges[i] == false){
			edgeIndices.push_back(i);
		}
		else {
			edgeIndices.push_back(-1);
		}
	}



	auto it = remove_if(edgeIndices.begin(), edgeIndices.end(),  [](const int i) {return i < 0; });
	edgeIndices.erase(it, edgeIndices.end());
	
	std::vector<int> row2 = {35 ,   76 ,   79 ,  111 ,  113 ,  151 ,  153 ,  360 ,  361 ,  362 ,  363 ,  364 ,  365 ,  505 ,  506 ,  515 ,  516 ,  593 ,  632};
	//std::vector<int> nodes_to_center;
	//generalParams.nodes_in_upperhem.resize(generalParams.maxNodeCount,-1);

	for (int i = 0; i < generalParams.maxNodeCount; i++){
		generalParams.nodes_in_upperhem[i] = -1;
		generalParams.nodes_in_tip[i] = -1;
	}

	for (int i = 0; i < row2.size(); i++){
		generalParams.nodes_in_upperhem[row2[i]] = 1;
		generalParams.nodes_in_tip[row2[i]] = 1;
	//	std::cout<<"nodes "<<i<<" "<<generalParams.nodes_in_upperhem[i]<<std::endl;		
	}
	// for (int i = 0; i < generalParams.maxNodeCount; i++){
	// 	if (coordInfoVecs.nodeLocZ[i] > (generalParams.centerZ + weakened)){
	// 		generalParams.nodes_in_upperhem[i] = 1;
	// 	}
	// 	else{
	// 		generalParams.nodes_in_upperhem[i] = -1;
	// 	}
	// //	std::cout<<"nodes "<<i<<" "<<generalParams.nodes_in_upperhem[i]<<std::endl;		
	// }

	//std::vector<int> nodes_to_center;
	//std::vector<int> nodes_in_tip;
	for (int i = 0; i < coordInfoVecs.num_triangles; i++){
		int aaa = generalParams.nodes_in_tip[coordInfoVecs.triangles2Nodes_1[i]];
		//std::cout<<aaa<<std::endl;
		int bbb = generalParams.nodes_in_tip[coordInfoVecs.triangles2Nodes_2[i]];
		//std::cout<<bbb<<std::endl;
		int ccc = generalParams.nodes_in_tip[coordInfoVecs.triangles2Nodes_3[i]];
		//std::cout<<ccc<<std::endl;
		if ((aaa+bbb+ccc) >= 2){
			generalParams.triangles_in_tip[i] = 1;
			//triangles_in_upperhem.push_back(i);
		}
		//else if ((aaa+bbb+ccc)==1){
		//	generalParams.triangles_in_upperhem[i] = 0;
			//triangles_in_upperhem.push_back(i);
		//}
		else{
			generalParams.triangles_in_tip[i] = -1;
		}
	//	std::cout<<"triangle "<<i<<" "<<generalParams.triangles_in_upperhem[i]<<std::endl;		
	}

	//std::vector<int> edges_in_upperhem;
//	generalParams.edges_in_upperhem.resize(coordInfoVecs.num_edges);
	for (int i = 0; i < coordInfoVecs.num_edges; i++){
		int aaa = generalParams.triangles_in_tip[coordInfoVecs.edges2Triangles_1[i]];//generalParams.nodes_in_upperhem[coordInfoVecs.edges2Nodes_1[i]];
		int bbb = generalParams.triangles_in_tip[coordInfoVecs.edges2Triangles_2[i]];//generalParams.nodes_in_upperhem[coordInfoVecs.edges2Nodes_2[i]];
		if (aaa == 1 && bbb == 1){
			generalParams.edges_in_tip[i] = 1;
			//generalParams.edges_in_upperhem_list.push_back(i);
		}
		else if (aaa == 1 || bbb == 1){
			generalParams.edges_in_tip[i] = 1;
		}
		else{
			generalParams.edges_in_tip[i] = -1;
		}
		
	}

	//generalParams.triangles_in_upperhem.resize(coordInfoVecs.num_triangles);
	for (int i = 0; i < coordInfoVecs.num_triangles; i++){
		int aaa = generalParams.nodes_in_upperhem[coordInfoVecs.triangles2Nodes_1[i]];
		//std::cout<<aaa<<std::endl;
		int bbb = generalParams.nodes_in_upperhem[coordInfoVecs.triangles2Nodes_2[i]];
		//std::cout<<bbb<<std::endl;
		int ccc = generalParams.nodes_in_upperhem[coordInfoVecs.triangles2Nodes_3[i]];
		//std::cout<<ccc<<std::endl;
		if ((aaa+bbb+ccc)==3){
			generalParams.triangles_in_upperhem[i] = 1;
			//triangles_in_upperhem.push_back(i);
		}
		//else if ((aaa+bbb+ccc)==1){
		//	generalParams.triangles_in_upperhem[i] = 0;
			//triangles_in_upperhem.push_back(i);
		//}
		else{
			generalParams.triangles_in_upperhem[i] = -1;
		}
	//	std::cout<<"triangle "<<i<<" "<<generalParams.triangles_in_upperhem[i]<<std::endl;		
	}

	//std::vector<int> edges_in_upperhem;
//	generalParams.edges_in_upperhem.resize(coordInfoVecs.num_edges);
	int edges_in_upperhem_COUNT = 0;
	for (int i = 0; i < coordInfoVecs.num_edges; i++){
		int aaa = generalParams.triangles_in_upperhem[coordInfoVecs.edges2Triangles_1[i]];//generalParams.nodes_in_upperhem[coordInfoVecs.edges2Nodes_1[i]];
		int bbb = generalParams.triangles_in_upperhem[coordInfoVecs.edges2Triangles_2[i]];//generalParams.nodes_in_upperhem[coordInfoVecs.edges2Nodes_2[i]];
		if (aaa == 1 && bbb == 1){
			generalParams.edges_in_upperhem[i] = 1;
			//generalParams.edges_in_upperhem_list.push_back(i);
			generalParams.edges_in_upperhem_list[i] = i;
			edges_in_upperhem_COUNT += 1;
		}
		else if (aaa == 1 || bbb == 1){
			generalParams.edges_in_upperhem[i] = 1;
			generalParams.edges_in_upperhem_list[i] = -INT_MAX;
			edges_in_upperhem_COUNT += 1;
		}
		else{
			generalParams.edges_in_upperhem[i] = -1;
			generalParams.edges_in_upperhem_list[i] = -INT_MAX;
		}
		
	}
	std::cout<<"INITIAL EDGES IN UPPERHEM = "<<edges_in_upperhem_COUNT<<std::endl;

	int COUNTING_EDGE = 0;
	for (int y = 0; y < coordInfoVecs.num_edges; y++){
		if (generalParams.edges_in_upperhem_list[y] >= 0){
			COUNTING_EDGE += 1;
		}
		generalParams.edges_in_upperhem_list_length = COUNTING_EDGE;
	}
	

	//Find the boundary of the nodes_in_upperhem region
	//generalParams.boundaries_in_upperhem.resize(coordInfoVecs.num_edges);
	std::vector<int> boundary_node_list;
	std::vector<int> boundary_edge_list;
	for (int i = 0; i < coordInfoVecs.num_edges; i++){
		double T1 = coordInfoVecs.edges2Triangles_1[i];
		double T2 = coordInfoVecs.edges2Triangles_2[i];
		if (generalParams.triangles_in_upperhem[T1] == 1 && generalParams.triangles_in_upperhem[T2] != 1){
			generalParams.boundaries_in_upperhem[i] = 1;
			//std::cout<<generalParams.boundaries_in_upperhem[i]<<std::endl;
		//	generalParams.triangles_in_upperhem[T1] = 0;
		//	generalParams.triangles_in_upperhem[T2] = 0;
			double bdry_node1 = coordInfoVecs.edges2Nodes_1[i];
			double bdry_node2 = coordInfoVecs.edges2Nodes_2[i];
			boundary_node_list.push_back(bdry_node1);
			boundary_node_list.push_back(bdry_node2);
			boundary_edge_list.push_back(i);
			//generalParams.nodes_in_upperhem[bdry_node1] = 0;
			//generalParams.nodes_in_upperhem[bdry_node2] = 0;
			//coordInfoVecs.isNodeFixed[bdry_node1] = true;
			//coordInfoVecs.isNodeFixed[bdry_node2] = true;
		}
		else if (generalParams.triangles_in_upperhem[T1] != 1 && generalParams.triangles_in_upperhem[T2] == 1){
			generalParams.boundaries_in_upperhem[i] = 1;
			//std::cout<<generalParams.boundaries_in_upperhem[i]<<std::endl;
		//	generalParams.triangles_in_upperhem[T1] = 0;
		//	generalParams.triangles_in_upperhem[T2] = 0;
			double bdry_node1 = coordInfoVecs.edges2Nodes_1[i];
			double bdry_node2 = coordInfoVecs.edges2Nodes_2[i];
			boundary_node_list.push_back(bdry_node1);
			boundary_node_list.push_back(bdry_node2);
			boundary_edge_list.push_back(i);
			//generalParams.nodes_in_upperhem[bdry_node1] = 0;
			//generalParams.nodes_in_upperhem[bdry_node2] = 0;
			//coordInfoVecs.isNodeFixed[bdry_node1] = true;
			//coordInfoVecs.isNodeFixed[bdry_node2] = true;
		}
		else {
			generalParams.boundaries_in_upperhem[i] = -1;
			//std::cout<<generalParams.boundaries_in_upperhem[i]<<std::endl;
		}
	}
	std::cout<<"size of boundary_node_list (this is double-counted) = "<<boundary_node_list.size()<<std::endl;
	//generalParams.eq_total_boundary_length = generalParams.boundaries_in_upperhem.size()*generalParams.Rmin;

	/*for (int i = 0; i < coordInfoVecs.num_edges; i++){
		int aaa = coordInfoVecs.edges2Nodes_1[i];
		int bbb = coordInfoVecs.edges2Nodes_2[i];
		if (aaa == 1 && bbb == 1){
			generalParams.edges_in_upperhem[i] = 1;
			generalParams.edges_in_upperhem_list.push_back(i);
		}
		else if (aaa == 1 || bbb == 1){
			generalParams.edges_in_upperhem[i] = 0;
		}
		else{
			generalParams.edges_in_upperhem[i] = -1;
		}
		
	}*/
	
	

	int true_num_edges_in_upperhem = 0;
	for (int i = 0; i < coordInfoVecs.num_edges; i++){
		if (generalParams.edges_in_upperhem_list[i] != INT_MAX && generalParams.edges_in_upperhem_list[i] >= 0){
		true_num_edges_in_upperhem += 1;
		}
	}
	

	//std::vector<int> edge_to_ljparticle;
	//generalParams.edge_to_ljparticle.reserve(coordInfoVecs.num_edges);
	for (int i = 0; i < coordInfoVecs.num_edges; i++){
		generalParams.edge_to_ljparticle.push_back(-1);
	};
	/////////////////////////////////////////////////////////////////////
	////////////// END OF IDENTIFYING REG. WITH DIFF. MECH PROP /////////
	/////////////////////////////////////////////////////////////////////


	//std::cout<<"ERROR HERE?"<<std::endl;
	ComputeVolume(
		generalParams,
		coordInfoVecs,
		linearSpringInfoVecs,
		ljInfoVecs
	);
	ComputeAreaBud(
		generalParams,
		coordInfoVecs
		);
	//std::cout<<"ERROR HERE 2?"<<std::endl;
	double initial_volume;
	initial_volume = generalParams.true_current_total_volume;
	generalParams.eq_total_volume = generalParams.true_current_total_volume*VOLUME_FACTOR;//This is for setting different equilibrium volume to mimic growth or shirnkage.
	std::cout<<"true_current_total_volume = "<<generalParams.true_current_total_volume<<std::endl;
	std::cout<<"eq_total_volume = "<<generalParams.eq_total_volume<<std::endl;
	std::cout<<"true_current_bud_area = "<<generalParams.current_bud_area<<std::endl;

	//////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////// START OF ACTUAL SIMULATION /////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////

	/* Build the initial gradient weakend scale */
// if (generalParams.SCALE_TYPE == 4){
// 	dtb = 0.0;//dtb := distance to boundary
// 	generalParams.septin_ring_z = 0.0;
// 	generalParams.boundary_z = 0.0;
// 	for (int k = 0; k < boundary_edge_list.size(); k++){
// 	//for (int k = 0; k < boundary_node_list.size(); k++){
// 		double n1 = coordInfoVecs.edges2Nodes_1[boundary_edge_list[k]];
// 		double n2 = coordInfoVecs.edges2Nodes_2[boundary_edge_list[k]];
// 		double cent_of_edge_x = (coordInfoVecs.nodeLocX[n1] + coordInfoVecs.nodeLocX[n2])/2.0;
// 		double cent_of_edge_y = (coordInfoVecs.nodeLocY[n1] + coordInfoVecs.nodeLocY[n2])/2.0;
// 		double cent_of_edge_z = (coordInfoVecs.nodeLocZ[n1] + coordInfoVecs.nodeLocZ[n2])/2.0;
// 		double dist_x = coordInfoVecs.nodeLocX[max_height_index] - cent_of_edge_x;
// 		double dist_y = coordInfoVecs.nodeLocY[max_height_index] - cent_of_edge_y;
// 		double dist_z = coordInfoVecs.nodeLocZ[max_height_index] - cent_of_edge_z;
// 		// double temp_dist = sqrt((coordInfoVecs.nodeLocX[max_height_index] - coordInfoVecs.nodeLocX[n1])*(coordInfoVecs.nodeLocX[max_height_index] - coordInfoVecs.nodeLocX[n1]) +
// 		// (coordInfoVecs.nodeLocY[max_height_index] - coordInfoVecs.nodeLocY[n1])*(coordInfoVecs.nodeLocY[max_height_index] - coordInfoVecs.nodeLocY[n1]) +
// 		// 	(coordInfoVecs.nodeLocZ[max_height_index] - coordInfoVecs.nodeLocZ[n1])*(coordInfoVecs.nodeLocZ[max_height_index] - coordInfoVecs.nodeLocZ[n1]));
// 		double temp_dist = sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);
// 		generalParams.septin_ring_z += (coordInfoVecs.nodeLocZ[n1] + coordInfoVecs.nodeLocZ[n2]);
// 		if (temp_dist >= dtb){
// 			dtb = temp_dist;
// 			/* "dtb" will be used to identify where the septin ring is located, and used to determine the Hill coefficient*/
// 		}
// 	}

// 	//std::cout<<"dtb = "<<dtb<<std::endl;
// 	generalParams.septin_ring_z = generalParams.septin_ring_z/(1.0*boundary_node_list.size());
// // 	//generalParams.boundary_z = generalParams.septin_ring_z - generalParams.Rmin;
// 	/* dtb will be only calculated once so we can effectively keep the Hill eqn curve consistent with only horizontal shift */
// 	//dtb_max = dtb + (generalParams.Rmin);
	
// 	//std::cout<<"initial distance between cell tip and the boundary of weakened area = "<<dtb<<std::endl;
// 	//std::cout<<"Notice that here, the distance from the tip to the boundary is slightly extended by half of the equilibrium length of an edge"<<std::endl;
// 	//std::cout<<"If this message is present, we are forcing a fixed portion of the bud tip to be occupied by the max concentration"<<std::endl;
// 	//generalParams.hilleqnconst = (dtb + generalParams.Rmin/4.0)/dtb_max;
// 	//generalParams.hilleqnconst = dtb/dtb_max;
// 	//generalParams.hilleqnpow = 35.0;
// 	//std::cout<<"hill equation constant K = "<<generalParams.hilleqnconst<<std::endl;
// 	// std::cout<<"hill (equation) coefficient = "<<generalParams.hilleqnpow<<std::endl;
// 	// std::cout<<"NOTE: IN THIS SIMULATION, THE LOCATION WHERE 50% WEAKENING IS EXPERIENCED IS LOCATED SLIGHTLY AWAY FROM THE SEPTIN RING, "<<std::endl;
// 	// std::cout<<"THIS IS DUE TO THE FACT THAT IN ISOTROPIC CASE, SEPTIN RING LOCATION MUST BE SUFFICIENTLY WEAKENED TO INDUCE BUDDING"<<std::endl;
// 	// std::cout<<" "<<std::endl;
// 	// std::cout<<" "<<std::endl;
// 	// std::cout<<" "<<std::endl;
// 	// std::cout<<" "<<std::endl;
// 	// std::cout<<" "<<std::endl;
// 	// std::cout<<" "<<std::endl;
// 	// std::cout<<" "<<std::endl;
// 	// std::cout<<" "<<std::endl;

// }
	edgeswap_ptr->transferDtoH(generalParams, coordInfoVecs, build_ptr->hostSetInfoVecs);
	// edgeswap_ptr->gradient_weakening_update_host_vecs(sigma,
	// 	//max_height_index,
	// 	coordInfoVecs.nodeLocX[max_height_index],
	// 	coordInfoVecs.nodeLocY[max_height_index],
	// 	coordInfoVecs.nodeLocZ[max_height_index],
	// 	dtb,
	// 	dtb_max,
	// 	generalParams,
	// 	coordInfoVecs,
	// 	build_ptr->hostSetInfoVecs);

	generalParams.septin_ring_z = 0.0;
	generalParams.boundary_z = 0.0;
	for (int k = 0; k < boundary_edge_list.size(); k++){
	//for (int k = 0; k < boundary_node_list.size(); k++){
		double n1 = coordInfoVecs.edges2Nodes_1[boundary_edge_list[k]];
		double n2 = coordInfoVecs.edges2Nodes_2[boundary_edge_list[k]];
		double cent_of_edge_x = (coordInfoVecs.nodeLocX[n1] + coordInfoVecs.nodeLocX[n2])/2.0;
		double cent_of_edge_y = (coordInfoVecs.nodeLocY[n1] + coordInfoVecs.nodeLocY[n2])/2.0;
		double cent_of_edge_z = (coordInfoVecs.nodeLocZ[n1] + coordInfoVecs.nodeLocZ[n2])/2.0;
		double dist_x = coordInfoVecs.nodeLocX[max_height_index] - cent_of_edge_x;
		double dist_y = coordInfoVecs.nodeLocY[max_height_index] - cent_of_edge_y;
		double dist_z = coordInfoVecs.nodeLocZ[max_height_index] - cent_of_edge_z;
		// double temp_dist = sqrt((coordInfoVecs.nodeLocX[max_height_index] - coordInfoVecs.nodeLocX[n1])*(coordInfoVecs.nodeLocX[max_height_index] - coordInfoVecs.nodeLocX[n1]) +
		// (coordInfoVecs.nodeLocY[max_height_index] - coordInfoVecs.nodeLocY[n1])*(coordInfoVecs.nodeLocY[max_height_index] - coordInfoVecs.nodeLocY[n1]) +
		// 	(coordInfoVecs.nodeLocZ[max_height_index] - coordInfoVecs.nodeLocZ[n1])*(coordInfoVecs.nodeLocZ[max_height_index] - coordInfoVecs.nodeLocZ[n1]));
		double temp_dist = sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);
		generalParams.septin_ring_z += (coordInfoVecs.nodeLocZ[n1] + coordInfoVecs.nodeLocZ[n2]);
		if (temp_dist >= dtb){
			dtb = temp_dist;
			/* "dtb" will be used to identify where the septin ring is located, and used to determine the Hill coefficient*/
		}
	}
	generalParams.septin_ring_z = generalParams.septin_ring_z/(1.0*boundary_node_list.size());
	std::cout<<"septin_ring_z = "<<generalParams.septin_ring_z<<std::endl;
	std::cout<<"dtb = "<<dtb<<std::endl;
	generalParams.boundary_z = generalParams.septin_ring_z;
	edgeswap_ptr->gradient_weakening_update_host_vecs_tip(sigma,
		//max_height_index,
		coordInfoVecs.nodeLocX[max_height_index],
		coordInfoVecs.nodeLocY[max_height_index],
		coordInfoVecs.nodeLocZ[max_height_index],
		dtb,
		distance_uniform_weak,
		generalParams,
		coordInfoVecs,
		build_ptr->hostSetInfoVecs);
	for (int u = 0; u < generalParams.maxNodeCount; u++){
		int BETA = edgeswap_ptr->nodes2Triangles_host_vecs(
			u,
			build_ptr->hostSetInfoVecs,
			coordInfoVecs,
			generalParams,
			auxVecs);
	}
	edgeswap_ptr->transferHtoD(generalParams, coordInfoVecs, build_ptr->hostSetInfoVecs);
	/*for (int h = 0; h < coordInfoVecs.num_edges; h++){
		std::cout<<coordInfoVecs.scaling_per_edge[h]<<std::endl;
		double scaling = 0.0;//(spring_constant_weak/spring_constant);
		double what_spring_constant = bendingTriangleInfoVecs.spring_constant*((1.0/(1.0+pow(0.5/coordInfoVecs.scaling_per_edge[h], 6.0)))*(1-scaling) + scaling);
		if (what_spring_constant < bendingTriangleInfoVecs.spring_constant_weak){what_spring_constant = bendingTriangleInfoVecs.spring_constant_weak;}
		std::cout<<"bend_constant = "<<what_spring_constant<<std::endl;
	}
	std::cout<<"end of scaling_per_edge printout"<<std::endl;*/
	

	
	/*for (int i = 0; i < 10; i++){
		std::cout<<"nodess2Triangles_1[ "<<i<<" ] = "<<coordInfoVecs.nodes2Triangles_1[i]<<std::endl;
		std::cout<<"nodess2Triangles_2[ "<<i<<" ] = "<<coordInfoVecs.nodes2Triangles_2[i]<<std::endl;
		std::cout<<"nodess2Triangles_3[ "<<i<<" ] = "<<coordInfoVecs.nodes2Triangles_3[i]<<std::endl;
		std::cout<<"nodess2Triangles_4[ "<<i<<" ] = "<<coordInfoVecs.nodes2Triangles_4[i]<<std::endl;
		std::cout<<"nodess2Triangles_5[ "<<i<<" ] = "<<coordInfoVecs.nodes2Triangles_5[i]<<std::endl;
		std::cout<<"nodess2Triangles_6[ "<<i<<" ] = "<<coordInfoVecs.nodes2Triangles_6[i]<<std::endl;
		std::cout<<"nodess2Triangles_7[ "<<i<<" ] = "<<coordInfoVecs.nodes2Triangles_7[i]<<std::endl;
		std::cout<<"nodess2Triangles_8[ "<<i<<" ] = "<<coordInfoVecs.nodes2Triangles_8[i]<<std::endl;
		std::cout<<"nodess2Triangles_9[ "<<i<<" ] = "<<coordInfoVecs.nodes2Triangles_9[i]<<std::endl;

	}*/
	while (runSim == true){
		//WHEN += 1;
		double current_time = 0.0;
		//nodenormal_1.resize(generalParams.maxNodeCount);
		//nodenormal_2.resize(generalParams.maxNodeCount);
		//nodenormal_3.resize(generalParams.maxNodeCount);
		//std::fill(nodenormal_1.begin(), nodenormal_1.end(), 0.0);
		//std::fill(nodenormal_2.begin(), nodenormal_2.end(), 0.0);
		//std::fill(nodenormal_3.begin(), nodenormal_3.end(), 0.0);
		/* for (int k = 0; k < coordInfoVecs.num_triangles; k++){
			if (coordInfoVecs.triangles2Nodes_1[k] != INT_MAX || coordInfoVecs.triangles2Nodes_3[k] != INT_MAX || coordInfoVecs.triangles2Nodes_3[k] != INT_MAX){
				double x1 = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[k]];
				double y1 = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[k]];
				double z1 = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[k]];
				double x2 = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_2[k]];
				double y2 = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_2[k]];
				double z2 = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_2[k]];
				double x3 = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_3[k]];
				double y3 = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_3[k]];
				double z3 = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_3[k]];
				double nx = (y2 - y1)*(z3 - z1) - (y3 - y1)*(z2 - z1);
				double ny = -((x2 - x1)*(z3 - z1) - (x3 - x1)*(z2 - z1));
				double nz = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1);
				nodenormal_1[coordInfoVecs.triangles2Nodes_1[k]] += nx;
				nodenormal_2[coordInfoVecs.triangles2Nodes_1[k]] += ny;
				nodenormal_3[coordInfoVecs.triangles2Nodes_1[k]] += nz;
				nodenormal_1[coordInfoVecs.triangles2Nodes_2[k]] += nx;
				nodenormal_2[coordInfoVecs.triangles2Nodes_2[k]] += ny;
				nodenormal_3[coordInfoVecs.triangles2Nodes_2[k]] += nz;
				nodenormal_1[coordInfoVecs.triangles2Nodes_3[k]] += nx;
				nodenormal_2[coordInfoVecs.triangles2Nodes_3[k]] += ny;
				nodenormal_3[coordInfoVecs.triangles2Nodes_3[k]] += nz;
			}
			else{continue;}
		}
		for (int k = 0; k < nodenormal_1.size(); k++){
			double UN = sqrt((nodenormal_1[k]*nodenormal_1[k]) + (nodenormal_2[k]*nodenormal_2[k]) + (nodenormal_3[k]*nodenormal_3[k]));
			nodenormal_1[k] = nodenormal_1[k]/UN;
			nodenormal_2[k] = nodenormal_2[k]/UN;
			nodenormal_3[k] = nodenormal_3[k]/UN;
		} */


		
		//generalParams.kT = 1.0;//reset kT before simulations starts.
		//Max_Runtime = 0.0;//2.5;
		int translate_counter = 0;
		
			while (current_time < 0.0*(Max_Runtime)){
					translate_counter += 1;
					Solve_Forces();
					//std::cout<<"ERROR 1"<<std::endl;
					// if (generalParams.true_current_total_volume/initial_volume >= LINE_TENSION_THRESHOLD){
					// 	ComputeLineTensionSprings(
					// 		generalParams,
					// 		coordInfoVecs,
					// 		linearSpringInfoVecs);
					// 	}
				
					energy_rep =
					ComputeMemRepulsionEnergy(
						coordInfoVecs,
						linearSpringInfoVecs, 
						capsidInfoVecs,
						generalParams,
						auxVecs);

					//now forces are computed, move nodes.
					
					

					

					double beta;
					

				/*for (int k = 0; k < generalParams.maxNodeCount; k++){
					if (generalParams.nodes_in_upperhem[k] == 1 && generalParams.boundaries_in_upperhem[k] != 1){
						coordInfoVecs.nodeForceX[k] = 0.001*coordInfoVecs.nodeForceX[k];
						coordInfoVecs.nodeForceY[k] = 0.001*coordInfoVecs.nodeForceY[k];
						coordInfoVecs.nodeForceZ[k] = 0.001*coordInfoVecs.nodeForceZ[k];
					}
				}*/

				/*for (int k = 0; k < generalParams.maxNodeCount; k++){
					coordInfoVecs.nodeForceX[k] += generalParams.volume_spring_constant*coordInfoVecs.SurfaceNormalX[k];
					coordInfoVecs.nodeForceY[k] += generalParams.volume_spring_constant*coordInfoVecs.SurfaceNormalY[k];
					coordInfoVecs.nodeForceZ[k] += generalParams.volume_spring_constant*coordInfoVecs.SurfaceNormalZ[k];
				}*/

				/*ComputeSurfaceNormal(
					coordInfoVecs,
					generalParams,
					auxVecs
				);*/
				/*for (int u = 0; u < generalParams.maxNodeCount; u++){
					int GAMMA = edgeswap_ptr->surfaceNormal_device_vecs(
						u,
						coordInfoVecs,
						generalParams
					);
				}*/
				AdvancePositions(
					coordInfoVecs,
					generalParams,
					domainParams);

/*				if (translate_counter % translate_frequency == 1){

					newcenterX = 0.0;
					newcenterY = 0.0;
					newcenterZ = 0.0;
					for (int i = 0; i < generalParams.maxNodeCount; i++){//for (int i = 0; i < coordInfoVecs.nodeLocX.size(); i++){
						newcenterX += coordInfoVecs.nodeLocX[i];
						newcenterY += coordInfoVecs.nodeLocY[i];
						newcenterZ += coordInfoVecs.nodeLocZ[i];
					}
					newcenterX = newcenterX/generalParams.maxNodeCount; //coordInfoVecs.nodeLocX.size();
					newcenterY = newcenterY/generalParams.maxNodeCount; //coordInfoVecs.nodeLocX.size();
					newcenterZ = newcenterZ/generalParams.maxNodeCount; //coordInfoVecs.nodeLocX.size();
					displacementX = newcenterX - generalParams.centerX;
					displacementY = newcenterY - generalParams.centerY;
					displacementZ = newcenterZ - generalParams.centerZ;

					for (int i = 0; i < generalParams.maxNodeCount; i++){
						coordInfoVecs.nodeLocX[i] += -displacementX;
						coordInfoVecs.nodeLocY[i] += -displacementY;
						coordInfoVecs.nodeLocZ[i] += -displacementZ;
					}
					for (int i = 0; i < ljInfoVecs.LJ_PosX_all.size(); i++){
						ljInfoVecs.LJ_PosX_all[i] += -displacementX;
						ljInfoVecs.LJ_PosY_all[i] += -displacementY;
						ljInfoVecs.LJ_PosZ_all[i] += -displacementZ;
					}

					//Here we re-establish the new filament base according to the current location of nuclei nodes
					 int maxElementIndex = std::max_element(coordInfoVecs.nodeLocZ.begin(),coordInfoVecs.nodeLocZ.end()) - coordInfoVecs.nodeLocZ.begin();
					for (int i = 0; i < generalParams.maxNodeCountLJ; i++){
						if (i == 0){
							filament_base[i] = maxElementIndex;
							continue;
						}
						for (int j = 0; j < generalParams.maxNodeCount; j++){
							double xsquared = (ljInfoVecs.LJ_PosX_all[i] - coordInfoVecs.nodeLocX[j])*
												(ljInfoVecs.LJ_PosX_all[i] - coordInfoVecs.nodeLocX[j]);
							double ysquared = (ljInfoVecs.LJ_PosY_all[i] - coordInfoVecs.nodeLocY[j])*
												(ljInfoVecs.LJ_PosY_all[i] - coordInfoVecs.nodeLocY[j]);
							double zsquared = (ljInfoVecs.LJ_PosZ_all[i] - coordInfoVecs.nodeLocZ[j])*
												(ljInfoVecs.LJ_PosZ_all[i] - coordInfoVecs.nodeLocZ[j]);
							double R = sqrt(xsquared + ysquared + zsquared);
							if (R < (max_height - min_height)/2.0 && j != maxElementIndex){
								filament_base[i] = j;
								break;
							}
							else{filament_base[i] = -1;}
						}
					} 
				}*/
							
					new_total_energy = linearSpringInfoVecs.linear_spring_energy + 
						areaTriangleInfoVecs.area_triangle_energy + 
						bendingTriangleInfoVecs.bending_triangle_energy + 
						0.5*energy_rep;// + 
						//ljInfoVecs.lj_energy_M +
						//ljInfoVecs.lj_energy_LJ +
						//generalParams.volume_energy;

				energy_gradient = sqrt((new_total_energy - old_total_energy)*(new_total_energy - old_total_energy))/old_total_energy;
				if (current_time >= Max_Runtime*0.25 && energy_gradient < energy_gradient_threshold){
					break;
				}
				old_total_energy = new_total_energy;
				current_time+=generalParams.dt;
				

			}
			
		   
			/*max_height = -10000.0;
			min_height = 10000.0;
			for (int k = 0; k < generalParams.maxNodeCount; k++){
				if (coordInfoVecs. nodeLocZ[k] >= max_height){
					max_height = coordInfoVecs. nodeLocZ[k];
				}
				if (coordInfoVecs.nodeLocZ[k] <= min_height){
					min_height = coordInfoVecs.nodeLocZ[k];
				}
			}*/

		std::cout<<"current time (1st iter before edgeswap): "<< current_time << std::endl;
		std::cout<<"current total energy (1st iter before edgeswap) = "<<new_total_energy<<std::endl;
		std::cout<<"LINEAR ENERGY = "<<linearSpringInfoVecs.linear_spring_energy<<std::endl;
		std::cout<<"BEND ENERGY = "<<bendingTriangleInfoVecs.bending_triangle_energy<<std::endl;
		std::cout<<"AREA ENERGY = "<<areaTriangleInfoVecs.area_triangle_energy<<std::endl;
		std::cout<<"REPULSION ENERGY = "<<energy_rep<<std::endl;
		std::cout<<"VOLUME ENERGY = "<<generalParams.volume_energy<<std::endl;
		std::cout<<"true_current_total_volume = "<<generalParams.true_current_total_volume<<std::endl;
		std::cout<<"eq_total_volume = "<<generalParams.eq_total_volume<<std::endl;
		std::cout<<"true_current_bud_area = "<<generalParams.current_bud_area<<std::endl;
		std::cout<<"current KBT = "<<generalParams.kT<<std::endl;
		if (isnan(new_total_energy)==1){
			std::cout<<"Nan or Inf position update !!!!"<<std::endl;
			runSim = false;
			break;
		}
	
		//edgeswap_ptr->transferDtoH(coordInfoVecs, build_ptr->hostSetInfoVecs);//Currently this is treated as a backup of coordInfoVecs
		storage->print_VTK_File();
		//storage->storeVariables();
		//storage->storeVariables();
		//runSim = false;
		//break;

		int edgeswap_iteration = 0;
		//double preswap_energy = new_total_energy;
		//double postswap_energy;
		//double Ediff = 0.0;
		//initial_kT = generalParams.kT;
		num_edge_loop = 0;//round(true_num_edges_in_upperhem*SAMPLE_SIZE);
		//if (num_edge_loop == 0){
		//	num_edge_loop = 1;
		//}	
		
		std::cout<<"if SAMPLE_SIZE = 0, this implies that all edges in the budding site will be tested for edgeswap"<<std::endl;
		

		

		int LINE_TENSION_START = 0;
		
		bool WEAKENED_START = false;
		bool EDGESWAP_ALGORITHM_TRIGGERED;
 		while (initial_kT > 0){
			
			/*if (generalParams.true_current_total_volume/initial_volume >= VOLUME_THRESHOLD && WEAKENED_START == false){
				linearSpringInfoVecs.spring_constant_weak = scale_linear;
				bendingTriangleInfoVecs.spring_constant_weak = scale_bend;
				areaTriangleInfoVecs.spring_constant_weak = scale_area;
				std::cout<<"membrane weakening initiated"<<std::endl;
				WEAKENED_START = true;
			}*/
 					////////////////////NOW RELAX THE ATTEMPTED EDGESWAP//////////////////////
					 current_time = 0.0;
					 translate_counter = 0;
					 double VOLUME_RATIO = generalParams.true_current_total_volume/generalParams.eq_total_volume;
					//if (VOLUME_RATIO > 0.75 && VOLUME_FACTOR <= 2.5){
					//	VOLUME_FACTOR += 0.2;
					//	generalParams.eq_total_volume = initial_volume*VOLUME_FACTOR;
					//};
					if (edgeswap_iteration >= 10){
					//if (generalParams.true_current_total_volume/initial_volume >= LINE_TENSION_THRESHOLD){
						if (LINE_TENSION_START < 1){
							double DIST = 0.0;
							double COUNT = 0.0;
							for (int t = 0; t < coordInfoVecs.num_edges; t++){
								if (generalParams.boundaries_in_upperhem[t] == 1){
									COUNT += 1.0;
									int node1 = coordInfoVecs.edges2Nodes_1[t];
									int node2 = coordInfoVecs.edges2Nodes_2[t];
									DIST += sqrt((coordInfoVecs.nodeLocX[node2] - coordInfoVecs.nodeLocX[node1])*(coordInfoVecs.nodeLocX[node2] - coordInfoVecs.nodeLocX[node1]) +
									(coordInfoVecs.nodeLocY[node2] - coordInfoVecs.nodeLocY[node1])*(coordInfoVecs.nodeLocY[node2] - coordInfoVecs.nodeLocY[node1]) + 
									(coordInfoVecs.nodeLocZ[node2] - coordInfoVecs.nodeLocZ[node1])*(coordInfoVecs.nodeLocZ[node2] - coordInfoVecs.nodeLocZ[node1]));
								}
							}
							generalParams.length_scale = (DIST/COUNT)/generalParams.Rmin;
							std::cout<<"equilibrium length of each segment of the septin ring = "<<generalParams.length_scale*generalParams.Rmin<<std::endl;
							generalParams.eq_total_boundary_length = COUNT*generalParams.length_scale* generalParams.Rmin;
							std::cout<<"equilibrium length of the septin ring = "<<generalParams.eq_total_boundary_length<<std::endl;
							LINE_TENSION_START += 1;
						}
						
					}

					if (edgeswap_iteration >= 10){
						linearSpringInfoVecs.spring_constant_weak = scale_linear;
						bendingTriangleInfoVecs.spring_constant_weak = scale_bend;
						areaTriangleInfoVecs.spring_constant_weak = scale_area;
					}
					else {
						linearSpringInfoVecs.spring_constant_weak = linearSpringInfoVecs.spring_constant;
						bendingTriangleInfoVecs.spring_constant_weak = bendingTriangleInfoVecs.spring_constant;
						areaTriangleInfoVecs.spring_constant_weak = areaTriangleInfoVecs.spring_constant;
					}
					//std::cout<<"start relaxation step"<<std::endl;
					EDGESWAP_ALGORITHM_TRIGGERED = false;
					bool end_of_relaxation = false;
 					while (current_time < Max_Runtime){
						
						 if (Max_Runtime <= 0.0){
							 std::cout<<"Max_Runtime is set to be 0 or negative! "<<std::endl;
							 break;
						 }
						 
						 translate_counter += 1;
						 //std::cout<<"STOPPED BEFORE Solve_Forces"<<std::endl;
						 Solve_Forces();
						 //std::cout<<"ERROR HERE 2"<<std::endl;
						 //if (generalParams.true_current_total_volume/initial_volume >= LINE_TENSION_THRESHOLD){
						if (LINE_TENSION_START >= 1){
							ComputeLineTensionSprings(
								generalParams,
								coordInfoVecs,
								linearSpringInfoVecs);
							}
						//std::cout<<"STOPPED BEFORE MemRepul"<<std::endl;
 						energy_rep =
 						ComputeMemRepulsionEnergy(
 							coordInfoVecs,
 							linearSpringInfoVecs, 
 							capsidInfoVecs,
 							generalParams,
							 auxVecs);
					if ((generalParams.true_current_total_volume/initial_volume) < 0.6 || generalParams.true_current_total_volume/initial_volume >= MAX_VOLUME_RATIO){
						generalParams.true_num_edges = 0;
						for (int i = 0; i < coordInfoVecs.num_edges; i++){
							if (coordInfoVecs.edges2Nodes_1[i] != INT_MAX && coordInfoVecs.edges2Nodes_2[i] != INT_MAX){
								generalParams.true_num_edges += 1;
							}
						}
						storage-> print_VTK_File();
						//storage-> storeVariables();
						runSim = false;
						initial_kT = -0.00000000000000001;
						if (generalParams.true_current_total_volume/initial_volume < 0.6){
							std::cout<<"Cell over compression 60%"<<std::endl;
						}
						else if (generalParams.true_current_total_volume/initial_volume >= MAX_VOLUME_RATIO){
							std::cout<<"Target volume ratio exceeded. Current volume ratio = "<<generalParams.true_current_total_volume/initial_volume<<std::endl;
							std::cout<<"true_current_bud_area = "<<generalParams.current_bud_area<<std::endl;
						}
						std::cout<<"Total iterative steps: "<<edgeswap_iteration<<" x Max_Runtime"<<std::endl;
						Max_Runtime = 0.0;
						runSim = false;
						initial_kT = -0.00000001;
						break;

						}
						//std::cout<<"ERROR HERE 3"<<std::endl;
 						//now forces are computed, move nodes.
						 double beta;
						
						 /*for (int k = 0; k < generalParams.maxNodeCount; k++){
							if (generalParams.nodes_in_upperhem[k] == 1 && generalParams.boundaries_in_upperhem[k] != 1){
								coordInfoVecs.nodeForceX[k] = 0.001*coordInfoVecs.nodeForceX[k];
								coordInfoVecs.nodeForceY[k] = 0.001*coordInfoVecs.nodeForceY[k];
								coordInfoVecs.nodeForceZ[k] = 0.001*coordInfoVecs.nodeForceZ[k];
							}
						}*/
						
						/*for (int k = 0; k < generalParams.maxNodeCount; k++){
						
							coordInfoVecs.nodeForceX[k] += generalParams.volume_spring_constant*coordInfoVecs.SurfaceNormalX[k];
							coordInfoVecs.nodeForceY[k] += generalParams.volume_spring_constant*coordInfoVecs.SurfaceNormalY[k];
							coordInfoVecs.nodeForceZ[k] += generalParams.volume_spring_constant*coordInfoVecs.SurfaceNormalZ[k];
						  }*/
						  
						  
						  
						//  std::cout<<"STOPPED BEFORE surfacenormal"<<std::endl;
						//std::cout<<"IS IT ADVANCE POSITION PROBLEM?"<<std::endl;
						/*for (int u = 0; u < generalParams.maxNodeCount; u++){
							int GAMMA = edgeswap_ptr->surfaceNormal_device_vecs(
								u,
								coordInfoVecs,
								generalParams
							);
						}*/
						//std::cout<<"STOPPED BEFORE AdvancePos"<<std::endl;
 						AdvancePositions(
 							coordInfoVecs,
 							generalParams,
							 domainParams);
							// std::cout<<"ERROR HERE 4"<<std::endl;
						new_total_energy = linearSpringInfoVecs.linear_spring_energy + 
								areaTriangleInfoVecs.area_triangle_energy + 
								bendingTriangleInfoVecs.bending_triangle_energy +
								0.5*energy_rep;// +
					//ljInfoVecs.lj_energy_M +  
					// ljInfoVecs.lj_energy_LJ +
						//generalParams.volume_energy;
				//std::cout<<"new_total_energy = "<<new_total_energy<<std::endl;

				energy_gradient = sqrt((new_total_energy - old_total_energy)*(new_total_energy - old_total_energy))/old_total_energy;
				old_total_energy = new_total_energy;
				current_time+=generalParams.dt;
				//std::cout<<"ERROR HERE 5"<<std::endl;
						
						

					
					// if (generalParams.SCALE_TYPE != 3){
						if (translate_counter % (translate_frequency*1) == 0 || end_of_relaxation == true){
							// max_height = -10000.0;
							// double current_center_x = 0.0;
							// double current_center_y = 0.0;
							
							// for (int k = 0; k < generalParams.maxNodeCount; k++){
							// 	if (generalParams.nodes_in_upperhem[k] != INT_MAX){
							// 		current_center_x += coordInfoVecs.nodeLocX[k];
							// 		current_center_y += coordInfoVecs.nodeLocY[k];
							// 	}
								
							// 	if (coordInfoVecs. nodeLocZ[k] >= max_height){
							// 		max_height = coordInfoVecs.nodeLocZ[k];
							// 		max_height_index = k;
							// 	}
						
							// }
							// current_center_x = current_center_x/generalParams.maxNodeCount;
							// current_center_y = current_center_y/generalParams.maxNodeCount;
							
							// //std::cout<<"max_height_index = "<<max_height_index<<std::endl;
							// dtb = 0.0;//dtb := distance to boundary
							// generalParams.septin_ring_z = 0.0;
							// generalParams.boundary_z = 0.0;
							// for (int k = 0; k < boundary_edge_list.size(); k++){
							// //for (int k = 0; k < boundary_node_list.size(); k++){
							// 	double n1 = coordInfoVecs.edges2Nodes_1[boundary_edge_list[k]];
							// 	double n2 = coordInfoVecs.edges2Nodes_2[boundary_edge_list[k]];
							// 	double cent_of_edge_x = (coordInfoVecs.nodeLocX[n1] + coordInfoVecs.nodeLocX[n2])/2.0;
							// 	double cent_of_edge_y = (coordInfoVecs.nodeLocY[n1] + coordInfoVecs.nodeLocY[n2])/2.0;
							// 	double cent_of_edge_z = (coordInfoVecs.nodeLocZ[n1] + coordInfoVecs.nodeLocZ[n2])/2.0;
							// 	double dist_x = current_center_x - cent_of_edge_x;
							// 	double dist_y = current_center_y - cent_of_edge_y;
							// 	double dist_z = max_height - cent_of_edge_z;
							// 	// double temp_dist = sqrt((coordInfoVecs.nodeLocX[max_height_index] - coordInfoVecs.nodeLocX[n1])*(coordInfoVecs.nodeLocX[max_height_index] - coordInfoVecs.nodeLocX[n1]) +
							// 	// (coordInfoVecs.nodeLocY[max_height_index] - coordInfoVecs.nodeLocY[n1])*(coordInfoVecs.nodeLocY[max_height_index] - coordInfoVecs.nodeLocY[n1]) +
							// 	// 	(coordInfoVecs.nodeLocZ[max_height_index] - coordInfoVecs.nodeLocZ[n1])*(coordInfoVecs.nodeLocZ[max_height_index] - coordInfoVecs.nodeLocZ[n1]));
							// 	// generalParams.septin_ring_z += coordInfoVecs.nodeLocZ[n1];
							// 	double temp_dist = sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);
							// 	if (temp_dist >= dtb){
							// 		dtb = temp_dist;
							// 		/* "dtb" will be used to identify where the septin ring is located, and used to determine the Hill coefficient*/
							// 	}
							// }
							// //std::cout<<"dtb = "<<dtb<<std::endl;
							// generalParams.septin_ring_z = generalParams.septin_ring_z/boundary_node_list.size();
							// generalParams.boundary_z = generalParams.septin_ring_z - generalParams.Rmin;
							// /* dtb will be only calculated once so we can effectively keep the Hill eqn curve consistent with only horizontal shift */
							// dtb_max = dtb + (generalParams.Rmin);
							// // generalParams.septin_ring_z = 0.0;
							// // generalParams.boundary_z = 0.0;
							// // //for (int k = 0; k < boundary_edge_list.size(); k++){
							// // for (int k = 0; k < boundary_node_list.size(); k++){
							// // 	double n1 = boundary_node_list[k];//coordInfoVecs.edges2Nodes_1[boundary_edge_list[k]];
							// // 	generalParams.septin_ring_z += coordInfoVecs.nodeLocZ[n1];
							// // }
							// //generalParams.septin_ring_z = generalParams.septin_ring_z/boundary_node_list.size();
							// //generalParams.boundary_z = generalParams.septin_ring_z - generalParams.Rmin;
							// /* dtb will be only calculated once so we can effectively keep the Hill eqn curve consistent with only horizontal shift */
					
							// //generalParams.hilleqnconst = (dtb + generalParams.Rmin/4.0)/dtb_max;
							// generalParams.hilleqnconst = dtb/dtb_max;

							// edgeswap_ptr->transferDtoH(generalParams, coordInfoVecs, build_ptr->hostSetInfoVecs);//Currently this is treated as a backup of coordInfoVecs
							// edgeswap_ptr->gradient_weakening_update_host_vecs(sigma,
							// 	current_center_x,
							// 	current_center_y,
							// 	max_height,
							// 	dtb,
							// 	dtb_max,
							// 	generalParams,
							// 	coordInfoVecs,
							// 	build_ptr->hostSetInfoVecs);
							// edgeswap_ptr->transferHtoD(generalParams, coordInfoVecs, build_ptr->hostSetInfoVecs);//Currently this is treated as a backup of coordInfoVecs
							}
					// }	

				}
					//std::cout<<"current_time (# of relaxation step) = "<<current_time<<std::endl;
					//if (EDGESWAP_ALGORITHM_TRIGGERED == false){
					if (1 > 0){
						//EDGE_SWAP IS TRIGGERED HERE IF THE RELAXATION IN THE PREVIOUS SECTION DID NOT HIT THE THRESHOLD VALUE TO TRIGGER
						//EDGESWAP NORMALLY.
						//std::cout<<"ERROR HERE 6"<<std::endl;
						edgeswap_ptr->transferDtoH(generalParams, coordInfoVecs, build_ptr->hostSetInfoVecs);
							//std::cout<<"ERROR 1.5"<<std::endl;
							VectorShuffleForEdgeswapLoop.clear();
							for (int i = 0; i < coordInfoVecs.num_edges; i++){
								if (//generalParams.edges_in_tip[i] == 1 &&
									generalParams.edges_in_upperhem_list[i] >= 0 && 
									generalParams.edges_in_upperhem_list[i] < (INT_MAX-100) &&
									generalParams.boundaries_in_upperhem[i] != 1)
									VectorShuffleForEdgeswapLoop.push_back(generalParams.edges_in_upperhem_list[i]);
								}	
						//	std::cout<<"STOPPED BEFORE edgeswap"<<std::endl;
							//std::random_device rand_dev;
							//std::mt19937 generator_edgeswap(rand_dev());
							num_edge_loop = round(true_num_edges_in_upperhem*SAMPLE_SIZE);
							if (num_edge_loop <= min_num_edge_loop){
								num_edge_loop = min_num_edge_loop;
							}
							//generalParams.kT = generalParams.kT*2.0;
						//	double kT_reduction = generalParams.kT/5.0;
						//while (generalParams.kT > 0.15){
							std::shuffle(std::begin(VectorShuffleForEdgeswapLoop), std::end(VectorShuffleForEdgeswapLoop), generator_edgeswap);
							//	std::shuffle(std::begin(generalParams.edges_in_upperhem_list), std::end(generalParams.edges_in_upperhem_list), generator_edgeswap);
								//for (int edge_loop = 0; edge_loop < VectorShuffleForEdgeswapLoop.size(); edge_loop++){
							for (int edge_loop = 0; edge_loop < num_edge_loop; edge_loop++) {
								bool bdry_associated_edge = false;
									//std::cout<<"edge_loop = "<<edge_loop<<std::endl;
									
									//std::random_device rand_dev;
									//std::mt19937 generator(rand_dev());
								
								std::uniform_int_distribution<int> distribution(1,VectorShuffleForEdgeswapLoop.size());
								
								int dice_roll = distribution(generator_edgeswap);
								
								int edge = VectorShuffleForEdgeswapLoop[dice_roll - 1];
								// int nd1 = coordInfoVecs.edges2Nodes_1[edge];
								// int nd2 = coordInfoVecs.edges2Nodes_2[edge];
								// for (int v = 0; v < septin_node.size(); v++){
								// 	if (nd1 == septin_node[v]){
								// 		bdry_associated_edge = true;
								// 		break;
								// 	}
								// 	else if (nd2 == septin_node[v]){
								// 		bdry_associated_edge = true;
								// 		break;
								// 	}
								// }
								//int edge = dice_roll -1;
								while (bdry_associated_edge == true || generalParams.boundaries_in_upperhem[edge] == 1 || edge == INT_MAX || edge < 0 ){
								//while (generalParams.boundaries_in_upperhem[edge] == 1 || edge == INT_MAX || edge < 0 ){
										bdry_associated_edge = false;
										dice_roll = distribution(generator_edgeswap);
										edge = VectorShuffleForEdgeswapLoop[dice_roll - 1];
										// nd1 = coordInfoVecs.edges2Nodes_1[edge];
										// nd2 = coordInfoVecs.edges2Nodes_2[edge];
										// for (int v = 0; v < septin_node.size(); v++){
										// 	if (nd1 == septin_node[v]){
										// 		bdry_associated_edge = true;
										// 		break;
										// 	}
										// 	else if (nd2 == septin_node[v]){
										// 		bdry_associated_edge = true;
										// 		break;
										// 	}
										// }
										//edge =  generalParams.edges_in_upperhem_list[dice_roll - 1];
										//edge = dice_roll -1;
									 }
									//int edge = generalParams.edges_in_upperhem_list[edge_loop];
									//int edge = VectorShuffleForEdgeswapLoop[edge_loop];
									//std::cout<<"edge = "<<edge<<std::endl;
									

									int ALPHA = edgeswap_ptr->edge_swap_host_vecs(
										edge,
										generalParams,
										build_ptr->hostSetInfoVecs,
										linearSpringInfoVecs,
										bendingTriangleInfoVecs,
										areaTriangleInfoVecs);
									
								}

								max_height = -10000.0;
								current_center_x = 0.0;
								current_center_y = 0.0;
								
								for (int k = 0; k < generalParams.maxNodeCount; k++){
									if (generalParams.nodes_in_upperhem[k] < (INT_MAX-100)){
										current_center_x += coordInfoVecs.nodeLocX[k];
										current_center_y += coordInfoVecs.nodeLocY[k];
									}
								}
								current_center_x = current_center_x/generalParams.maxNodeCount;
								current_center_y = current_center_y/generalParams.maxNodeCount;

								for (int k = 0; k < generalParams.maxNodeCount; k++){								
									if (coordInfoVecs. nodeLocZ[k] >= max_height){
										max_height = coordInfoVecs.nodeLocZ[k];
										max_height_index = k;
									}
								}

								generalParams.septin_ring_z = 0.0;
								for (int k = 0; k < boundary_edge_list.size(); k++){
								//for (int k = 0; k < boundary_node_list.size(); k++){
									double n1 = coordInfoVecs.edges2Nodes_1[boundary_edge_list[k]];
									double n2 = coordInfoVecs.edges2Nodes_2[boundary_edge_list[k]];
									double cent_of_edge_x = (coordInfoVecs.nodeLocX[n1] + coordInfoVecs.nodeLocX[n2])/2.0;
									double cent_of_edge_y = (coordInfoVecs.nodeLocY[n1] + coordInfoVecs.nodeLocY[n2])/2.0;
									double cent_of_edge_z = (coordInfoVecs.nodeLocZ[n1] + coordInfoVecs.nodeLocZ[n2])/2.0;
									double dist_x = coordInfoVecs.nodeLocX[max_height_index] - cent_of_edge_x;
									double dist_y = coordInfoVecs.nodeLocY[max_height_index] - cent_of_edge_y;
									double dist_z = coordInfoVecs.nodeLocZ[max_height_index] - cent_of_edge_z;
									// double temp_dist = sqrt((coordInfoVecs.nodeLocX[max_height_index] - coordInfoVecs.nodeLocX[n1])*(coordInfoVecs.nodeLocX[max_height_index] - coordInfoVecs.nodeLocX[n1]) +
									// (coordInfoVecs.nodeLocY[max_height_index] - coordInfoVecs.nodeLocY[n1])*(coordInfoVecs.nodeLocY[max_height_index] - coordInfoVecs.nodeLocY[n1]) +
									// 	(coordInfoVecs.nodeLocZ[max_height_index] - coordInfoVecs.nodeLocZ[n1])*(coordInfoVecs.nodeLocZ[max_height_index] - coordInfoVecs.nodeLocZ[n1]));
									double temp_dist = sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);
									generalParams.septin_ring_z += (coordInfoVecs.nodeLocZ[n1] + coordInfoVecs.nodeLocZ[n2]);
									if (temp_dist >= dtb){
										dtb = temp_dist;
										/* "dtb" will be used to identify where the septin ring is located, and used to determine the Hill coefficient*/
									}
								}
								generalParams.septin_ring_z = generalParams.septin_ring_z/(1.0*boundary_node_list.size());
								generalParams.boundary_z = generalParams.septin_ring_z;
								edgeswap_ptr->gradient_weakening_update_host_vecs_tip(sigma,
								current_center_x,
								current_center_y,
								max_height,
								dtb,
								distance_uniform_weak,
								generalParams,
								coordInfoVecs,
								build_ptr->hostSetInfoVecs);
								
							//	std::cout<<"STOPPED after edgeswap"<<std::endl;
						//		generalParams.kT -= kT_reduction;
						//	}
							//generalParams.kT = initial_kT;

						//std::cout<<"IS IT NODES2TRIANGLES PROBLEM?"<<std::endl;
						/*for (int u = 0; u < generalParams.maxNodeCount; u++){
							int BETA = edgeswap_ptr->nodes2Triangles_host_vecs(
								u,
								build_ptr->hostSetInfoVecs,
								coordInfoVecs,
								generalParams,
								auxVecs);
						}*/
						//std::cout<<"IT IS NOT"<<std::endl;
							//NOTE: EDGESWAP ALGORITHM CURRENTLY IS WRITTEN TO ALLOW AT MOST 8 NEIGHBORING NODES PER NODE.
							//std::cout<<"edgeswap done!"<<std::endl;
							edgeswap_ptr->transferHtoD(generalParams, coordInfoVecs, build_ptr->hostSetInfoVecs);
							//std::cout<<"ERROR HERE 7"<<std::endl;
						if (bud_tip_necessary == true){	
							double current_dist;
							//std::cout<<"ERROR HERE 8"<<std::endl;
								for (int i = 0; i < generalParams.maxNodeCount; i++){
									current_dist = sqrt((coordInfoVecs.nodeLocX[i] - current_center_x)*(coordInfoVecs.nodeLocX[i] - current_center_x) +
									(coordInfoVecs.nodeLocY[i] - current_center_y)*(coordInfoVecs.nodeLocY[i] - current_center_y) +
									(coordInfoVecs.nodeLocZ[i] - coordInfoVecs.nodeLocZ[max_height_index])*(coordInfoVecs.nodeLocZ[i] - coordInfoVecs.nodeLocZ[max_height_index]));
									//if (coordInfoVecs.nodeLocZ[i] > (coordInfoVecs.nodeLocZ[max_height_index] - 1.5*generalParams.Rmin)){
									if (current_dist <= 1.0){// 2.5*generalParams.Rmin){
										generalParams.nodes_in_tip[i] = 1;
									}
									else{
										generalParams.nodes_in_tip[i] = -1;
									}
									//std::cout<<"nodes "<<i<<" "<<generalParams.nodes_in_upperhem[i]<<std::endl;		
								}
								//std::cout<<"ERROR HERE 9"<<std::endl;
								for (int i = 0; i < coordInfoVecs.num_triangles; i++){
									//std::cout<<"i = "<<i<<std::endl;
									//std::cout<<coordInfoVecs.triangles2Nodes_1[i]<<" "<<coordInfoVecs.triangles2Nodes_2[i]<<" "<<coordInfoVecs.triangles2Nodes_3[i]<<std::endl;
									if (coordInfoVecs.triangles2Nodes_1[i] < 0 || coordInfoVecs.triangles2Nodes_1[i] > (INT_MAX-100)){
										generalParams.triangles_in_tip[i] = -1;
										continue;
									}
									else if (coordInfoVecs.triangles2Nodes_2[i] < 0 || coordInfoVecs.triangles2Nodes_2[i] > (INT_MAX-100)){
										generalParams.triangles_in_tip[i] = -1;
										continue;
									}
									else if (coordInfoVecs.triangles2Nodes_3[i] < 0 || coordInfoVecs.triangles2Nodes_3[i] > (INT_MAX-100)){
										generalParams.triangles_in_tip[i] = -1;
										continue;
									}

									int aaa = generalParams.nodes_in_tip[coordInfoVecs.triangles2Nodes_1[i]];
									//std::cout<<aaa<<std::endl;
									int bbb = generalParams.nodes_in_tip[coordInfoVecs.triangles2Nodes_2[i]];
									//std::cout<<bbb<<std::endl;
									int ccc = generalParams.nodes_in_tip[coordInfoVecs.triangles2Nodes_3[i]];
									//std::cout<<ccc<<std::endl;
									if ((aaa+bbb+ccc) >= 2){
										generalParams.triangles_in_tip[i] = 1;
										//triangles_in_upperhem.push_back(i);
									}
									//else if ((aaa+bbb+ccc)==1){
									//	generalParams.triangles_in_upperhem[i] = 0;
										//triangles_in_upperhem.push_back(i);
									//}
									else{
										generalParams.triangles_in_tip[i] = -1;
									}
								//	std::cout<<"triangle "<<i<<" "<<generalParams.triangles_in_upperhem[i]<<std::endl;		
								}
								//std::cout<<"ERROR HERE 10"<<std::endl;
								//std::vector<int> edges_in_upperhem;
							//	generalParams.edges_in_upperhem.resize(coordInfoVecs.num_edges);
								
								for (int i = 0; i < coordInfoVecs.num_edges; i++){
									if (coordInfoVecs.edges2Triangles_1[i] < 0 || coordInfoVecs.edges2Triangles_1[i] > (INT_MAX-100)){
										generalParams.edges_in_tip[i] = -1;
										continue;
									}
									else if (coordInfoVecs.edges2Triangles_2[i] < 0 || coordInfoVecs.edges2Triangles_2[i] > (INT_MAX-100)){
										generalParams.edges_in_tip[i] = -1;
										continue;
									}
									
									int aaa = generalParams.triangles_in_tip[coordInfoVecs.edges2Triangles_1[i]];//generalParams.nodes_in_upperhem[coordInfoVecs.edges2Nodes_1[i]];
									int bbb = generalParams.triangles_in_tip[coordInfoVecs.edges2Triangles_2[i]];//generalParams.nodes_in_upperhem[coordInfoVecs.edges2Nodes_2[i]];
									if (aaa == 1 && bbb == 1){
										generalParams.edges_in_tip[i] = 1;
										//generalParams.edges_in_upperhem_list.push_back(i);
									}
									else if (aaa == 1 || bbb == 1){
										generalParams.edges_in_tip[i] = 1;
									}
									else{
										generalParams.edges_in_tip[i] = -1;
									}
									
								}
							}
							else{
								generalParams.nodes_in_tip = generalParams.nodes_in_upperhem;
								generalParams.edges_in_tip = generalParams.edges_in_upperhem;
								generalParams.triangles_in_tip = generalParams.triangles_in_upperhem;
							}
							edgeswap_iteration += 1;
					}
					//std::cout<<"energy_gradient = "<<energy_gradient<<std::endl;
					//std::cout<<"end relaxation step"<<std::endl;
					 
						
			
 					/*if (edgeswap_iteration % (2*RECORD_TIME) == 0){
						if (reduce_counter*0.05 < 0.89){
							reduce_counter = reduce_counter + 1;
						linearSpringInfoVecs.spring_constant_weak = linearSpringInfoVecs.spring_constant - linearSpringInfoVecs.spring_constant*(reduce_counter*0.05);
						bendingTriangleInfoVecs.spring_constant_weak = bendingTriangleInfoVecs.spring_constant - bendingTriangleInfoVecs.spring_constant*(reduce_counter*0.05);
						areaTriangleInfoVecs.spring_constant_weak = areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant*(reduce_counter*0.05);
						std::cout<<"current weakened linear spring coeff = "<<linearSpringInfoVecs.spring_constant_weak<<std::endl;
						std::cout<<"current weakened bending spring ceoff = "<<bendingTriangleInfoVecs.spring_constant_weak<<std::endl;
						std::cout<<"current weakened area spring coeff = "<<areaTriangleInfoVecs.spring_constant_weak<<std::endl;
						}
						else{}
					 }*/		
					//std::cout<<"edgeswap_iteration = "<<edgeswap_iteration<<std::endl;
 					if (edgeswap_iteration % 200 == 0){
						
						double r1,r2,r3,r1x, r1y, r1z, r2x, r2y, r2z, r3x, r3y, r3z, norm_r1r2, norm_r2r3, norm_r3r1, bud_s, bud_area;
						generalParams.current_bud_area = 0.0;
						for (int counter = 0; counter < coordInfoVecs.num_triangles; counter++){
							r1 = coordInfoVecs.triangles2Nodes_1[counter];
							r2 = coordInfoVecs.triangles2Nodes_2[counter];
							r3 = coordInfoVecs.triangles2Nodes_3[counter];
							if (generalParams.triangles_in_upperhem[counter] == 1 && r1 != INT_MAX && r2 != INT_MAX && r3 != INT_MAX){
								r1x = coordInfoVecs.nodeLocX[r1];
								r1y = coordInfoVecs.nodeLocY[r1];
								r1z = coordInfoVecs.nodeLocZ[r1];
								r2x = coordInfoVecs.nodeLocX[r2];
								r2y = coordInfoVecs.nodeLocY[r2];
								r2z = coordInfoVecs.nodeLocZ[r2];
								r3x = coordInfoVecs.nodeLocX[r3];
								r3y = coordInfoVecs.nodeLocY[r3];
								r3z = coordInfoVecs.nodeLocZ[r3];
								norm_r1r2 = sqrt((r2x-r1x)*(r2x-r1x) + (r2y-r1y)*(r2y-r1y) + (r2z-r1z)*(r2z-r1z));
								norm_r2r3 = sqrt((r3x-r2x)*(r3x-r2x) + (r3y-r2y)*(r3y-r2y) + (r3z-r2z)*(r3z-r2z));
								norm_r3r1 = sqrt((r3x-r1x)*(r3x-r1x) + (r3y-r1y)*(r3y-r1y) + (r3z-r1z)*(r3z-r1z));
								bud_s = (norm_r1r2 + norm_r2r3 + norm_r3r1)/2.0;
								bud_area = sqrt(bud_s*(bud_s-norm_r1r2)*(bud_s-norm_r2r3)*(bud_s-norm_r3r1));
								generalParams.current_bud_area += bud_area;
							}
						}
						 
					//if (GROWTH_COUNTER % RECORD_TIME == 0 && GROWTH_COUNTER != 0){

						for (int v = 0; v < coordInfoVecs.num_edges; v++){
							double ev1 = coordInfoVecs.edges2Nodes_1[v];
							double ev2 = coordInfoVecs.edges2Nodes_2[v];
							if (ev1 == INT_MAX || ev2 == INT_MAX){
								continue;
							}
							double ed = sqrt((coordInfoVecs.nodeLocX[ev2] - coordInfoVecs.nodeLocX[ev1])*(coordInfoVecs.nodeLocX[ev2] - coordInfoVecs.nodeLocX[ev1]) +
										(coordInfoVecs.nodeLocY[ev2] - coordInfoVecs.nodeLocY[ev1])*(coordInfoVecs.nodeLocY[ev2] - coordInfoVecs.nodeLocY[ev1]) +
										(coordInfoVecs.nodeLocZ[ev2] - coordInfoVecs.nodeLocZ[ev1])*(coordInfoVecs.nodeLocZ[ev2] - coordInfoVecs.nodeLocZ[ev1]));
							if (ed >= 2.0){
								std::cout<<"Edge over extension, possibly instability occuring"<<std::endl;
								runSim = false;
								initial_kT = -0.00000000001;
								break;
							}
						}
						generalParams.angle_per_edge.clear();
						//generalParams.angle_per_edge.resize(coordInfoVecs.num_edges);
						int j = 0;
						for (int j = 0; j < coordInfoVecs.num_edges; j++){
							if (coordInfoVecs.edges2Nodes_1[j] == INT_MAX || coordInfoVecs.edges2Nodes_2[j] == INT_MAX){
								generalParams.angle_per_edge.push_back(-INT_MAX);
								continue;
							}
							
							double T1 = coordInfoVecs.edges2Triangles_1[j];
							double T2 = coordInfoVecs.edges2Triangles_2[j];
							double T1v1x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_2[T1]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[T1]];
							double T1v1y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_2[T1]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[T1]];
							double T1v1z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_2[T1]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[T1]];
							double T1v2x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_3[T1]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[T1]];
							double T1v2y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_3[T1]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[T1]];
							double T1v2z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_3[T1]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[T1]];
							double T2v1x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_2[T2]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[T2]];
							double T2v1y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_2[T2]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[T2]];
							double T2v1z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_2[T2]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[T2]];
							double T2v2x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_3[T2]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[T2]];
							double T2v2y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_3[T2]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[T2]];
							double T2v2z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_3[T2]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[T2]];
							double N1x = T1v1y*T1v2z - T1v2y*T1v1z;
							double N1y = -(T1v1x*T1v2z - T1v2x*T1v1z);
							double N1z = T1v1x*T1v2y - T1v2x*T1v1y;
							double N2x = T2v1y*T2v2z - T2v2y*T2v1z;
							double N2y = -(T2v1x*T2v2z - T2v2x*T2v1z);
							double N2z = T2v1x*T2v2y - T2v2x*T2v1y;
							
							double nN1 = sqrt(N1x*N1x + N1y*N1y + N1z*N1z);
							double nN2 = sqrt(N2x*N2x + N2y*N2y + N2z*N2z);
							double cosAngle = (N1x*N2x + N1y*N2y + N1z*N2z)/ (nN1*nN2);

							double direction_check_x = N1y*N2z - N2y*N1z;
							double direction_check_y = -(N1x*N2z - N2x*N1z);
							double direction_check_z = (N1x*N2y - N2x*N1y);
							double edge_direction_x = coordInfoVecs.nodeLocX[coordInfoVecs.edges2Nodes_2[j]] - coordInfoVecs.nodeLocX[coordInfoVecs.edges2Nodes_1[j]];
							double edge_direction_y = coordInfoVecs.nodeLocY[coordInfoVecs.edges2Nodes_2[j]] - coordInfoVecs.nodeLocY[coordInfoVecs.edges2Nodes_1[j]];
							double edge_direction_z = coordInfoVecs.nodeLocZ[coordInfoVecs.edges2Nodes_2[j]] - coordInfoVecs.nodeLocZ[coordInfoVecs.edges2Nodes_1[j]];
							double direction_check = direction_check_x*edge_direction_x + direction_check_y*edge_direction_y + direction_check_z*edge_direction_z;

							if (direction_check >= 0){
								generalParams.angle_per_edge.push_back( acos(cosAngle));
							}
							else{
								generalParams.angle_per_edge.push_back( -acos(cosAngle));
							}
							//j += 1;
							
						}
						generalParams.true_num_edges = 0;
						for (int i = 0; i < coordInfoVecs.num_edges; i++){
							if (coordInfoVecs.edges2Nodes_1[i] != INT_MAX && coordInfoVecs.edges2Nodes_2[i] != INT_MAX){
								generalParams.true_num_edges += 1;
							}
						 }
						 //std::cout<<"push_nodes_down size = "<<push_nodes_down.size()<<std::endl;
						 //std::cout<<"push_nodes_up size = "<<push_nodes_up.size()<<std::endl;
						// for (int i = 0; i < push_nodes_down.size(); i++){
						 //	std::cout<<"push_nodes_down "<<push_nodes_down[i]<<std::endl;
						 //}
						 //for (int i = 0; i < push_nodes_up.size(); i++){
				//			std::cout<<"push_nodes_up "<<push_nodes_up[i]<<std::endl;
						//}
						if (true){
							max_height = -10000.0;
							min_height = 10000.0;
							for (int k = 0; k < generalParams.maxNodeCount; k++){
								if (coordInfoVecs. nodeLocZ[k] >= max_height){
									max_height = coordInfoVecs. nodeLocZ[k];
								}
								if (coordInfoVecs.nodeLocZ[k] <= min_height){
									min_height = coordInfoVecs.nodeLocZ[k];
								}
							}
							std::cout<<"cell diameter = "<<max_height - min_height<<std::endl;
						}
						 storage->print_VTK_File();
						 //storage->storeVariables();
						//  std::cout<<"current Hill equation constant = "<<generalParams.hilleqnconst<<std::endl;
						//  //storage->storeVariables();
						//  std::cout<<"current total energy = "<< new_total_energy<<std::endl;
						//  std::cout<<"LINEAR ENERGY = "<<linearSpringInfoVecs.linear_spring_energy<<std::endl;
						// std::cout<<"BEND ENERGY = "<<bendingTriangleInfoVecs.bending_triangle_energy<<std::endl;
						// std::cout<<"AREA ENERGY = "<<areaTriangleInfoVecs.area_triangle_energy<<std::endl;
						// std::cout<<"REPULSION ENERGY = "<<energy_rep<<std::endl;
						// std::cout<<"VOLUME ENERGY = "<<generalParams.volume_energy<<std::endl;
						 std::cout<<"energy_gradient = "<<energy_gradient<<std::endl;
						 std::cout<<"true current total volume = "<<generalParams.true_current_total_volume<<std::endl;
						//std::cout<<"equilibrium total volume = "<<generalParams.eq_total_volume<<std::endl;
						std::cout<<"true_current_bud_area = "<<generalParams.current_bud_area<<std::endl;
						std::cout<<"current iterative steps: "<<edgeswap_iteration<<" x Max_Runtime"<<std::endl;
					}

					if (edgeswap_iteration >= NKBT){
 					//if (edgeswap_iteration == NKBT){
						generalParams.true_num_edges = 0;
						for (int i = 0; i < coordInfoVecs.num_edges; i++){
							if (coordInfoVecs.edges2Nodes_1[i] != INT_MAX && coordInfoVecs.edges2Nodes_2[i] != INT_MAX){
								generalParams.true_num_edges += 1;
							}
						}
						storage-> print_VTK_File();
						//storage-> storeVariables();
 						generalParams.kT = -1.0;//generalParams.kT - 0.072;
						 std::cout<<"Current kBT = "<<generalParams.kT<<std::endl;
						 std::cout<<"MAX Runtime reached"<<std::endl;
 						edgeswap_iteration = 0;
 					}
 					if (generalParams.kT < min_kT){
 						initial_kT = -1.0;
					runSim = false;
					break;
					 }

//std::cout<<"ERROR BEFORE GROWTH"<<std::endl;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// GROWTH OF THE CELL (MEMBRANE) ////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
//GROWTH_COUNTER = 0;
ComputeVolume(
	generalParams,
	coordInfoVecs,
	linearSpringInfoVecs,
	ljInfoVecs);
int DYN_GROWTH_FREQUENCY = GROWTH_FREQUENCY;
// if (generalParams.true_current_total_volume/initial_volume <= 1.15){DYN_GROWTH_FREQUENCY = GROWTH_FREQUENCY;}
// else if (generalParams.true_current_total_volume/initial_volume <= 1.30){DYN_GROWTH_FREQUENCY = 2*GROWTH_FREQUENCY;}
// else if (generalParams.true_current_total_volume/initial_volume <= 1.45){DYN_GROWTH_FREQUENCY = 4*GROWTH_FREQUENCY;}
// else if (generalParams.true_current_total_volume/initial_volume <= 1.60){DYN_GROWTH_FREQUENCY = 6*GROWTH_FREQUENCY;}
// else if (generalParams.true_current_total_volume/initial_volume <= 1.75){DYN_GROWTH_FREQUENCY = 8*GROWTH_FREQUENCY;}
// else{}

if (edgeswap_iteration % DYN_GROWTH_FREQUENCY == 0){
//if (edgeswap_iteration % GROWTH_FREQUENCY == 0){
//if (1 > 0){
	GROWTH_COUNTER += 1;
	//std::cout<<"initiaze growth algorithm"<<std::endl;
	
// 	double delayed_stiffening = 1.128;
// 	if (edgeswap_iteration <= GROWTH_FREQUENCY){
// 		std::cout<<"delayed_stiffening = "<<delayed_stiffening<<std::endl;
// 		std::cout<<"This determines when stiffening actually start taking effect after the cell volume gain (X*100)% increase"<<std::endl;
// 	}
// // linearSpringInfoVecs.spring_constant_weak = linearSpringInfoVecs.spring_constant_weak +
// // 								((linearSpringInfoVecs.spring_constant - scale_linear)/500.0);//(MAX_GROWTH_EVENT_NUMBER));
// if (generalParams.true_current_total_volume >= initial_volume*delayed_stiffening){
// 	linearSpringInfoVecs.spring_constant_weak = scale_linear +
// 									((linearSpringInfoVecs.spring_constant - scale_linear)*(generalParams.true_current_total_volume - initial_volume*delayed_stiffening)/(initial_volume*MAX_VOLUME_RATIO_stiffening - initial_volume*delayed_stiffening));

// 	if (linearSpringInfoVecs.spring_constant_weak >= linearSpringInfoVecs.spring_constant){
// 		linearSpringInfoVecs.spring_constant_weak = linearSpringInfoVecs.spring_constant;
// 	}
// }

// // bendingTriangleInfoVecs.spring_constant_weak = bendingTriangleInfoVecs.spring_constant_weak +
// // 								((bendingTriangleInfoVecs.spring_constant - scale_bend)/500.0);///MAX_GROWTH_EVENT_NUMBER);;

// if (generalParams.true_current_total_volume >= initial_volume*delayed_stiffening){
// 	bendingTriangleInfoVecs.spring_constant_weak = scale_bend +
// 									((bendingTriangleInfoVecs.spring_constant - scale_bend)*(generalParams.true_current_total_volume - initial_volume*delayed_stiffening)/(initial_volume*MAX_VOLUME_RATIO_stiffening - initial_volume*delayed_stiffening));///MAX_GROWTH_EVENT_NUMBER);;

// 	if (bendingTriangleInfoVecs.spring_constant_weak >= bendingTriangleInfoVecs.spring_constant){
// 		bendingTriangleInfoVecs.spring_constant_weak = bendingTriangleInfoVecs.spring_constant;
// 	}
// }
// // areaTriangleInfoVecs.spring_constant_weak = areaTriangleInfoVecs.spring_constant_weak +
// // 								((areaTriangleInfoVecs.spring_constant - scale_area)/500.0);///MAX_GROWTH_EVENT_NUMBER);;

// if (generalParams.true_current_total_volume >= initial_volume*delayed_stiffening){
// 	areaTriangleInfoVecs.spring_constant_weak = scale_area +
// 									((areaTriangleInfoVecs.spring_constant - scale_area)*(generalParams.true_current_total_volume - initial_volume*delayed_stiffening)/(initial_volume*MAX_VOLUME_RATIO_stiffening - initial_volume*delayed_stiffening));///MAX_GROWTH_EVENT_NUMBER);;

// 	if (areaTriangleInfoVecs.spring_constant_weak >= areaTriangleInfoVecs.spring_constant){
// 		areaTriangleInfoVecs.spring_constant_weak = areaTriangleInfoVecs.spring_constant;
// 	}
// }
// // generalParams.kT = generalParams.kT + (0.07 - 0.035)/500.0;//MAX_GROWTH_EVENT_NUMBER;
// if (generalParams.true_current_total_volume >= initial_volume*delayed_stiffening){
// 	generalParams.kT = initial_kT + (initial_kT*2.0 - initial_kT)*((generalParams.true_current_total_volume - initial_volume*delayed_stiffening)/(initial_volume*MAX_VOLUME_RATIO_stiffening - initial_volume*delayed_stiffening));
// }
// std::cout<<"kT = "<<generalParams.kT<<std::endl;
// std::cout<<"linear spring coeff weak = "<<linearSpringInfoVecs.spring_constant_weak<<std::endl;
// std::cout<<"bending spring coeff weak = "<<bendingTriangleInfoVecs.spring_constant_weak<<std::endl;
// std::cout<<"area spring coeff weak = "<<areaTriangleInfoVecs.spring_constant_weak<<std::endl;
// std::cout<<"true current total volume = "<<generalParams.true_current_total_volume<<std::endl;
// std::cout<<"initial volume = "<<initial_volume<<std::endl;
// std::cout<<"initial volume * MAX_VOLUME RATIO = "<<initial_volume*MAX_VOLUME_RATIO<<std::endl;
// std::cout<<"current volume ratio = "<<(generalParams.true_current_total_volume)/(initial_volume*MAX_VOLUME_RATIO)<<std::endl;
// if (generalParams.true_current_total_volume >= initial_volume*delayed_stiffening){
// 	generalParams.strain_threshold = 0.05 +
// 								((0.4 - 0.05)*pow((generalParams.true_current_total_volume - initial_volume*delayed_stiffening)/(initial_volume*MAX_VOLUME_RATIO_stiffening - initial_volume*delayed_stiffening), 1.0));///MAX_GROWTH_EVENT_NUMBER)
// }


VectorShuffleForGrowthLoop.clear();
int VectorShuffleForGrowthLoop_COUNT = 0;
for (int y = 0; y < coordInfoVecs.num_edges; y++){
	if (generalParams.edges_in_upperhem_list[y] >= 0 &&
		generalParams.edges_in_upperhem_list[y] != INT_MAX &&
		generalParams.boundaries_in_upperhem[y] != 1){
		VectorShuffleForGrowthLoop.push_back(y);
		VectorShuffleForGrowthLoop_COUNT += 1;
	}
	/*if (generalParams.edges_in_upperhem_list[y] >= 0 &&
		generalParams.edges_in_upperhem_list[y] != INT_MAX &&
		generalParams.boundaries_in_upperhem[y] != 1 &&
		edges_in_growth[y] == 1){
		VectorShuffleForGrowthLoop.push_back(y);
	}*/
	
	
}
//std::cout<<VectorShuffleForGrowthLoop_COUNT<<std::endl;
//double budcenter_z = 0.0;
//int budcenter_z_count = 0;
max_height = -10000.0;
for (int i = 0; i < generalParams.maxNodeCount; i++){
	if (generalParams.nodes_in_upperhem[i] == 1){
		//budcenter_z += coordInfoVecs.nodeLocZ[i];
		//budcenter_z_count += 1;
		if (coordInfoVecs.nodeLocZ[i] >= max_height){
			max_height = coordInfoVecs.nodeLocZ[i];
		}
	}
}
//budcenter_z = budcenter_z/budcenter_z_count;

//std::random_device rand_dev;
//std::mt19937 generator2(rand_dev());
std::shuffle(std::begin(VectorShuffleForGrowthLoop), std::end(VectorShuffleForGrowthLoop), generator2);
//MAX_GROWTH_TEST = 4;//VectorShuffleForGrowthLoop.size();
//std::cout<<"number of growth attempt for each growth event = "<<MAX_GROWTH_TEST<<std::endl;
bool triggered = false;
int true_DELTA = 0;
double edge_midpt_z;
//std::cout<<"BEGIN GROWTH ALGORITHM"<<std::endl;
edgeswap_ptr->transferDtoH(generalParams, coordInfoVecs, build_ptr->hostSetInfoVecs);
int GROWTH_COUNT = 0;
for (int p = 0; p < VectorShuffleForGrowthLoop.size(); p++){
	if (coordInfoVecs.edges2Nodes_1[VectorShuffleForGrowthLoop[p]] < 0 || coordInfoVecs.edges2Nodes_1[VectorShuffleForGrowthLoop[p]] == INT_MAX){
		continue;
	}
	else if (coordInfoVecs.edges2Nodes_2[VectorShuffleForGrowthLoop[p]] < 0 || coordInfoVecs.edges2Nodes_2[VectorShuffleForGrowthLoop[p]] == INT_MAX){
		continue;
	}

	//edge_midpt_z = (coordInfoVecs.nodeLocZ[coordInfoVecs.edges2Nodes_1[VectorShuffleForGrowthLoop[p]]] + coordInfoVecs.nodeLocZ[coordInfoVecs.edges2Nodes_2[VectorShuffleForGrowthLoop[p]]])/2.0;
	// The commented out portion below this line restrict the growth to only happen at the tip region.
	// if (edge_midpt_z < (max_height - 1.5*generalParams.Rmin)){
	// 	continue;
	// }
	//std::cout<<"begin growth test"<<std::endl;
	int DELTA = edgeswap_ptr->growth_host_vecs(
		VectorShuffleForGrowthLoop[p],
		generalParams,
		build_ptr->hostSetInfoVecs,
		coordInfoVecs,
		linearSpringInfoVecs,
		bendingTriangleInfoVecs,
		areaTriangleInfoVecs);
	GROWTH_COUNT += DELTA;
	if (GROWTH_COUNT >= MAX_GROWTH_TEST){
		break;
	}
	//IN THIS CODE, THE GROWTH IS DETERMINISTIC SUCH THAT A CHOSEN EDGE WILL ALWAYS UNDERGO GROWTH!!!!!!!!!!!!!!!!!!!!!
	// std::cout<<"chosen edge mid point x = "<<(coordInfoVecs.nodeLocX[coordInfoVecs.edges2Nodes_1[VectorShuffleForGrowthLoop[p]]] +
	// coordInfoVecs.nodeLocX[coordInfoVecs.edges2Nodes_2[VectorShuffleForGrowthLoop[p]]])/2.0<<std::endl;
	// std::cout<<"chosen edge mid point y = "<<(coordInfoVecs.nodeLocY[coordInfoVecs.edges2Nodes_1[VectorShuffleForGrowthLoop[p]]] +
	// coordInfoVecs.nodeLocY[coordInfoVecs.edges2Nodes_2[VectorShuffleForGrowthLoop[p]]])/2.0<<std::endl;
	// std::cout<<"chosen edge mid point z = "<<(coordInfoVecs.nodeLocZ[coordInfoVecs.edges2Nodes_1[VectorShuffleForGrowthLoop[p]]] +
	// coordInfoVecs.nodeLocZ[coordInfoVecs.edges2Nodes_2[VectorShuffleForGrowthLoop[p]]])/2.0<<std::endl;
}
edgeswap_ptr->transferHtoD(generalParams, coordInfoVecs, build_ptr->hostSetInfoVecs);
//std::cout<<"END GROWTH ALGORITHM"<<std::endl;
if (GROWTH_COUNT >= 1){

	triggered = true;
	max_height = -10000.0;
	current_center_x = 0.0;
	current_center_y = 0.0;
	for (int k = 0; k < generalParams.maxNodeCount; k++){
		if (generalParams.nodes_in_upperhem[k] < (INT_MAX-100)){
			current_center_x += coordInfoVecs.nodeLocX[k];
			current_center_y += coordInfoVecs.nodeLocY[k];
		}
	}
	current_center_x = current_center_x/generalParams.maxNodeCount;
	current_center_y = current_center_y/generalParams.maxNodeCount;

	generalParams.septin_ring_z = 0.0;
	generalParams.boundary_z = 0.0;
	for (int k = 0; k < boundary_edge_list.size(); k++){
	//for (int k = 0; k < boundary_node_list.size(); k++){
		double n1 = coordInfoVecs.edges2Nodes_1[boundary_edge_list[k]];
		double n2 = coordInfoVecs.edges2Nodes_2[boundary_edge_list[k]];
		double cent_of_edge_x = (coordInfoVecs.nodeLocX[n1] + coordInfoVecs.nodeLocX[n2])/2.0;
		double cent_of_edge_y = (coordInfoVecs.nodeLocY[n1] + coordInfoVecs.nodeLocY[n2])/2.0;
		double cent_of_edge_z = (coordInfoVecs.nodeLocZ[n1] + coordInfoVecs.nodeLocZ[n2])/2.0;
		double dist_x = coordInfoVecs.nodeLocX[max_height_index] - cent_of_edge_x;
		double dist_y = coordInfoVecs.nodeLocY[max_height_index] - cent_of_edge_y;
		double dist_z = coordInfoVecs.nodeLocZ[max_height_index] - cent_of_edge_z;
		// double temp_dist = sqrt((coordInfoVecs.nodeLocX[max_height_index] - coordInfoVecs.nodeLocX[n1])*(coordInfoVecs.nodeLocX[max_height_index] - coordInfoVecs.nodeLocX[n1]) +
		// (coordInfoVecs.nodeLocY[max_height_index] - coordInfoVecs.nodeLocY[n1])*(coordInfoVecs.nodeLocY[max_height_index] - coordInfoVecs.nodeLocY[n1]) +
		// 	(coordInfoVecs.nodeLocZ[max_height_index] - coordInfoVecs.nodeLocZ[n1])*(coordInfoVecs.nodeLocZ[max_height_index] - coordInfoVecs.nodeLocZ[n1]));
		double temp_dist = sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);
		generalParams.septin_ring_z += (coordInfoVecs.nodeLocZ[n1] + coordInfoVecs.nodeLocZ[n2]);
		if (temp_dist >= dtb){
			dtb = temp_dist;
			/* "dtb" will be used to identify where the septin ring is located, and used to determine the Hill coefficient*/
		}
	}

	
	generalParams.septin_ring_z = generalParams.septin_ring_z/(1.0*boundary_node_list.size());
	generalParams.boundary_z = generalParams.septin_ring_z;
	edgeswap_ptr->transferDtoH(generalParams, coordInfoVecs, build_ptr->hostSetInfoVecs);//Currently this is treated as a backup of coordInfoVecs
	edgeswap_ptr->gradient_weakening_update_host_vecs_tip(sigma,
		current_center_x,
		current_center_y,
		max_height,
		dtb,
		3.32*generalParams.Rmin,
		generalParams,
		coordInfoVecs,
		build_ptr->hostSetInfoVecs);
	edgeswap_ptr->transferHtoD(generalParams, coordInfoVecs, build_ptr->hostSetInfoVecs);

	for (int k = 0; k < generalParams.maxNodeCount; k++){								
		if (coordInfoVecs. nodeLocZ[k] >= max_height){
			max_height = coordInfoVecs.nodeLocZ[k];
			max_height_index = k;
		}
		
	}
	if (bud_tip_necessary == true){
		double current_dist;
		//std::cout<<"ERROR HERE 8"<<std::endl;
		for (int i = 0; i < generalParams.maxNodeCount; i++){
			current_dist = sqrt((coordInfoVecs.nodeLocX[i] - current_center_x)*(coordInfoVecs.nodeLocX[i] - current_center_x) +
			(coordInfoVecs.nodeLocY[i] - current_center_y)*(coordInfoVecs.nodeLocY[i] - current_center_y) +
			(coordInfoVecs.nodeLocZ[i] - coordInfoVecs.nodeLocZ[max_height_index])*(coordInfoVecs.nodeLocZ[i] - coordInfoVecs.nodeLocZ[max_height_index]));
			//if (coordInfoVecs.nodeLocZ[i] > (coordInfoVecs.nodeLocZ[max_height_index] - 1.5*generalParams.Rmin)){
			if (current_dist <= 1.0){//2.5*generalParams.Rmin){
				generalParams.nodes_in_tip[i] = 1;
			}
			else{
				generalParams.nodes_in_tip[i] = -1;
			}
			//std::cout<<"nodes "<<i<<" "<<generalParams.nodes_in_upperhem[i]<<std::endl;		
		}
		//std::cout<<"ERROR HERE 9"<<std::endl;
		for (int i = 0; i < coordInfoVecs.num_triangles; i++){
			//std::cout<<"i = "<<i<<std::endl;
			//std::cout<<coordInfoVecs.triangles2Nodes_1[i]<<" "<<coordInfoVecs.triangles2Nodes_2[i]<<" "<<coordInfoVecs.triangles2Nodes_3[i]<<std::endl;
			if (coordInfoVecs.triangles2Nodes_1[i] < 0 || coordInfoVecs.triangles2Nodes_1[i] > (INT_MAX-100)){
				generalParams.triangles_in_tip[i] = -1;
				continue;
			}
			else if (coordInfoVecs.triangles2Nodes_2[i] < 0 || coordInfoVecs.triangles2Nodes_2[i] > (INT_MAX-100)){
				generalParams.triangles_in_tip[i] = -1;
				continue;
			}
			else if (coordInfoVecs.triangles2Nodes_3[i] < 0 || coordInfoVecs.triangles2Nodes_3[i] > (INT_MAX-100)){
				generalParams.triangles_in_tip[i] = -1;
				continue;
			}

			int aaa = generalParams.nodes_in_tip[coordInfoVecs.triangles2Nodes_1[i]];
			//std::cout<<aaa<<std::endl;
			int bbb = generalParams.nodes_in_tip[coordInfoVecs.triangles2Nodes_2[i]];
			//std::cout<<bbb<<std::endl;
			int ccc = generalParams.nodes_in_tip[coordInfoVecs.triangles2Nodes_3[i]];
			//std::cout<<ccc<<std::endl;
			if ((aaa+bbb+ccc) >= 2){
				generalParams.triangles_in_tip[i] = 1;
				//triangles_in_upperhem.push_back(i);
			}
			//else if ((aaa+bbb+ccc)==1){
			//	generalParams.triangles_in_upperhem[i] = 0;
				//triangles_in_upperhem.push_back(i);
			//}
			else{
				generalParams.triangles_in_tip[i] = -1;
			}
		//	std::cout<<"triangle "<<i<<" "<<generalParams.triangles_in_upperhem[i]<<std::endl;		
		}
		//std::cout<<"ERROR HERE 10"<<std::endl;
		//std::vector<int> edges_in_upperhem;
	//	generalParams.edges_in_upperhem.resize(coordInfoVecs.num_edges);
		
		for (int i = 0; i < coordInfoVecs.num_edges; i++){
			if (coordInfoVecs.edges2Triangles_1[i] < 0 || coordInfoVecs.edges2Triangles_1[i] > (INT_MAX-100)){
				generalParams.edges_in_tip[i] = -1;
				continue;
			}
			else if (coordInfoVecs.edges2Triangles_2[i] < 0 || coordInfoVecs.edges2Triangles_2[i] > (INT_MAX-100)){
				generalParams.edges_in_tip[i] = -1;
				continue;
			}
			
			int aaa = generalParams.triangles_in_tip[coordInfoVecs.edges2Triangles_1[i]];//generalParams.nodes_in_upperhem[coordInfoVecs.edges2Nodes_1[i]];
			int bbb = generalParams.triangles_in_tip[coordInfoVecs.edges2Triangles_2[i]];//generalParams.nodes_in_upperhem[coordInfoVecs.edges2Nodes_2[i]];
			if (aaa == 1 && bbb == 1){
				generalParams.edges_in_tip[i] = 1;
				//generalParams.edges_in_upperhem_list.push_back(i);
			}
			else if (aaa == 1 || bbb == 1){
				generalParams.edges_in_tip[i] = 1;
			}
			else{
				generalParams.edges_in_tip[i] = -1;
			}
			
		}
	}
	else{
		generalParams.nodes_in_tip = generalParams.nodes_in_upperhem;
		generalParams.edges_in_tip = generalParams.edges_in_upperhem;
		generalParams.triangles_in_tip = generalParams.triangles_in_upperhem;
	}
}

				if (triggered == true){	
					true_num_edges_in_upperhem = 0;
					for (int i = 0; i < coordInfoVecs.num_edges; i++){
						if (generalParams.edges_in_upperhem_list[i] != INT_MAX && generalParams.edges_in_upperhem_list[i] >= 0){
							true_num_edges_in_upperhem += 1;
							
						}
					}
					//std::cout<<"WHERE iS THE PROBLEM 3"<<std::endl;
				}
			}
			
			
			
			
//std::cout<<"GROWTH DONE!"<<std::endl;
 ////storage->print_VTK_File();
////storage->storeVariables();

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////// END OF GROWTH SECTION //////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					

ComputeVolume(
	generalParams,
	coordInfoVecs,
	linearSpringInfoVecs,
	ljInfoVecs);
					
					
 			}
		
		}
		

	};
	
	





void System::assignStorage(std::shared_ptr<Storage> _storage) {
	storage = _storage;
};
void System::set_weak_builder(std::weak_ptr<SystemBuilder> _weak_bld_ptr) {
	weak_bld_ptr = _weak_bld_ptr;
};



//initialize memory for thrust vectors and set coordInfoVecs vals from input. 
void System::initializeSystem(HostSetInfoVecs& hostSetInfoVecs) {
	std::cout<<"Initializing"<<std::endl;

	generalParams.maxNodeCount = hostSetInfoVecs.nodeLocX.size();
	coordInfoVecs.num_edges = hostSetInfoVecs.edges2Nodes_1.size();
	coordInfoVecs.num_triangles = hostSetInfoVecs.triangles2Nodes_1.size();

	std::cout<<"num nodes: "<< generalParams.maxNodeCount << std::endl;
	std::cout<<"num edges: "<< coordInfoVecs.num_edges << std::endl;
	std::cout<<"num elems: "<< coordInfoVecs.num_triangles << std::endl;
	//allocate memory
	int mem_prealloc = 3;
	coordInfoVecs.scaling_per_edge.resize(mem_prealloc*coordInfoVecs.num_edges, 0.0);
	hostSetInfoVecs.scaling_per_edge.resize(coordInfoVecs.scaling_per_edge.size(), 0.0);

	coordInfoVecs.isNodeFixed.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size(),false);
	coordInfoVecs.prevNodeLocX.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size());
	coordInfoVecs.prevNodeLocY.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size());
	coordInfoVecs.prevNodeLocZ.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size());

	coordInfoVecs.prevNodeForceX.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size());
	coordInfoVecs.prevNodeForceY.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size());
	coordInfoVecs.prevNodeForceZ.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size());
	
	coordInfoVecs.nodeLocX.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size());
	coordInfoVecs.nodeLocY.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size());
	coordInfoVecs.nodeLocZ.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size());

	coordInfoVecs.nodeForceX.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size(), 0.0);
	coordInfoVecs.nodeForceY.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size(), 0.0);
	coordInfoVecs.nodeForceZ.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size(), 0.0);

	coordInfoVecs.triangles2Nodes_1.resize( mem_prealloc*coordInfoVecs.num_triangles );
	coordInfoVecs.triangles2Nodes_2.resize( mem_prealloc*coordInfoVecs.num_triangles );
	coordInfoVecs.triangles2Nodes_3.resize( mem_prealloc*coordInfoVecs.num_triangles );
	
	coordInfoVecs.triangles2Edges_1.resize( mem_prealloc*coordInfoVecs.num_triangles );
	coordInfoVecs.triangles2Edges_2.resize( mem_prealloc*coordInfoVecs.num_triangles );
	coordInfoVecs.triangles2Edges_3.resize( mem_prealloc*coordInfoVecs.num_triangles );

	coordInfoVecs.edges2Nodes_1.resize( mem_prealloc*coordInfoVecs.num_edges );
	coordInfoVecs.edges2Nodes_2.resize( mem_prealloc*coordInfoVecs.num_edges );
	
	coordInfoVecs.edges2Triangles_1.resize( mem_prealloc*coordInfoVecs.num_edges );
	coordInfoVecs.edges2Triangles_2.resize( mem_prealloc*coordInfoVecs.num_edges );

	coordInfoVecs.nndata1.resize( mem_prealloc*generalParams.maxNodeCount);
	coordInfoVecs.nndata2.resize( mem_prealloc*generalParams.maxNodeCount);
	coordInfoVecs.nndata3.resize( mem_prealloc*generalParams.maxNodeCount);
	coordInfoVecs.nndata4.resize( mem_prealloc*generalParams.maxNodeCount);
	coordInfoVecs.nndata5.resize( mem_prealloc*generalParams.maxNodeCount);
	coordInfoVecs.nndata6.resize( mem_prealloc*generalParams.maxNodeCount);
	coordInfoVecs.nndata7.resize( mem_prealloc*generalParams.maxNodeCount);
	coordInfoVecs.nndata8.resize( mem_prealloc*generalParams.maxNodeCount);
	coordInfoVecs.nndata9.resize( mem_prealloc*generalParams.maxNodeCount);
	//coordInfoVecs.nndata10.resize( mem_prealloc*generalParams.maxNodeCount);
	//coordInfoVecs.nndata11.resize( mem_prealloc*generalParams.maxNodeCount);
	//coordInfoVecs.nndata12.resize( mem_prealloc*generalParams.maxNodeCount);

	coordInfoVecs.SurfaceNormalX.resize( mem_prealloc*generalParams.maxNodeCount);
	coordInfoVecs.SurfaceNormalY.resize( mem_prealloc*generalParams.maxNodeCount);
	coordInfoVecs.SurfaceNormalZ.resize( mem_prealloc*generalParams.maxNodeCount);

	generalParams.nodes_in_upperhem.resize(mem_prealloc*generalParams.maxNodeCount);
	generalParams.triangles_in_upperhem.resize(mem_prealloc*coordInfoVecs.num_triangles);
	generalParams.edges_in_upperhem.resize(mem_prealloc*coordInfoVecs.num_edges);
	generalParams.edges_in_upperhem_list.resize(mem_prealloc*coordInfoVecs.num_edges);
	generalParams.boundaries_in_upperhem.resize(mem_prealloc*coordInfoVecs.num_edges, -1);
	generalParams.nodes_in_tip.resize(mem_prealloc*generalParams.maxNodeCount);
	generalParams.edges_in_tip.resize(mem_prealloc*coordInfoVecs.num_edges);
	generalParams.triangles_in_tip.resize(mem_prealloc*coordInfoVecs.num_triangles);

	hostSetInfoVecs.nodes_in_upperhem.resize(generalParams.nodes_in_upperhem.size());
	hostSetInfoVecs.triangles_in_upperhem.resize(generalParams.triangles_in_upperhem.size());
	hostSetInfoVecs.edges_in_upperhem.resize(generalParams.edges_in_upperhem.size());
	hostSetInfoVecs.edges_in_upperhem_list.resize(mem_prealloc*coordInfoVecs.num_edges);
	hostSetInfoVecs.boundaries_in_upperhem.resize(mem_prealloc*coordInfoVecs.num_edges, -1);
	hostSetInfoVecs.nodes_in_tip.resize(generalParams.nodes_in_tip.size());
	hostSetInfoVecs.edges_in_tip.resize(generalParams.edges_in_tip.size());
	hostSetInfoVecs.triangles_in_tip.resize(generalParams.triangles_in_tip.size());

	hostSetInfoVecs.nodes2Triangles_1.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	hostSetInfoVecs.nodes2Triangles_2.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	hostSetInfoVecs.nodes2Triangles_3.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	hostSetInfoVecs.nodes2Triangles_4.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	hostSetInfoVecs.nodes2Triangles_5.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	hostSetInfoVecs.nodes2Triangles_6.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	hostSetInfoVecs.nodes2Triangles_7.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	hostSetInfoVecs.nodes2Triangles_8.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	hostSetInfoVecs.nodes2Triangles_9.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	

	coordInfoVecs.nodes2Triangles_1.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	coordInfoVecs.nodes2Triangles_2.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	coordInfoVecs.nodes2Triangles_3.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	coordInfoVecs.nodes2Triangles_4.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	coordInfoVecs.nodes2Triangles_5.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	coordInfoVecs.nodes2Triangles_6.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	coordInfoVecs.nodes2Triangles_7.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	coordInfoVecs.nodes2Triangles_8.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	coordInfoVecs.nodes2Triangles_9.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	

	thrust::copy(coordInfoVecs.nodes2Triangles_1.begin(), coordInfoVecs.nodes2Triangles_1.end(), hostSetInfoVecs.nodes2Triangles_1.begin() );
	thrust::copy(coordInfoVecs.nodes2Triangles_2.begin(), coordInfoVecs.nodes2Triangles_2.end(), hostSetInfoVecs.nodes2Triangles_2.begin() );
	thrust::copy(coordInfoVecs.nodes2Triangles_3.begin(), coordInfoVecs.nodes2Triangles_3.end(), hostSetInfoVecs.nodes2Triangles_3.begin() );
	thrust::copy(coordInfoVecs.nodes2Triangles_4.begin(), coordInfoVecs.nodes2Triangles_4.end(), hostSetInfoVecs.nodes2Triangles_4.begin() );
	thrust::copy(coordInfoVecs.nodes2Triangles_5.begin(), coordInfoVecs.nodes2Triangles_5.end(), hostSetInfoVecs.nodes2Triangles_5.begin() );
	thrust::copy(coordInfoVecs.nodes2Triangles_6.begin(), coordInfoVecs.nodes2Triangles_6.end(), hostSetInfoVecs.nodes2Triangles_6.begin() );
	thrust::copy(coordInfoVecs.nodes2Triangles_7.begin(), coordInfoVecs.nodes2Triangles_7.end(), hostSetInfoVecs.nodes2Triangles_7.begin() );
	thrust::copy(coordInfoVecs.nodes2Triangles_8.begin(), coordInfoVecs.nodes2Triangles_8.end(), hostSetInfoVecs.nodes2Triangles_8.begin() );
	thrust::copy(coordInfoVecs.nodes2Triangles_9.begin(), coordInfoVecs.nodes2Triangles_9.end(), hostSetInfoVecs.nodes2Triangles_9.begin() );
	//thrust::copy(coordInfoVecs.nodes2Triangles_10.begin(), coordInfoVecs.nodes2Triangles_10.end(), hostInfoVecs.nodes2Triangles_10.begin() );
	//thrust::copy(coordInfoVecs.nodes2Triangles_11.begin(), coordInfoVecs.nodes2Triangles_11.end(), hostInfoVecs.nodes2Triangles_11.begin() );
	//thrust::copy(coordInfoVecs.nodes2Triangles_12.begin(), coordInfoVecs.nodes2Triangles_12.end(), hostInfoVecs.nodes2Triangles_12.begin() );

	//copy info to GPU
	std::cout<<"Copying"<<std::endl;
	thrust::copy(hostSetInfoVecs.isNodeFixed.begin(),hostSetInfoVecs.isNodeFixed.end(), coordInfoVecs.isNodeFixed.begin());
	
	std::cout<<"fixed_node_in_host: "<<std::endl;
	for (int k = 0; k < hostSetInfoVecs.isNodeFixed.size(); k++){
		//std::cout<<hostSetInfoVecs.isNodeFixed[k]<<std::endl;
	}
	std::cout<<"end_of_fixed_node_host_printout"<<std::endl;
	std::cout<<"fixed_node_in_device: "<<std::endl;
	for (int k = 0; k < coordInfoVecs.isNodeFixed.size(); k++){
		//std::cout<<coordInfoVecs.isNodeFixed[k]<<std::endl;
	}
	std::cout<<"end_of_fixed_node_device_printout"<<std::endl;
std::cout<<"size of host fixed "<< hostSetInfoVecs.isNodeFixed.size()<<std::endl;
std::cout<<"size of device fixed "<< coordInfoVecs.isNodeFixed.size()<<std::endl;

	/*for (int k = 0; k < coordInfoVecs.isNodeFixed.size(); k++){
		bool isFixedHost = hostSetInfoVecs.isNodeFixed[k];
		bool isFixedDevice = coordInfoVecs.isNodeFixed[k];
		if (isFixedDevice != isFixedHost){

			std::cout<<"pos "<< k << " dev val = " << coordInfoVecs.isNodeFixed[k]
				<< " host val = " <<  hostSetInfoVecs.isNodeFixed[k] <<std::endl;
		}
	}*/
	thrust::fill(coordInfoVecs.nodeForceX.begin(), coordInfoVecs.nodeForceX.end(), 0.0);
	thrust::fill(coordInfoVecs.nodeForceY.begin(), coordInfoVecs.nodeForceY.end(), 0.0);
	thrust::fill(coordInfoVecs.nodeForceZ.begin(), coordInfoVecs.nodeForceZ.end(), 0.0);

	thrust::fill(coordInfoVecs.prevNodeForceX.begin(), coordInfoVecs.prevNodeForceX.end(), 0.0);
	thrust::fill(coordInfoVecs.prevNodeForceY.begin(), coordInfoVecs.prevNodeForceY.end(), 0.0);
	thrust::fill(coordInfoVecs.prevNodeForceZ.begin(), coordInfoVecs.prevNodeForceZ.end(), 0.0);
	
	thrust::copy(hostSetInfoVecs.nodeLocX.begin(), hostSetInfoVecs.nodeLocX.end(), coordInfoVecs.prevNodeLocX.begin() );
	thrust::copy(hostSetInfoVecs.nodeLocY.begin(), hostSetInfoVecs.nodeLocY.end(), coordInfoVecs.prevNodeLocY.begin() );
	thrust::copy(hostSetInfoVecs.nodeLocZ.begin(), hostSetInfoVecs.nodeLocZ.end(), coordInfoVecs.prevNodeLocZ.begin() );
	
	thrust::copy(hostSetInfoVecs.nodeLocX.begin(), hostSetInfoVecs.nodeLocX.end(), coordInfoVecs.nodeLocX.begin() );
	thrust::copy(hostSetInfoVecs.nodeLocY.begin(), hostSetInfoVecs.nodeLocY.end(), coordInfoVecs.nodeLocY.begin() );
	thrust::copy(hostSetInfoVecs.nodeLocZ.begin(), hostSetInfoVecs.nodeLocZ.end(), coordInfoVecs.nodeLocZ.begin() );
	
	thrust::copy(hostSetInfoVecs.triangles2Nodes_1.begin(), hostSetInfoVecs.triangles2Nodes_1.end(), coordInfoVecs.triangles2Nodes_1.begin() );
	thrust::copy(hostSetInfoVecs.triangles2Nodes_2.begin(), hostSetInfoVecs.triangles2Nodes_2.end(), coordInfoVecs.triangles2Nodes_2.begin() );
	thrust::copy(hostSetInfoVecs.triangles2Nodes_3.begin(), hostSetInfoVecs.triangles2Nodes_3.end(), coordInfoVecs.triangles2Nodes_3.begin() );
	
	thrust::copy(hostSetInfoVecs.triangles2Edges_1.begin(), hostSetInfoVecs.triangles2Edges_1.end(), coordInfoVecs.triangles2Edges_1.begin() );
	thrust::copy(hostSetInfoVecs.triangles2Edges_2.begin(), hostSetInfoVecs.triangles2Edges_2.end(), coordInfoVecs.triangles2Edges_2.begin() );
	thrust::copy(hostSetInfoVecs.triangles2Edges_3.begin(), hostSetInfoVecs.triangles2Edges_3.end(), coordInfoVecs.triangles2Edges_3.begin() );

	thrust::copy(hostSetInfoVecs.edges2Nodes_1.begin(), hostSetInfoVecs.edges2Nodes_1.end(), coordInfoVecs.edges2Nodes_1.begin() );
	thrust::copy(hostSetInfoVecs.edges2Nodes_2.begin(), hostSetInfoVecs.edges2Nodes_2.end(), coordInfoVecs.edges2Nodes_2.begin() );
	
	thrust::copy(hostSetInfoVecs.edges2Triangles_1.begin(), hostSetInfoVecs.edges2Triangles_1.end(), coordInfoVecs.edges2Triangles_1.begin() );
	thrust::copy(hostSetInfoVecs.edges2Triangles_2.begin(), hostSetInfoVecs.edges2Triangles_2.end(), coordInfoVecs.edges2Triangles_2.begin() );

	thrust::copy(hostSetInfoVecs.nndata1.begin(), hostSetInfoVecs.nndata1.end(), coordInfoVecs.nndata1.begin() );
	thrust::copy(hostSetInfoVecs.nndata2.begin(), hostSetInfoVecs.nndata2.end(), coordInfoVecs.nndata2.begin() );
	thrust::copy(hostSetInfoVecs.nndata3.begin(), hostSetInfoVecs.nndata3.end(), coordInfoVecs.nndata3.begin() );
	thrust::copy(hostSetInfoVecs.nndata4.begin(), hostSetInfoVecs.nndata4.end(), coordInfoVecs.nndata4.begin() );
	thrust::copy(hostSetInfoVecs.nndata5.begin(), hostSetInfoVecs.nndata5.end(), coordInfoVecs.nndata5.begin() );
	thrust::copy(hostSetInfoVecs.nndata6.begin(), hostSetInfoVecs.nndata6.end(), coordInfoVecs.nndata6.begin() );
	thrust::copy(hostSetInfoVecs.nndata7.begin(), hostSetInfoVecs.nndata7.end(), coordInfoVecs.nndata7.begin() );
	thrust::copy(hostSetInfoVecs.nndata8.begin(), hostSetInfoVecs.nndata8.end(), coordInfoVecs.nndata8.begin() );
	thrust::copy(hostSetInfoVecs.nndata9.begin(), hostSetInfoVecs.nndata9.end(), coordInfoVecs.nndata9.begin() );
	//thrust::copy(hostSetInfoVecs.nndata10.begin(), hostSetInfoVecs.nndata10.end(), coordInfoVecs.nndata10.begin() );
	//thrust::copy(hostSetInfoVecs.nndata11.begin(), hostSetInfoVecs.nndata11.end(), coordInfoVecs.nndata11.begin() );
	//thrust::copy(hostSetInfoVecs.nndata12.begin(), hostSetInfoVecs.nndata12.end(), coordInfoVecs.nndata12.begin() );


 
	//allocate memory for other data structures.   

	//area triangle info vec
	//number of area springs is the number of triangles
	std::cout<<"Mem"<<std::endl;
	areaTriangleInfoVecs.tempNodeIdUnreduced.resize(mem_prealloc*areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceXUnreduced.resize(mem_prealloc*areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceYUnreduced.resize(mem_prealloc*areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceZUnreduced.resize(mem_prealloc*areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	
	areaTriangleInfoVecs.tempNodeIdReduced.resize(mem_prealloc*areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceXReduced.resize(mem_prealloc*areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceYReduced.resize(mem_prealloc*areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceZReduced.resize(mem_prealloc*areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);

	//beinding triangle info vec
	//num bending springs is the number of times each edge is between two triangles. 
	bendingTriangleInfoVecs.numBendingSprings = coordInfoVecs.num_edges;//coordInfoVecs.edges2Triangles_1.size();

	bendingTriangleInfoVecs.tempNodeIdUnreduced.resize(mem_prealloc*bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceXUnreduced.resize(mem_prealloc*bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceYUnreduced.resize(mem_prealloc*bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceZUnreduced.resize(mem_prealloc*bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	
	bendingTriangleInfoVecs.tempNodeIdReduced.resize(mem_prealloc*bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceXReduced.resize(mem_prealloc*bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceYReduced.resize(mem_prealloc*bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceZReduced.resize(mem_prealloc*bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);

	//linear springs
	
	linearSpringInfoVecs.tempNodeIdUnreduced.resize(mem_prealloc*linearSpringInfoVecs.factor*coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceXUnreduced.resize(mem_prealloc*linearSpringInfoVecs.factor*coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceYUnreduced.resize(mem_prealloc*linearSpringInfoVecs.factor*coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceZUnreduced.resize(mem_prealloc*linearSpringInfoVecs.factor*coordInfoVecs.num_edges);
	
	linearSpringInfoVecs.tempNodeIdReduced.resize(mem_prealloc*linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceXReduced.resize(mem_prealloc*linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceYReduced.resize(mem_prealloc*linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceZReduced.resize(mem_prealloc*linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
	
	linearSpringInfoVecs.edge_initial_length.clear();
	//linearSpringInfoVecs.edge_initial_length.resize(mem_prealloc*coordInfoVecs.num_edges,1.0);
	
	//thrust::copy(hostSetInfoVecs.edge_initial_length.begin(), hostSetInfoVecs.edge_initial_length.end(), linearSpringInfoVecs.edge_initial_length.begin() );

	//Resize the hostSetInfoVecs so that we can copy data back and forth between hostSetinfoVecs and coordInfoVecs without problem.
	hostSetInfoVecs.isNodeFixed.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size());
	//hostSetInfoVecs.prevNodeLocX.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size());
	//hostSetInfoVecs.prevNodeLocY.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size());
	//hostSetInfoVecs.prevNodeLocZ.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size());

	//hostSetInfoVecs.prevNodeForceX.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size());
	//hostSetInfoVecs.prevNodeForceY.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size());
	//hostSetInfoVecs.prevNodeForceZ.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size());
	
	hostSetInfoVecs.nodeLocX.resize(coordInfoVecs.nodeLocX.size());
	hostSetInfoVecs.nodeLocY.resize(coordInfoVecs.nodeLocX.size());
	hostSetInfoVecs.nodeLocZ.resize(coordInfoVecs.nodeLocX.size());
	std::cout<<"Host_nodeLocX size = "<<hostSetInfoVecs.nodeLocX.size()<<std::endl;

	hostSetInfoVecs.nodeForceX.resize(coordInfoVecs.nodeLocX.size());
	hostSetInfoVecs.nodeForceY.resize(coordInfoVecs.nodeLocX.size());
	hostSetInfoVecs.nodeForceZ.resize(coordInfoVecs.nodeLocX.size());
	std::cout<<"Host_nodeForceX size = "<<hostSetInfoVecs.nodeLocX.size()<<std::endl;

	hostSetInfoVecs.triangles2Nodes_1.resize( coordInfoVecs.triangles2Nodes_1.size() );
	hostSetInfoVecs.triangles2Nodes_2.resize( coordInfoVecs.triangles2Nodes_2.size() );
	hostSetInfoVecs.triangles2Nodes_3.resize( coordInfoVecs.triangles2Nodes_3.size() );
	std::cout<<"Host_triangles2Nodes size = "<<hostSetInfoVecs.triangles2Nodes_1.size()<<std::endl;
	
	hostSetInfoVecs.triangles2Edges_1.resize( coordInfoVecs.triangles2Edges_1.size() );
	hostSetInfoVecs.triangles2Edges_2.resize( coordInfoVecs.triangles2Edges_2.size() );
	hostSetInfoVecs.triangles2Edges_3.resize( coordInfoVecs.triangles2Edges_3.size() );
	std::cout<<"Host_triangles2Edges size = "<<hostSetInfoVecs.triangles2Edges_1.size()<<std::endl;

	hostSetInfoVecs.edges2Nodes_1.resize( coordInfoVecs.edges2Nodes_1.size() );
	hostSetInfoVecs.edges2Nodes_2.resize( coordInfoVecs.edges2Nodes_2.size() );
	std::cout<<"Host_edges2Nodes size = "<<hostSetInfoVecs.edges2Nodes_1.size()<<std::endl;
	
	hostSetInfoVecs.edges2Triangles_1.resize( coordInfoVecs.edges2Triangles_1.size() );
	hostSetInfoVecs.edges2Triangles_2.resize( coordInfoVecs.edges2Triangles_2.size() );
	std::cout<<"Host_edges2Triangles size = "<<hostSetInfoVecs.edges2Triangles_1.size()<<std::endl;

	hostSetInfoVecs.nndata1.resize( mem_prealloc*generalParams.maxNodeCount);
	hostSetInfoVecs.nndata2.resize( mem_prealloc*generalParams.maxNodeCount);
	hostSetInfoVecs.nndata3.resize( mem_prealloc*generalParams.maxNodeCount);
	hostSetInfoVecs.nndata4.resize( mem_prealloc*generalParams.maxNodeCount);
	hostSetInfoVecs.nndata5.resize( mem_prealloc*generalParams.maxNodeCount);
	hostSetInfoVecs.nndata6.resize( mem_prealloc*generalParams.maxNodeCount);
	hostSetInfoVecs.nndata7.resize( mem_prealloc*generalParams.maxNodeCount);
	hostSetInfoVecs.nndata8.resize( mem_prealloc*generalParams.maxNodeCount);
	hostSetInfoVecs.nndata9.resize( mem_prealloc*generalParams.maxNodeCount);
	//hostSetInfoVecs.nndata10.resize( mem_prealloc*generalParams.maxNodeCount);
	//hostSetInfoVecs.nndata11.resize( mem_prealloc*generalParams.maxNodeCount);
	//hostSetInfoVecs.nndata12.resize( mem_prealloc*generalParams.maxNodeCount);

	//std::cout<<"initial lengths: "<< linearSpringInfoVecs.edge_initial_length.size()<<std::endl;

	std::cout<<"System Ready"<<std::endl;

	//Generate LJ particle list. and set LJ particle midpoint.
	//double maxX_lj = *(thrust::max_element(coordInfoVecs.nodeLocX.begin(),coordInfoVecs.nodeLocX.end()));
	//double minX_lj = *(thrust::min_element(coordInfoVecs.nodeLocX.begin(),coordInfoVecs.nodeLocX.end()));
	//double maxY_lj = *(thrust::max_element(coordInfoVecs.nodeLocY.begin(),coordInfoVecs.nodeLocY.end()));
	//double minY_lj = *(thrust::min_element(coordInfoVecs.nodeLocY.begin(),coordInfoVecs.nodeLocY.end()));
	
	//ljInfoVecs.LJ_PosX = (maxX_lj + minX_lj)/2.0;
	//ljInfoVecs.LJ_PosY = (maxY_lj + minY_lj)/2.0;


	//currently unused
	/*thrust::host_vector<int> tempIds;
	for (int i = 0; i < hostSetInfoVecs.nodeLocX.size(); i++ ) {
		double xLoc = hostSetInfoVecs.nodeLocX[i];
		double yLoc = hostSetInfoVecs.nodeLocY[i];
		double zLoc = hostSetInfoVecs.nodeLocZ[i];
		
		double xDist = ljInfoVecs.LJ_PosX - xLoc;
		double yDist = ljInfoVecs.LJ_PosY - yLoc;
		double zDist = ljInfoVecs.LJ_PosZ - zLoc;

		double dist = std::sqrt(xDist*xDist + yDist*yDist + zDist*zDist);
		//just test all poitns for now. Optimize later.
		if (dist < ljInfoVecs.Rcutoff) {
			tempIds.push_back(i);
		}
	}
	ljInfoVecs.node_id_close.resize( tempIds.size() );
	thrust::copy(tempIds.begin(), tempIds.end(), ljInfoVecs.node_id_close.begin());
	std::cout<<"lj nodes: "<< ljInfoVecs.node_id_close.size() << std::endl;*/






	//last, set memory foor buckets.
	auxVecs.id_bucket.resize(generalParams.maxNodeCount);
	auxVecs.id_value.resize(generalParams.maxNodeCount);
	auxVecs.id_bucket_expanded.resize(27 * (generalParams.maxNodeCount));
	auxVecs.id_value_expanded.resize(27 *( generalParams.maxNodeCount ));
 


};


