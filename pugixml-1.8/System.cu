#include "System.h"
#include "SystemStructures.h" 
#include "AreaTriangles.h"
#include "BendingTriangles.h"
#include "LinearSprings.h"
#include "LJSprings.h"
#include "NodeAdvance.h"
#include "Storage.h" 
#include "Edgeswap_test.h"

System::System() {};

void System::solveSystem(){

	Edgeswap edgeswap(coordInfoVecs);

	bool runSim = true; 
	while (runSim == true) { 
		
		for (unsigned i = 0; i < generalParams.solve_time; i++) {
			generalParams.iteration = i;
			thrust::fill(coordInfoVecs.nodeForceX.begin(), coordInfoVecs.nodeForceX.end(), 0.0);
			thrust::fill(coordInfoVecs.nodeForceY.begin(), coordInfoVecs.nodeForceY.end(), 0.0);
			thrust::fill(coordInfoVecs.nodeForceZ.begin(), coordInfoVecs.nodeForceZ.end(), 0.0);
	
			ComputeLinearSprings( 
				generalParams, 
				coordInfoVecs,
				linearSpringInfoVecs, 
				ljInfoVecs);
 
			ComputeAreaTriangleSprings(
				generalParams,
				coordInfoVecs,
				areaTriangleInfoVecs);
			
			ComputeCosTriangleSprings(
				generalParams,
				coordInfoVecs,  
				bendingTriangleInfoVecs); 
			
			ComputeLJSprings(
				coordInfoVecs,
				ljInfoVecs,
				generalParams);

			//now forces are computed, move nodes.
			AdvancePositions(
				coordInfoVecs,
				generalParams,
				domainParams);
				

			AdvanceLJParticle(
				generalParams,
				coordInfoVecs,
				ljInfoVecs);
			
		
		} 
		runSim = false;
		//storage->storeVariables();
		std::cout<<"lj points "<< ljInfoVecs.LJ_PosX<< " "<<  ljInfoVecs.LJ_PosY << " "<<  ljInfoVecs.LJ_PosZ << std::endl;
		std::cout<<"lj force "<< ljInfoVecs.forceX<< " "<<  ljInfoVecs.forceY << " "<<  ljInfoVecs.forceZ << std::endl;

		storage->print_VTK_File();
		
		for (unsigned edge = 2; edge < coordInfoVecs.num_edges; edge++  ) {
			edgeswap.edge_swap(
				edge,
				generalParams,
				coordInfoVecs,
				linearSpringInfoVecs,
				bendingTriangleInfoVecs);
		}
		storage->print_VTK_File();
		
		/*for (unsigned i = 0; i < bendingTriangleInfoVecs.tempNodeForceZUnreduced.size(); i++) {
			std::cout<< "1: "<< bendingTriangleInfoVecs.tempNodeForceZUnreduced[i]<<std::endl;
	std::cout<< "2 "<< bendingTriangleInfoVecs.tempNodeForceXUnreduced[i]<<std::endl;
	std::cout<< "2 "<< bendingTriangleInfoVecs.tempNodeForceYUnreduced[i]<<std::endl;
			std::cout<< "2 "<< bendingTriangleInfoVecs.tempNodeForceZUnreduced[i]<<std::endl;
		}*/
				//Test code first. 



		/*std::cout<<" node 11 loc: "<< coordInfoVecs.nodeLocX[11]<< " "<< coordInfoVecs.nodeLocY[11]<< " "<< coordInfoVecs.nodeLocZ[11]<< std::endl;
		std::cout<<" node 11 force: "<< coordInfoVecs.nodeForceX[11]<< " "<< coordInfoVecs.nodeForceY[11]<< " "<< coordInfoVecs.nodeForceZ[11]<< std::endl;
		std::cout<<" node 17 loc: "<< coordInfoVecs.nodeLocX[17]<< " "<< coordInfoVecs.nodeLocY[17]<< " "<< coordInfoVecs.nodeLocZ[17]<< std::endl;
		std::cout<<" node 17 force: "<< coordInfoVecs.nodeForceX[17]<< " "<< coordInfoVecs.nodeForceY[17]<< " "<< coordInfoVecs.nodeForceZ[17]<< std::endl;
		std::cout<<" node 39 loc: "<< coordInfoVecs.nodeLocX[39]<< " "<< coordInfoVecs.nodeLocY[39]<< " "<< coordInfoVecs.nodeLocZ[39]<< std::endl;
		std::cout<<" node 39 force: "<< coordInfoVecs.nodeForceX[39]<< " "<< coordInfoVecs.nodeForceY[39]<< " "<< coordInfoVecs.nodeForceZ[39]<< std::endl;
		std::cout<<" node 40 loc: "<< coordInfoVecs.nodeLocX[40]<< " "<< coordInfoVecs.nodeLocY[40]<< " "<< coordInfoVecs.nodeLocZ[40]<< std::endl;
		std::cout<<" node 40 force: "<< coordInfoVecs.nodeForceX[40]<< " "<< coordInfoVecs.nodeForceY[40]<< " "<< coordInfoVecs.nodeForceZ[40]<< std::endl;
		*/

	}

	
};





void System::assignStorage(std::shared_ptr<Storage> _storage) {
	storage = _storage;
}

//initialize memory for thrust vectors and set coordInfoVecs vals from input. 
void System::initializeSystem(HostSetInfoVecs& hostSetInfoVecs) {
	std::cout<<"Initializing"<<std::endl;

	generalParams.maxNodeCount = hostSetInfoVecs.hostNodeLocX.size();
	coordInfoVecs.num_edges = hostSetInfoVecs.hostEdges2Nodes_1.size();
	coordInfoVecs.num_triangles = hostSetInfoVecs.hostTriangles2Nodes_1.size();

	std::cout<<"num nodes: "<< generalParams.maxNodeCount << std::endl;
	std::cout<<"num edges: "<< coordInfoVecs.num_edges << std::endl;
	std::cout<<"num elems: "<< coordInfoVecs.num_triangles << std::endl;
	//allocate memory
	coordInfoVecs.isNodeFixed.resize(hostSetInfoVecs.hostNodeLocX.size());
	coordInfoVecs.prevNodeLocX.resize(hostSetInfoVecs.hostNodeLocX.size());
	coordInfoVecs.prevNodeLocY.resize(hostSetInfoVecs.hostNodeLocX.size());
	coordInfoVecs.prevNodeLocZ.resize(hostSetInfoVecs.hostNodeLocX.size());

	coordInfoVecs.prevNodeForceX.resize(hostSetInfoVecs.hostNodeLocX.size());
	coordInfoVecs.prevNodeForceY.resize(hostSetInfoVecs.hostNodeLocX.size());
	coordInfoVecs.prevNodeForceZ.resize(hostSetInfoVecs.hostNodeLocX.size());
	
	coordInfoVecs.nodeLocX.resize(hostSetInfoVecs.hostNodeLocX.size());
	coordInfoVecs.nodeLocY.resize(hostSetInfoVecs.hostNodeLocX.size());
	coordInfoVecs.nodeLocZ.resize(hostSetInfoVecs.hostNodeLocX.size());

	coordInfoVecs.nodeForceX.resize(hostSetInfoVecs.hostNodeLocX.size());
	coordInfoVecs.nodeForceY.resize(hostSetInfoVecs.hostNodeLocX.size());
	coordInfoVecs.nodeForceZ.resize(hostSetInfoVecs.hostNodeLocX.size());

	coordInfoVecs.triangles2Nodes_1.resize( coordInfoVecs.num_triangles );
	coordInfoVecs.triangles2Nodes_2.resize( coordInfoVecs.num_triangles );
	coordInfoVecs.triangles2Nodes_3.resize( coordInfoVecs.num_triangles );
	
	coordInfoVecs.triangles2Edges_1.resize( coordInfoVecs.num_triangles );
	coordInfoVecs.triangles2Edges_2.resize( coordInfoVecs.num_triangles );
	coordInfoVecs.triangles2Edges_3.resize( coordInfoVecs.num_triangles );

	coordInfoVecs.edges2Nodes_1.resize( coordInfoVecs.num_edges );
	coordInfoVecs.edges2Nodes_2.resize( coordInfoVecs.num_edges );
	
	coordInfoVecs.edges2Triangles_1.resize( coordInfoVecs.num_edges );
	coordInfoVecs.edges2Triangles_2.resize( coordInfoVecs.num_edges );



	//copy info to GPU
	std::cout<<"Copying"<<std::endl;
	thrust::copy(hostSetInfoVecs.hostIsNodeFixed.begin(),hostSetInfoVecs.hostIsNodeFixed.end(), coordInfoVecs.isNodeFixed.begin());
	thrust::fill(coordInfoVecs.nodeForceX.begin(), coordInfoVecs.nodeForceX.end(), 0.0);
	thrust::fill(coordInfoVecs.nodeForceY.begin(), coordInfoVecs.nodeForceY.end(), 0.0);
	thrust::fill(coordInfoVecs.nodeForceZ.begin(), coordInfoVecs.nodeForceZ.end(), 0.0);

	thrust::fill(coordInfoVecs.prevNodeForceX.begin(), coordInfoVecs.prevNodeForceX.end(), 0.0);
	thrust::fill(coordInfoVecs.prevNodeForceY.begin(), coordInfoVecs.prevNodeForceY.end(), 0.0);
	thrust::fill(coordInfoVecs.prevNodeForceZ.begin(), coordInfoVecs.prevNodeForceZ.end(), 0.0);
	
	thrust::copy(hostSetInfoVecs.hostNodeLocX.begin(), hostSetInfoVecs.hostNodeLocX.end(), coordInfoVecs.prevNodeLocX.begin() );
	thrust::copy(hostSetInfoVecs.hostNodeLocY.begin(), hostSetInfoVecs.hostNodeLocY.end(), coordInfoVecs.prevNodeLocY.begin() );
	thrust::copy(hostSetInfoVecs.hostNodeLocZ.begin(), hostSetInfoVecs.hostNodeLocZ.end(), coordInfoVecs.prevNodeLocZ.begin() );
	
	thrust::copy(hostSetInfoVecs.hostNodeLocX.begin(), hostSetInfoVecs.hostNodeLocX.end(), coordInfoVecs.nodeLocX.begin() );
	thrust::copy(hostSetInfoVecs.hostNodeLocY.begin(), hostSetInfoVecs.hostNodeLocY.end(), coordInfoVecs.nodeLocY.begin() );
	thrust::copy(hostSetInfoVecs.hostNodeLocZ.begin(), hostSetInfoVecs.hostNodeLocZ.end(), coordInfoVecs.nodeLocZ.begin() );
	
	thrust::copy(hostSetInfoVecs.hostTriangles2Nodes_1.begin(), hostSetInfoVecs.hostTriangles2Nodes_1.end(), coordInfoVecs.triangles2Nodes_1.begin() );
	thrust::copy(hostSetInfoVecs.hostTriangles2Nodes_2.begin(), hostSetInfoVecs.hostTriangles2Nodes_2.end(), coordInfoVecs.triangles2Nodes_2.begin() );
	thrust::copy(hostSetInfoVecs.hostTriangles2Nodes_3.begin(), hostSetInfoVecs.hostTriangles2Nodes_3.end(), coordInfoVecs.triangles2Nodes_3.begin() );
	
	thrust::copy(hostSetInfoVecs.hostTriangles2Edges_1.begin(), hostSetInfoVecs.hostTriangles2Edges_1.end(), coordInfoVecs.triangles2Edges_1.begin() );
	thrust::copy(hostSetInfoVecs.hostTriangles2Edges_2.begin(), hostSetInfoVecs.hostTriangles2Edges_2.end(), coordInfoVecs.triangles2Edges_2.begin() );
	thrust::copy(hostSetInfoVecs.hostTriangles2Edges_3.begin(), hostSetInfoVecs.hostTriangles2Edges_3.end(), coordInfoVecs.triangles2Edges_3.begin() );

	thrust::copy(hostSetInfoVecs.hostEdges2Nodes_1.begin(), hostSetInfoVecs.hostEdges2Nodes_1.end(), coordInfoVecs.edges2Nodes_1.begin() );
	thrust::copy(hostSetInfoVecs.hostEdges2Nodes_2.begin(), hostSetInfoVecs.hostEdges2Nodes_2.end(), coordInfoVecs.edges2Nodes_2.begin() );
	
	thrust::copy(hostSetInfoVecs.hostEdges2Triangles_1.begin(), hostSetInfoVecs.hostEdges2Triangles_1.end(), coordInfoVecs.edges2Triangles_1.begin() );
	thrust::copy(hostSetInfoVecs.hostEdges2Triangles_2.begin(), hostSetInfoVecs.hostEdges2Triangles_2.end(), coordInfoVecs.edges2Triangles_2.begin() );

 
	//allocate memory for other data structures.   

	//area triangle info vec
	//number of area springs is the number of triangles
	std::cout<<"Mem"<<std::endl;
	areaTriangleInfoVecs.tempNodeIdUnreduced.resize(areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceXUnreduced.resize(areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceYUnreduced.resize(areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceZUnreduced.resize(areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	
	areaTriangleInfoVecs.tempNodeIdReduced.resize(areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceXReduced.resize(areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceYReduced.resize(areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceZReduced.resize(areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);

	//beinding triangle info vec
	//num bending springs is the number of times each edge is between two triangles. 
	bendingTriangleInfoVecs.numBendingSprings = coordInfoVecs.edges2Triangles_1.size();

	bendingTriangleInfoVecs.tempNodeIdUnreduced.resize(bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceXUnreduced.resize(bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceYUnreduced.resize(bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceZUnreduced.resize(bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	
	bendingTriangleInfoVecs.tempNodeIdReduced.resize(bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceXReduced.resize(bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceYReduced.resize(bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceZReduced.resize(bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);

	//linear springs
	linearSpringInfoVecs.tempNodeIdUnreduced.resize(linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceXUnreduced.resize(linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceYUnreduced.resize(linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceZUnreduced.resize(linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
	
	linearSpringInfoVecs.tempNodeIdReduced.resize(linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceXReduced.resize(linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceYReduced.resize(linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceZReduced.resize(linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
	
	linearSpringInfoVecs.edge_initial_length.resize(coordInfoVecs.num_edges);
	
	thrust::copy(hostSetInfoVecs.hostEdge_initial_length.begin(), hostSetInfoVecs.hostEdge_initial_length.end(), linearSpringInfoVecs.edge_initial_length.begin() );
	std::cout<<"initial lengths: "<< linearSpringInfoVecs.edge_initial_length.size()<<std::endl;

	std::cout<<"System Ready"<<std::endl;

	//Generate LJ particle list. and set LJ particle midpoint.
	double maxX_lj = *(thrust::max_element(coordInfoVecs.nodeLocX.begin(),coordInfoVecs.nodeLocX.end()));
	double minX_lj = *(thrust::min_element(coordInfoVecs.nodeLocX.begin(),coordInfoVecs.nodeLocX.end()));
	double maxY_lj = *(thrust::max_element(coordInfoVecs.nodeLocY.begin(),coordInfoVecs.nodeLocY.end()));
	double minY_lj = *(thrust::min_element(coordInfoVecs.nodeLocY.begin(),coordInfoVecs.nodeLocY.end()));
	
	ljInfoVecs.LJ_PosX = (maxX_lj + minX_lj)/2.0;
	ljInfoVecs.LJ_PosY = (maxY_lj + minY_lj)/2.0;


	//currently unused
	thrust::host_vector<unsigned> tempIds;
	for (unsigned i = 0; i < hostSetInfoVecs.hostNodeLocX.size(); i++ ) {
		double xLoc = hostSetInfoVecs.hostNodeLocX[i];
		double yLoc = hostSetInfoVecs.hostNodeLocY[i];
		double zLoc = hostSetInfoVecs.hostNodeLocZ[i];
		
		double xDist = ljInfoVecs.LJ_PosX - xLoc;
		double yDist = ljInfoVecs.LJ_PosY - yLoc;
		double zDist = ljInfoVecs.LJ_PosZ - zLoc;

		double dist = std::sqrt(xDist*xDist + yDist*yDist + zDist*zDist);
		//just test all poitns for now. Optimize later.
		if (dist < 100 * ljInfoVecs.Rcutoff) {
			tempIds.push_back(i);
		}
	}
	ljInfoVecs.node_id_close.resize( tempIds.size() );
	thrust::copy(tempIds.begin(), tempIds.end(), ljInfoVecs.node_id_close.begin());
	std::cout<<"lj nodes: "<< ljInfoVecs.node_id_close.size() << std::endl;


};


