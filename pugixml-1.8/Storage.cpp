

#include "System.h"
#include "Storage.h"

Storage::Storage(
	std::weak_ptr<System> system_) {
	//std::cout << "FDM constructor" << std::endl;
	
	system = system_;

	std::ofstream statesOutput("Temp.sta");
	

	std::shared_ptr<System> SYSTEM = system.lock();//upgrades weak to shared

	
	if (SYSTEM) {
		
	 	statesOutput << "node_count " << SYSTEM->coordInfoVecs.nodeLocX.size() << '\n';
	 	statesOutput << "edge_count " << SYSTEM->coordInfoVecs.num_edges << '\n';
	 	statesOutput << "elem_count " << SYSTEM->coordInfoVecs.num_triangles << '\n';
	}


	statesOutput.close();
}


void Storage::print_VTK_File(void) {
	
	std::shared_ptr<System> SYSTEM = system.lock();

	if ((SYSTEM)) {

		iteration+=1;
		unsigned digits = ceil(log10(iteration + 1));
		std::string format = ".vtk";
		std::string Number;
		std::string initial = "AnimationTest/Vir_";
		std::ofstream ofs;
		if (digits == 1 || digits == 0) {
			Number = "0000" + std::to_string(iteration);
		}
		else if (digits == 2) {
			Number = "000" + std::to_string(iteration);
		}
		else if (digits == 3) {
			Number = "00" + std::to_string(iteration);
		}
		else if (digits == 4) {
			Number = "0" + std::to_string(iteration);
		}

		std::string Filename = initial + Number + format;

		ofs.open(Filename.c_str());
		
	
		unsigned numParticles = SYSTEM->generalParams.maxNodeCount;//one for lj particle

		
		

		
		ofs << "# vtk DataFile Version 3.0" << std::endl;
		ofs << "Point representing Sub_cellular elem model" << std::endl;
		ofs << "ASCII" << std::endl << std::endl;
		ofs << "DATASET UNSTRUCTURED_GRID" << std::endl;
		
		 
		ofs << "POINTS " <<numParticles + 1 << " float" << std::endl;
		for (unsigned i = 0; i< numParticles; i++) {
			double xPos = SYSTEM->coordInfoVecs.nodeLocX[i];
			double yPos = SYSTEM->coordInfoVecs.nodeLocY[i];
			double zPos = SYSTEM->coordInfoVecs.nodeLocZ[i];
			
			ofs << std::setprecision(5) <<std::fixed<< xPos << " " << yPos << " " << zPos << " " << '\n'<< std::fixed;
		}
		ofs << std::setprecision(5) <<std::fixed<< SYSTEM->ljInfoVecs.LJ_PosX << " " << SYSTEM->ljInfoVecs.LJ_PosY << " " << SYSTEM->ljInfoVecs.LJ_PosZ << " " << '\n'<< std::fixed;
		
		

		
		unsigned numEdges = SYSTEM->coordInfoVecs.num_edges;
		unsigned numCells = numEdges + 1;//one cell for LJ Particle, rest for edges of polymer
		unsigned numNumsInCells = (3 * numEdges) + (2);//add one for lj and one to list it.
		
		
		ofs << "CELLS " << numCells << " " << numNumsInCells << std::endl;
		//place edges as cells of type 2. 
		for (unsigned edge = 0; edge < numEdges; edge++ ){
			unsigned idA = SYSTEM->coordInfoVecs.edges2Nodes_1[edge];
			unsigned idB = SYSTEM->coordInfoVecs.edges2Nodes_2[edge];

			ofs<< 2 << " " << idA << " " << idB << std::endl;
			
		}

		ofs<< 1 << " " << SYSTEM->generalParams.maxNodeCount << std::endl;
		
		ofs << "CELL_TYPES " << numCells << std::endl;  
		//set edges and last set scattered points(dpd)
		for (unsigned i = 0; i < (numEdges); i++) {
				
			ofs << 3 << std::endl;//edge (2d line)
		}
		ofs << 1 << std::endl;//edge (2d line)
		
	
		ofs << "CELL_DATA " << numCells << std::endl;
		ofs << "SCALARS Strain double " << std::endl;
		ofs << "LOOKUP_TABLE default "  << std::endl;
		//set strain for each edge
		for (unsigned edge = 0; edge < numEdges; edge++ ){

			unsigned idA = SYSTEM->coordInfoVecs.edges2Nodes_1[edge];
			unsigned idB = SYSTEM->coordInfoVecs.edges2Nodes_2[edge];
			if (idA >= SYSTEM->generalParams.maxNodeCount)
				std::cout<<idA<<std::endl;
			if (idB >= SYSTEM->generalParams.maxNodeCount)
				std::cout<<idB<<std::endl;
			double L0 = SYSTEM->linearSpringInfoVecs.edge_initial_length[edge];
			double xL = SYSTEM->coordInfoVecs.nodeLocX[idA];
			double yL = SYSTEM->coordInfoVecs.nodeLocY[idA];
			double zL = SYSTEM->coordInfoVecs.nodeLocZ[idA];
			double xR = SYSTEM->coordInfoVecs.nodeLocX[idB];
			double yR = SYSTEM->coordInfoVecs.nodeLocY[idB];
			double zR = SYSTEM->coordInfoVecs.nodeLocZ[idB];
			
			double L1 = std::sqrt( (xL - xR)*(xL - xR)+(yL - yR)*(yL - yR)+(zL - zR)*(zL - zR));
			double strain = (L1 - L0) / L0;
			ofs << std::fixed << strain   << std::endl;
				
			
		}
		ofs << std::fixed << 0.1   << std::endl;
			
		ofs.close();
	
	}
};

void Storage::storeVariables(void) {
	std::shared_ptr<System> SYSTEM = system.lock();
	if (SYSTEM) {

		//first create a new file using the current network strain
		
		std::string format = ".sta";
		std::string lj_z =  std::to_string(SYSTEM->ljInfoVecs.LJ_PosZ);
		std::string initial = "Variables/File_";
		std::ofstream ofs;
		std::string Filename = initial + lj_z + format;
		ofs.open(Filename.c_str());




		double total_energy =  SYSTEM->linearSpringInfoVecs.linear_spring_energy + 
        						SYSTEM->areaTriangleInfoVecs.area_triangle_energy + 
        						SYSTEM->bendingTriangleInfoVecs.bending_triangle_energy + 
        						SYSTEM->ljInfoVecs.lj_energy;
								
		ofs << std::setprecision(5) <<std::fixed<< "total_energy=" << total_energy<<std::endl;
		
		ofs << std::setprecision(5) <<std::fixed<< "lj_x_pos=" << SYSTEM->ljInfoVecs.LJ_PosX<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "lj_y_pos=" << SYSTEM->ljInfoVecs.LJ_PosY<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "lj_z_pos=" << SYSTEM->ljInfoVecs.LJ_PosZ<<std::endl;
		

		//place nodes
		for (unsigned i = 0; i < SYSTEM->coordInfoVecs.nodeLocX.size(); i++) {
			double x = SYSTEM->coordInfoVecs.nodeLocX[i];
			double y = SYSTEM->coordInfoVecs.nodeLocY[i];
			double z = SYSTEM->coordInfoVecs.nodeLocZ[i];
			ofs << std::setprecision(5) <<std::fixed<< "node=" << x << " " << y << " " << z <<std::endl;
		
		}

		for (unsigned i = 0; i < SYSTEM->coordInfoVecs.nodeLocX.size(); i++) {
			unsigned t2n_1 = SYSTEM->coordInfoVecs.triangles2Nodes_1[i];
			unsigned t2n_2 = SYSTEM->coordInfoVecs.triangles2Nodes_2[i];
			unsigned t2n_3 = SYSTEM->coordInfoVecs.triangles2Nodes_3[i];
			ofs << std::setprecision(5) <<std::fixed<< "t2n=" << t2n_1 << " " << t2n_2 << " " << t2n_3 <<std::endl;
		
		}
	

	}
}
