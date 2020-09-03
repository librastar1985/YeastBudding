# YeastBudding
This code simulates the yeast budding process.
The instruction provided here is for the code to run on UCR HPCC cluster. Modification is likely needed if running on a different server.

To modify parameter sets, please refer to System.cu and Datastructure.xml.
Datastructure.xml contains the parameter set that initialize the elasticity of the model cell and data structure detailing information such as the neighboring triangles.
System.cu contains modifiable parameter set such as the elasticity of the bud site, chitin and septin ring stiffness.

To ensure that the simulation output can be recorded correctly, please make modification in Storage.cpp and Sbatch.sh.

Before running the code, make sure the current directory points to the the root folder containing these scripts.
To run the code on HPCC, please use the following commands:
1. module load extra; module load GCC; module load cuda/9.1
2. make
3. sbatch -p gpu --gres=gpu:1 --mem=2g --time=48:00:00 Sbatch.sh
