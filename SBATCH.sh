#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --output=YB_cellwall4_l1d0b0d05_st0d05_dampedring400_longsim_nontipdoublecoefRandge1d0_SCALETYPE1_maxbendscale0d2_theta0d0# This affects the print out of the "std::cout" in the script, make sure this is changed for different jobs.
#SBATCH --mail-user=useremail@address.com
#SBATCH --mail-type=ALL
#SBATCH --job-name="tip"
#SBATCH -p gpu # This is the default partition, you can use any of the following; intel, batch, highmem, gpu

module load extra
module load GCC
module load cuda/9.1


srun -p gpu --gres=gpu:1 ./virus-model -dt=0.001 Data_structure.xml

