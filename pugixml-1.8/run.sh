#!/bin/csh

#to run without script you must type this line into the terminal:
setenv LD_LIBRARY_PATH "${LD_LIBRARY_PATH}:/afs/crc.nd.edu/x86_64_linux/m/matlab/R2017b/extern/bin/glnxa64/"

#$ -M 	brit004@ucr.edu # Email address for job notification
#$ -m  abe		 # Send mail when job begins, ends and aborts
#$ -q  gpu 	 # Specify queue
#$ -l gpu_card=1
#s -pe smp 4         #specifies threads??? maybe
#$ -N  DPD	 # Specify job name
#$ -t 1       #specify number of data input files


set data = ( Data_structure.xml )

module purge
module load gcc/6.2.0
module load cuda
echo -n "It is currently: ";date
echo -n "I am logged on as ";whoami
echo -n "This computer is called ";hostname
echo -n "I am currently in the directory ";pwd


./virus-model  -solve_time=200 -dt=0.00001 $data[${SGE_TASK_ID}]

