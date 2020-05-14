################################################################################
# Automatically-generated file. Do not edit!
################################################################################
CXX := gcc
NVCC := nvcc
CFLAGS := -pthread -std=c++11 -Wall -Wextra
MFLAGS := -f
NVCCFLAGS := -std=c++11
LIBS11 := -L/afs/crc.nd.edu/x86_64_linux/m/matlab/R2017b/bin/glnxa64/
LIBS1 := -L/afs/crc.nd.edu/x86_64_linux/m/matlab/R2017b/extern/bin/glnxa64/
LIBSFINAL := -lMatlabEngine -lMatlabDataArray 

ILIBS1 := -I/afs/crc.nd.edu/x86_64_linux/c/cuda/8.0/include/ 
ILIBS2 := -I/afs/crc.nd.edu/x86_64_linux/m/matlab/R2017b/extern/include/

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../engdemo.cpp


# this is a variable
OBJS += \
./engdemo.o 


CPP_DEPS += \
./engdemo.d 

#need o have ILIBS2
#cpp files
%.o : ./%.cpp 
	 $(CXX) $(CFLAGS) $(ILIBS2) $(LIBS11) $(LIBS1) -o $@ -c $^ $(LIBSFINAL)

	
#cuda files
%.o : ./%.cu 
	$(NVCC) $(NVCCFLAGS) -dc -o $@ $^ 



#g++ -std=c++11 -I/afs/crc.nd.edu/x86_64_linux/m/matlab/R2018a/extern/include/ -L/afs/crc.nd.edu/x86_64_linux/m/matlab/R2018a/bin/glnxa64/ -pthread engdemo.cpp -lMatlabDataArray -lMatlabEngine
