################################################################################
# Automatically-generated file. Do not edit!
################################################################################
CXX := gcc
MXX := mex
NVCC := nvcc
CFLAGS := -fPIC -std=c++11 -Wall -Wextra
MFLAGS := -f
NVCCFLAGS := -std=c++11
IFLAGS := -I/afs/crc.nd.edu/x86_64_linux/c/cuda/8.0/include/ -I/afs/crc.nd.edu/x86_64_linux/m/matlab/R2018a/extern/include/ 
LIBS := -L/afs/crc.nd.edu/x86_64_linux/m/matlab/R2018a/bin/glnxa64/ -L/afs/crc.nd.edu/x86_64_linux/m/matlab/R2018a/extern/bin/glnxa64/


# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../engdemo.cpp


# this is a variable
OBJS += \
./engdemo.o 

 
CPP_DEPS += \
./engdemo.d 


#cpp files
%.o : ./%.cpp 
	$(CXX) $(CFLAGS) $(IFLAGS) -o  $@ -c $^

	
#cuda files
%.o : ./%.cu 
	$(NVCC) $(NVCCFLAGS) -dc -o $@ $^ 



