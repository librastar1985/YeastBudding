################################################################################
# Automatically-generated file. Do not edit!
################################################################################
CXX := gcc
NVCC := nvcc
CFLAGS := -pthread -std=c++11 -Wall -Wextra
MFLAGS := -f
NVCCFLAGS := -std=c++11

current_dir := $(shell pwd)
LIBS:=  -lpugixml -L/$(current_dir)/pugixml/lib64
#-lgsl -lgslcblas

ILIBS_cuda8 = -I/opt/linux/centos/7.x/x86_64/pkgs/cuda/8.0/include/
ILIBS_cuda9 := -I/opt/linux/centos/7.x/x86_64/pkgs/cuda/9.1/include/

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../TurgorForce.cu \
../AreaTriangles.cu \
../BendingTriangles.cu \
../MemRepulsionSprings_universal.cu\
../MemRepulsionSprings_local.cu\
../LinearSprings.cu \
../AreaCompBud.cu \
../VolumeComp.cu \
../VolumeSprings.cu \
../LineTensionSprings.cu \
../NodeAdvance.cu \
../MemRepulsionEnergy.cu \
../System.cu \
../Edgeswap_test.cpp \
../SystemBuilder.cpp \
../Storage.cpp \
../main.cpp


# this is a variable
OBJS += \
./TurgorForce.o \
./AreaTriangles.o \
./BendingTriangles.o \
./MemRepulsionSprings_universal.o\
./MemRepulsionSprings_local.o\
./LinearSprings.o \
./AreaCompBud.o \
./VolumeComp.o \
./VolumeSprings.o \
./LineTensionSprings.o \
./NodeAdvance.o \
./MemRepulsionEnergy.o \
./System.o \
./Edgeswap_test.o \
./SystemBuilder.o \
./Storage.o \
./main.o


CPP_DEPS += \
./TurgorForce.d \
./AreaTriangles.d \
./BendingTriangles.d \
./MemRepulsionSprings_universal.d\
./MemRepulsionSprings_local.d\
./LinearSprings.d \
./AreaCompBud.d \
./VolumeComp.d \
./VolumeSprings.d \
./LineTensionSprings.d \
./NodeAdvance.d \
./MemRepulsionEnergy.d \
./System.d \
./Edgeswap_test.d \
./SystemBuilder.d \
./Storage.d \
./main.d

#need o have ILIBS2
#cpp files
%.o : ./%.cpp 
	 $(CXX) $(CFLAGS) $(ILIBS_cuda8) $(LIBS) -o $@ -c $^ 

	
#cuda files
%.o : ./%.cu 
	$(NVCC) $(NVCCFLAGS) $(ILIBS_cuda9) -dc -o $@ $^