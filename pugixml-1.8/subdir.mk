################################################################################
# Automatically-generated file. Do not edit!
################################################################################
CXX := gcc
NVCC := nvcc
CFLAGS := -pthread -std=c++11 -Wall -Wextra
MFLAGS := -f
NVCCFLAGS := -std=c++11

LIBS := -lpugixml -L/afs/crc.nd.edu/user/s/sbritto2/Kevin_Code_CAll/pugixml/lib64
LIBS11 := -L/afs/crc.nd.edu/x86_64_linux/m/matlab/R2017b/bin/glnxa64/
LIBS1 := -L/afs/crc.nd.edu/x86_64_linux/m/matlab/R2017b/extern/bin/glnxa64/
LIBSFINAL := -lMatlabEngine -lMatlabDataArray 

ILIBSCPP := -I/afs/crc.nd.edu/x86_64_linux/c/cuda/8.0/include/
ILIBS1 := -I/afs/crc.nd.edu/x86_64_linux/c/cuda/9.1/include/
ILIBS2 := -I/afs/crc.nd.edu/x86_64_linux/m/matlab/R2017b/extern/include/

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../AreaTriangles.cu \
../BendingTriangles.cu \
../LinearSprings.cu \
../LJSprings.cu \
../NodeAdvance.cu \
../System.cu \
../Edgeswap_test.cpp \
../SystemBuilder.cpp \
../Storage.cpp \
../main.cpp


# this is a variable
OBJS += \
./AreaTriangles.o \
./BendingTriangles.o \
./LinearSprings.o \
./LJSprings.o \
./NodeAdvance.o \
./System.o \
./Edgeswap_test.o \
./SystemBuilder.o \
./Storage.o \
./main.o


CPP_DEPS += \
./AreaTriangles.d \
./BendingTriangles.d \
./LinearSprings.d \
./LJSprings.d \
./NodeAdvance.d \
./System.d \
./Edgeswap_test.d \
./SystemBuilder.d \
./Storage.d \
./main.d

#need o have ILIBS2
#cpp files
%.o : ./%.cpp 
	 $(CXX) $(CFLAGS) $(ILIBSCPP) $(ILIBS2) $(LIBS) $(LIBS1) $(LIBS11) -o $@ -c $^ $(LIBSFINAL)

	
#cuda files
%.o : ./%.cu 
	$(NVCC) $(NVCCFLAGS) $(ILIBS1) $(ILIBS2) $(LIBS11) $(LIBS1) -dc -o $@ $^ $(LIBSFINAL)
