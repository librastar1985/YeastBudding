

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include subdir.mk


ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
endif

-include ../makefile.defs


# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: dpd-model

#Flags 
CXXFLAGS=-O3 -std=c++0x -pg -g -c -Wall
NVCCFLAGS=-std=c++11 

LIBS := -L/afs/crc.nd.edu/x86_64_linux/m/matlab/R2018a/bin/glnxa64/ -L/afs/crc.nd.edu/x86_64_linux/m/matlab/R2018a/extern/bin/glnxa64/
ILIBS := -I/afs/crc.nd.edu/x86_64_linux/c/cuda/8.0/include/ -I/afs/crc.nd.edu/x86_64_linux/m/matlab/R2018a/extern/include/
LIBS1 := -lMatlabEngine -lMatlabDataArray -lmat -leng -lmx -lmex

# Tool invocations 
dpd-model: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: NVCC C++ Linker'
	nvcc -shared -o "dpd-model" $(OBJS) $(USER_OBJS) $(ILIBS) $(LIBS) $(NVCCFLAGS) 
	@echo 'Finished building target: $@'
	@echo ' '
	$(MAKE) --no-print-directory post-build 

# Other Targets
clean:
	-$(RM) $(OBJS)$(C++_DEPS)$(C_DEPS)$(CC_DEPS)$(CPP_DEPS)$(EXECUTABLES)$(CXX_DEPS)$(C_UPPER_DEPS) dpd-model
	-@echo ' '

post-build:
	-mkdir --parents ../../dpd-0.0.6/Debug; cd ../../dpd-0.0.6/Debug; cp ../../dpd-model/Debug/dpd-model dpd-model
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY: post-build

-include ../makefile.targets
# DO NOT DELETE
