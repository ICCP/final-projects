################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../CCoordSet.cpp \
../CPlanet.cpp \
../CSVector.cpp \
../CSimulation.cpp \
../CSolarSystem.cpp \
../Parameters.cpp \
../main.cpp 

OBJS += \
./CCoordSet.o \
./CPlanet.o \
./CSVector.o \
./CSimulation.o \
./CSolarSystem.o \
./Parameters.o \
./main.o 

CPP_DEPS += \
./CCoordSet.d \
./CPlanet.d \
./CSVector.d \
./CSimulation.d \
./CSolarSystem.d \
./Parameters.d \
./main.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


