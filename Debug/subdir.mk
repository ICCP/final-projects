################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../CChromosome.cpp \
../CCoordSet.cpp \
../CPlanet.cpp \
../CSVector.cpp \
../CSatellite.cpp \
../CSimulation.cpp \
../CSolarSystem.cpp \
../Parameters.cpp \
../main.cpp 

OBJS += \
./CChromosome.o \
./CCoordSet.o \
./CPlanet.o \
./CSVector.o \
./CSatellite.o \
./CSimulation.o \
./CSolarSystem.o \
./Parameters.o \
./main.o 

CPP_DEPS += \
./CChromosome.d \
./CCoordSet.d \
./CPlanet.d \
./CSVector.d \
./CSatellite.d \
./CSimulation.d \
./CSolarSystem.d \
./Parameters.d \
./main.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -std=c++11 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


