/*
 * Parameters.cpp
 *
 *  Created on: Apr 8, 2015
 *      Author: Mark Wittekind and Drew Murray
 */

#ifndef CPARAMETERS_H_
#define CPARAMETERS_H_

#include <math.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <string>

using namespace std;

struct Par{
	//Astronomical units (Me = 1, Re = 1);
	float G = 6.6738480 * pow(10,-11); //m3*kg-1*s-2
	float SolarMass = 1.98855*pow(10.0,30.0);
	float AU = 149597870700; //meters
	double boundaryDistanceSquared = pow(7311000000000*3,2);//must be greater than pluto's

	//Single simulation integration parameters
	float dt = 1*pow(10,4);
	float maxT = 5*23740*dt;
	unsigned long int maxTimeSteps = (maxT/dt);
	unsigned short integration_method = 4;
	unsigned short print_freq = 1;
	//0 = errorflag
	//1 = Euler
	//2 = 2nd order Runga-Kutta
	//4 = 4th order Runga-Kutta

	//Realtime for 2.4*1000000*10^4 of Mars with RK4: 1345017204 nanoseconds
	//Realtime for 20 orbits of mars with Euler: 8:44.6

	unsigned short runMode = 4;
	//0 = errorflag
	//1 = GA
	//2 = plot positions
	//3 = Energy (RK vs Euler)
	//4 = Energy (Kinetic vs Potential)

	//GA Params
	float mutateChance = .2;
	unsigned int satsPerCore = 0;
	unsigned int randSatsPerGen = 0;
	unsigned int maxGenerations = 1;
	float breedingNum = 0;//# will breed
	float randomBreedingNum = 0;//# of breeding pop will be selected without bias for fitness

	// Satellite params
	float maxThrusterVel = 15000;
	unsigned int thrustQuant = 4;
	float satMass = 722;
	string homePlanetName = "Earth";

	//fitness params
	string destinationPlanetName = "Mars";
	float goalDistanceMultiplier = 1;//6650;//multiplied by destiantion planet radius (must be >= 1)
	//multiplied by pluto's aphelion


	//TODO List:
	/*
	 * Mode 2 should pull satellites from success list (end)
	 * Fix collision printing so it doesn't print multiple instances before sat name & timesteps
	 * Parallelize loop of simulations
	 * Write loop for parallelization
	 * Parameter for max run time (realtime) & how often to check
	 * Run on HPCC
	 * Improve breeding algorithm
	 * Increase demands on GA with a more specific fitness function
	 */

};

struct SatName{
	unsigned short coreNum = 0;
	unsigned int genNum = 0;
	unsigned long long selfID = 0;
	unsigned long long momID = 0;
	unsigned long long dadID = 0;
};

#endif /* CPARAMETERS_H_ */
