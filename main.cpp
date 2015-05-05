/*
 * main.cpp
 *
 *  Created on: Apr 8, 2015
 *      Author: Mark Wittekind and Drew Murray
 */
#include "CSolarSystem.h"
#include "Parameters.cpp"
#include <algorithm>
#include <fstream>
#include "CSimulation.h"
#include <sys/time.h>
#include "CChromosome.h"
#include <random>
#include <chrono>
#include "CSatellite.h"
#include <string>

using namespace std;

int main() {
	streambuf * pos_buf;
	ofstream pos_of;
	//set to filename
	pos_of.open("data/raw_positions.txt");
	pos_buf = pos_of.rdbuf();
	ostream pos_stream(pos_buf);

	streambuf * en_buf;
	ofstream en_of;
	en_of.open("data/delta_E.txt");
	en_buf = en_of.rdbuf();
	ostream en_stream(en_buf);

	streambuf * param_buf;
	ofstream param_of;
	param_of.open("data/parameters.txt");
	param_buf = param_of.rdbuf();
	ostream param_stream(param_buf);

	streambuf * com_buf;
	ofstream com_of;
	com_of.open("data/center_of_mass.txt");
	com_buf = com_of.rdbuf();
	ostream com_stream(com_buf);

	streambuf * success_buf;
	ofstream success_of;
	success_of.open("data/success.txt");
	success_buf = success_of.rdbuf();
	ostream success_stream(success_buf);

	Par par;
	CSolarSystem newSS(&par);
	// 2nd int params are gen,corenum
	CSimulation RKSim(newSS, &par, "RK4", 0, 0);
	//RKSim.print_pos(pos_stream);

	Par par2;
	par2.integration_method = 1;
	par2.dt /= 4;
	par2.maxT *= 4;
	CSolarSystem newSS2(&par2);
	CSimulation ESim(newSS2, &par2, "Euler", 0, 0);

	double init_en = newSS.totalEnergy();

	//TODO: update parameter output
	param_stream << "-----" << RKSim.m_name << "-----\n";
	param_stream << "AU: " << RKSim.m_par->AU << "\n";
	param_stream << "G: " << RKSim.m_par->G << "\n";
	param_stream << "SolarMass: " << RKSim.m_par->SolarMass << "\n";
	param_stream << "dt: " << RKSim.m_par->dt << "\n";
	param_stream << "maxT: " << RKSim.m_par->maxT << "\n";
	param_stream << "integration_method: " << RKSim.m_par->integration_method
			<< "\n";
	param_stream << "print_freq: " << RKSim.m_par->print_freq << "\n\n";
	param_stream.flush();
	param_stream << "Boundary Distance Squared: "
			<< RKSim.m_par->boundaryDistanceSquared << "\n";
	param_stream << "runMode: " << RKSim.m_par->runMode << "\n";
	param_stream << "mutateChance: " << RKSim.m_par->mutateChance << "\n";
	param_stream << "satsPerCore: " << RKSim.m_par->satsPerCore << "\n";
	param_stream << "maxGenerations: " << RKSim.m_par->maxGenerations << "\n";

	param_stream.flush();

	if (par.runMode == 1) {
		for (unsigned int g = 0; g < par.maxGenerations; g++) {
			for (unsigned long int t = 0; t < RKSim.m_par->maxTimeSteps; t++) {
				if (int(t) % par.print_freq == 0) {
					cout << "Timestep " << t << " of "
							<< RKSim.m_par->maxTimeSteps << "\n";
					cout.flush();
				}
				RKSim.update(t);
			}
			cout << "Generation: " << RKSim.m_genNum << "\n";
			cout.flush();
			RKSim = RKSim * RKSim;
			RKSim.successPrint(success_stream);
			success_stream << "---------\n" << "Generation: " << RKSim.m_genNum
					<< "\n---------\n";
			success_stream.flush();
			RKSim.m_succList.clear();
			RKSim.m_breedList.clear();
			RKSim.m_deadList.clear();
		}
	} else if (par.runMode == 2) {
//		CCoordSet tempSet;
//		RKSim.m_SS.m_sats.clear();
//
//		SatName name1;
//		name1.coreNum = 666;
//		name1.genNum = 14;
//		name1.dadID = 2216910592929146480;
//		name1.momID = 7021797237217784752;
//		name1.selfID = 5453765456158656951;
//		CSatellite sat1(&par, &RKSim.m_SS, name1);
//		sat1.m_thrusts[0].m_dvx[0] = 28672;
//		sat1.m_thrusts[0].m_dvy[0] = -7328.31;
//		sat1.m_thrusts[0].m_dvz[0] = -3773.48;
//		sat1.m_thrusts[0].m_t[0] = 10387.2;
//		sat1.m_thrusts[1].m_dvx[0] = -19364.3;
//		sat1.m_thrusts[1].m_dvy[0] = 26883.8;
//		sat1.m_thrusts[1].m_dvz[0] = 485.551;
//		sat1.m_thrusts[1].m_t[0] = 8728.5;
//		sat1.m_thrusts[2].m_dvx[0] = -19364.3;
//		sat1.m_thrusts[2].m_dvy[0] = 26883.8;
//		sat1.m_thrusts[2].m_dvz[0] = 485.551;
//		sat1.m_thrusts[2].m_t[0] = 8728.5;
//		sat1.m_thrusts[3].m_dvx[0] = -2996.13;
//		sat1.m_thrusts[3].m_dvy[0] = 12409.2;
//		sat1.m_thrusts[3].m_dvz[0] = -17033.3;
//		sat1.m_thrusts[3].m_t[0] = 0;
//		tempSet = RKSim.m_SS.getPlanetDynamics(sat1.m_homePlanetName);
//		tempSet.m_position.m_y = -tempSet.m_position.m_y
//				+ RKSim.m_SS.getPlanetRadius(sat1.m_homePlanetName);
//
//		sat1.setDynamics(tempSet);
//		sat1.m_startDynamics = tempSet;
//
//		sat1.m_thrusts[sat1.m_thrusts.size() - 1].m_t[0] = 0;
//		RKSim.m_SS.m_sats.push_back(sat1);
//
//		SatName name2;
//		name2.coreNum = 666;
//		name2.genNum = 14;
//		name2.dadID = 2216910592929146480;
//		name2.momID = 7021797237217784752;
//		name2.selfID = 5001196753065513141;
//		CSatellite sat2(&par, &RKSim.m_SS, name2);
//		sat2.m_thrusts[0].m_dvx[0] = 28672;
//		sat2.m_thrusts[0].m_dvy[0] = -7328.31;
//		sat2.m_thrusts[0].m_dvz[0] = -3773.48;
//		sat2.m_thrusts[0].m_t[0] = 10387.2;
//		sat2.m_thrusts[1].m_dvx[0] = -19364.3;
//		sat2.m_thrusts[1].m_dvy[0] = 26883.8;
//		sat2.m_thrusts[1].m_dvz[0] = 485.551;
//		sat2.m_thrusts[1].m_t[0] = 8728.5;
//		sat2.m_thrusts[2].m_dvx[0] = -19364.3;
//		sat2.m_thrusts[2].m_dvy[0] = 26883.8;
//		sat2.m_thrusts[2].m_dvz[0] = 485.551;
//		sat2.m_thrusts[2].m_t[0] = 8728.5;
//		sat2.m_thrusts[3].m_dvx[0] = -2996.13;
//		sat2.m_thrusts[3].m_dvy[0] = 12409.2;
//		sat2.m_thrusts[3].m_dvz[0] = -17033.3;
//		sat2.m_thrusts[3].m_t[0] = 0;
//
//		sat2.setDynamics(tempSet);
//		sat2.m_startDynamics = tempSet;
//
//		sat2.m_thrusts[sat2.m_thrusts.size() - 1].m_t[0] = 0;
//		//RKSim.m_SS.m_sats.push_back(sat2);
//
//		SatName name3;
//		name3.coreNum = 666;
//		name3.genNum = 14;
//		name3.dadID = 6676077082895485240;
//		name3.momID = 6716117035617490574;
//		name3.selfID = 6024843162545690288;
//		CSatellite sat3(&par, &RKSim.m_SS, name3);
//		sat3.m_thrusts[0].m_dvx[0] = -27411.8;
//		sat3.m_thrusts[0].m_dvy[1] = 17314.2;
//		sat3.m_thrusts[0].m_dvz[2] = -11792.8;
//		sat3.m_thrusts[0].m_t[3] = 11959;
//		sat3.m_thrusts[1].m_dvx[0] = 16358;
//		sat3.m_thrusts[1].m_dvy[1] = 24760.8;
//		sat3.m_thrusts[1].m_dvz[2] = -6297.76;
//		sat3.m_thrusts[1].m_t[3] = 8950.29;
//		sat3.m_thrusts[2].m_dvx[0] = -13418.9;
//		sat3.m_thrusts[2].m_dvy[1] = -12908.4;
//		sat3.m_thrusts[2].m_dvz[2] = -21812.2;
//		sat3.m_thrusts[2].m_t[3] = 7425.84;
//		sat3.m_thrusts[3].m_dvx[0] = 10539.8;
//		sat3.m_thrusts[3].m_dvy[1] = -17380.6;
//		sat3.m_thrusts[3].m_dvz[2] = 21069.1;
//		sat3.m_thrusts[3].m_t[3] = 0;
//
//		sat3.setDynamics(tempSet);
//		sat3.m_startDynamics = tempSet;
//
//		sat3.m_thrusts[sat3.m_thrusts.size() - 1].m_t[0] = 0;
//		RKSim.m_SS.m_sats.push_back(sat3);
		//cout << RKSim.m_SS.m_sats[0].m_thrusts[0].m_dvx[0] << "\n";
		for (unsigned long int t = 0; t < RKSim.m_par->maxTimeSteps; t++) {
			if (int(t) % par.print_freq == 0) {
				cout << "Timestep " << t << " of " << RKSim.m_par->maxTimeSteps
						<< "\n";
				cout.flush();
				RKSim.print_pos(pos_stream);
			}
			RKSim.update(t);
		}
	} else if (par.runMode == 3) {
		unsigned init_time =
				chrono::high_resolution_clock::now().time_since_epoch().count();
		en_stream << 0 << ", " << 0 << "\n";
		for (unsigned long int t = 0; t < RKSim.m_par->maxTimeSteps; t++) {
			if (int(t) % par.print_freq == 0) {
				cout << "Timestep " << t << " of " << RKSim.m_par->maxTimeSteps
						<< "\n";
				cout.flush();
				RKSim.print_en(en_stream, init_en);
				en_stream << ", ";
				ESim.print_en(en_stream, init_en);
				en_stream.flush();
				en_stream << "\n";
			}
			RKSim.update(t);
			ESim.update(t);
		}
		unsigned final_time =
				chrono::high_resolution_clock::now().time_since_epoch().count();
		//en_stream << (final_time-init_time) << "\n";
		en_stream.flush();
	} else if (par.runMode == 4){
		for (unsigned long int t = 0; t < RKSim.m_par->maxTimeSteps; t++) {
					if (int(t) % par.print_freq == 0) {
						cout << "Timestep " << t << " of " << RKSim.m_par->maxTimeSteps
								<< "\n";
						cout.flush();
						//RKSim.print_pos(pos_stream);
						en_stream << RKSim.m_SS.kinetic() << ", " << RKSim.m_SS.potential() << "\n";
					}
					RKSim.update(t);
				}
	}
	else {
		cout << "Invalid runMode selected.\n";
	}

	param_of.close();
	en_of.close();
	pos_of.close();
	com_of.close();
	success_of.close();
	cout << "PROGRAM TERMINATING!\n";
	return 0;
}

