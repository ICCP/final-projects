/*
 * CSatellite.cpp
 *
 *  Created on: Apr 20, 2015
 *      Author: Mark Wittekind and Drew Murray
 */

#include "CSatellite.h"
using namespace std;

unsigned long long CSatellite::unsignedLongLongRand()
{
	//from stack overflow
	//http://stackoverflow.com/questions/3665257/generate-random-long-number-in-c
	unsigned seed = chrono::high_resolution_clock::now().time_since_epoch().count()+m_name.coreNum;
	srand(seed);
    unsigned long long rand1 = abs(rand());
    unsigned long long rand2 = abs(rand());
    rand1 = rand1 << (sizeof(int)*8);
    unsigned long long randULL = (rand1 | rand2);
    return randULL;
}

bool thrustGreaterCompare(CChromosome lhs, CChromosome rhs){
		return(lhs.m_t>rhs.m_t);
}


CSatellite::CSatellite(Par* par, CSolarSystem* ss, SatName name) : CPlanet() {

	m_thrustQuant = par->thrustQuant;
	for(int i = 0; i < m_thrustQuant; i++){
		CChromosome a(par, name.coreNum);
		m_thrusts.push_back(a);
	}
	sort(m_thrusts.begin(), m_thrusts.end(),thrustGreaterCompare);
	m_thrusts[m_thrusts.size()-1].m_t[0]=0;
	m_thrustsBackup = m_thrusts;

	m_homePlanetName = par->homePlanetName;
	m_fuel = 0;
	m_ss = ss;
	m_mass = par->satMass;
	m_perihelion = -1;
	m_aphelion = -1;
	CCoordSet m_dynamics;
	CCoordSet m_startDynamics;
	m_name = name;
	m_dynamicsSet = false;
	m_angle = -1;
	m_delFlag = false;
	m_radius = 0;
	m_succFlag = false;

}


CSatellite::CSatellite(const CSatellite &rhs){
	m_thrustQuant = rhs.m_thrustQuant;
	m_thrusts = rhs.m_thrusts;
	m_homePlanetName = rhs.m_homePlanetName;
	m_fuel = rhs.m_fuel;
	m_ss = rhs.m_ss;
	m_mass = rhs.m_mass;
	m_perihelion = rhs.m_perihelion;
	m_aphelion = rhs.m_aphelion;
	m_dynamics = rhs.m_dynamics;
	m_name = rhs.m_name;
	m_dynamicsSet = rhs.m_dynamicsSet;
	m_angle = rhs.m_angle;
	m_delFlag = false;
	m_radius = rhs.m_radius;
	m_startDynamics = rhs.m_startDynamics;
	m_succFlag = false;
	m_par = rhs.m_par;
	m_thrustsBackup = rhs.m_thrustsBackup;
}


CSatellite::CSatellite(CSatellite rhs, vector<bool> flag) {

	m_thrustQuant = rhs.m_thrustQuant;
	m_thrusts = rhs.m_thrusts;
	m_homePlanetName = rhs.m_homePlanetName;
	m_fuel = rhs.m_fuel;
	m_ss = rhs.m_ss;
	m_mass = rhs.m_mass;
	m_perihelion = rhs.m_perihelion;
	m_aphelion = rhs.m_aphelion;
	m_dynamics = rhs.m_dynamics;
	m_name = rhs.m_name;
	m_dynamicsSet = rhs.m_dynamicsSet;
	m_angle = rhs.m_angle;
	m_delFlag = false;
	m_radius = rhs.m_radius;
	m_startDynamics = rhs.m_startDynamics;
	m_succFlag = false;
	m_thrustsBackup = rhs.m_thrustsBackup;
	m_par = rhs.m_par;
}

CSatellite CSatellite::operator*(CSatellite& rhs)
{
	if(!m_dynamicsSet){
		cout << "WARNING: m_dynamics not set for satellite: " << m_name.selfID << ".\n";
	}
	vector<bool> flag;
	CSatellite newSat(rhs, flag);
	newSat.m_thrusts.clear();
	for(unsigned int i = 0; i<rhs.m_thrustsBackup.size();i++){
		newSat.m_thrusts.push_back(rhs.m_thrustsBackup[i]*(this->m_thrustsBackup[i]));
	}
	sort(newSat.m_thrusts.begin(), newSat.m_thrusts.end(),thrustGreaterCompare);
	newSat.m_thrusts[newSat.m_thrusts.size()-1].m_t[0]=0;
	newSat.m_thrustsBackup=newSat.m_thrusts;
	newSat.m_name.genNum++;
	newSat.m_name.momID = m_name.selfID;
	newSat.m_name.dadID = rhs.m_name.selfID;
	newSat.m_name.selfID = unsignedLongLongRand();


	newSat.setDynamics(m_startDynamics);
	m_fuel = 0;
	m_delFlag = false;
	m_succFlag = false;

	return newSat;
}

CSatellite CSatellite::operator=(const CSatellite& rhs)
{
	if(!m_dynamicsSet){
			cout << "WARNING: m_dynamics not set for satellite: " << m_name.selfID << ".\n";
		}
	if(this == &rhs)
	        return *this;
	m_fuel = rhs.m_fuel;
	m_thrustQuant = rhs.m_thrustQuant;
	m_homePlanetName = rhs.m_homePlanetName;
	m_thrusts = rhs.m_thrusts;
	m_ss = rhs.m_ss;
	m_mass = rhs.m_mass;
	m_perihelion = rhs.m_perihelion;
	m_aphelion = rhs.m_aphelion;
	m_dynamics = rhs.m_dynamics;
	m_name = rhs.m_name;
	m_dynamicsSet = rhs.m_dynamicsSet;
	m_startDynamics = rhs.m_startDynamics;
	m_angle = rhs.m_angle;
	m_delFlag = false;
	m_radius = rhs.m_radius;
	m_succFlag = false;
	m_par = rhs.m_par;
	m_thrustsBackup = rhs.m_thrustsBackup;

	return *this;
}

CSatellite::CSatellite() {
	m_fuel = -1;
	m_thrustQuant = -1;
	m_dynamicsSet = false;
	m_succFlag = false;
	cout << "ERROR: Default constructor of CSatellite called.\n";
}

CSatellite::~CSatellite() {
}

void CSatellite::printThrust(ostream &output){
	if(!m_dynamicsSet){
			cout << "WARNING: m_dynamics not set for satellite: " << m_name.selfID << ".\n";
		}
	for(int i = 0; i < m_thrustQuant; i++){
		m_thrustsBackup[i].printChrom(output);
	}
	output << "\n";
}

void CSatellite::setDynamics(CCoordSet dynamics){
	m_dynamics = dynamics;
	m_dynamicsSet = true;
	return;
}


