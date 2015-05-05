/*
 * CSimulation.cpp
 *
 *  Created on: Apr 8, 2015
 *      Author: Mark Wittekind and Drew Murray
 */

#include "CSimulation.h"

bool fuelLesserCompare(CSatellite lhs, CSatellite rhs) {
	return (lhs.m_fuel < rhs.m_fuel);
}

CSimulation::CSimulation(CSolarSystem solarSystem, Par* par, string name,
		unsigned short coreNum, unsigned int genNum) {
	m_par = par;
	m_name = name;
	m_coreNum = coreNum;
	m_genNum = genNum;
	m_SS = solarSystem;
	CCoordSet tempSet;
	for (unsigned int i = 0; i < m_par->satsPerCore; i++) {
		SatName newSatelliteName;
		newSatelliteName.coreNum = coreNum;
		newSatelliteName.genNum = genNum;
		CSatellite a(m_par, &m_SS, newSatelliteName);
		a.m_name.selfID = a.unsignedLongLongRand();

		tempSet = m_SS.getPlanetDynamics(a.m_homePlanetName);
		tempSet.m_position.m_y = -tempSet.m_position.m_y
				+ m_SS.getPlanetRadius(a.m_homePlanetName);
		a.setDynamics(tempSet);
		a.m_startDynamics = tempSet;

		a.m_thrusts[a.m_thrusts.size() - 1].m_t[0] = 0;

		m_SS.m_sats.push_back(a);
	}
}

CSimulation::CSimulation(CSolarSystem solarSystem, Par* par, string name,
		vector<CSatellite> sats, unsigned short coreNum, unsigned int genNum) {

	m_par = par;
	m_name = name;
	m_coreNum = coreNum;
	m_genNum = genNum;
	m_SS = solarSystem;
	vector<bool> flag;
	CCoordSet tempSet;

	m_SS.m_sats = sats;
}

CSimulation CSimulation::operator=(const CSimulation &rhs) {

	m_name = rhs.m_name;
	m_coreNum = rhs.m_coreNum;
	m_genNum = rhs.m_genNum;
	m_par = rhs.m_par;
	m_succList = rhs.m_succList;
	m_breedList = rhs.m_breedList;
	m_deadList = rhs.m_deadList;
	m_SS = rhs.m_SS;
	return *this;
}

CSimulation::CSimulation(const CSimulation &rhs) {
	m_name = rhs.m_name;
	m_coreNum = rhs.m_coreNum;
	m_genNum = rhs.m_genNum;
	m_par = rhs.m_par;
	m_succList = rhs.m_succList;
	m_breedList = rhs.m_breedList;
	m_deadList = rhs.m_deadList;
	m_SS = rhs.m_SS;
}

CSimulation CSimulation::operator*(const CSimulation &rhs) {
	CSolarSystem newSS(rhs.m_SS);
	CSimulation newSim(rhs);
	newSim.m_SS.m_sats.clear();
	newSim.m_genNum++;
	CCoordSet tempSet;
	vector<CSatellite> loserList;
	cout << "-----------------------------\n";
	for (unsigned int i = 0; i < m_deadList.size(); i++) {
		loserList.push_back(m_deadList[i]);
	}
	for (unsigned int i = 0; i < rhs.m_SS.m_sats.size(); i++) {
		loserList.push_back(rhs.m_SS.m_sats[i]);
	}

	unsigned seed =
			chrono::high_resolution_clock::now().time_since_epoch().count()
					+ m_coreNum;
	default_random_engine gen(seed);

	//fill the success portion
	unsigned int succUses = 0;
	sort(m_succList.begin(), m_succList.end(), fuelLesserCompare);

	if (m_succList.size() >= m_par->breedingNum - m_par->randomBreedingNum) {
		//enough successes
		succUses = m_par->breedingNum - m_par->randomBreedingNum;
	} else {
		//not enough successes
		succUses = m_succList.size();
	}
	if(m_succList.size() < .1){
		CSolarSystem newSS2(rhs.m_SS);
		CSimulation newSim2(newSS2, rhs.m_par, rhs.m_name,
		666, rhs.m_genNum+1);
		return newSim2;
	}
	else{
		for (unsigned int i = 0; i < succUses; i++) {
			m_breedList.push_back(m_succList[i]);
		}
		for (unsigned int i = succUses; i < m_succList.size(); i++) {
			loserList.push_back(m_succList[i]);
		}

		uniform_int_distribution<int> dist(0, loserList.size() - 1);
		for (unsigned int i = 0; i < m_par->randomBreedingNum; i++) {
			int r = dist(gen);
			m_breedList.push_back(loserList[r]);
		}
		uniform_int_distribution<int> dist2(0, m_breedList.size() - 1);
		for (unsigned int x = 0; x < (m_par->satsPerCore-m_par->randSatsPerGen); x++) {

			newSim.m_SS.m_sats.push_back(
					m_breedList[dist2(gen)] * m_breedList[dist2(gen)]);
		}
		for(unsigned int x = 0; x < m_par->randSatsPerGen; x++){
			CCoordSet tempSet;
			SatName name;
			name.coreNum = 666;
			name.genNum = 1+rhs.m_genNum;
			name.momID = 0;
			name.dadID = 0;
			CSatellite a(m_par, &newSS, name);
			tempSet = m_SS.getPlanetDynamics(a.m_homePlanetName);
			tempSet.m_position.m_y = -tempSet.m_position.m_y
					+ m_SS.getPlanetRadius(a.m_homePlanetName);
			tempSet = m_SS.getPlanetDynamics(a.m_homePlanetName);
					tempSet.m_position.m_y = -tempSet.m_position.m_y
							+ m_SS.getPlanetRadius(a.m_homePlanetName);
					a.setDynamics(tempSet);
					a.m_startDynamics = tempSet;

					a.m_thrusts[a.m_thrusts.size() - 1].m_t[0] = 0;
			newSim.m_SS.m_sats.push_back(a);
		}
		return newSim;
	}
}

CSimulation::~CSimulation() {
}

void CSimulation::Euler(CSVector &pos_i, CSVector &vel_i, string name) {
	float dt = m_par->dt;
	CSVector new_pos = pos_i + vel_i * dt;
	bool temp = false;
	vel_i += CalcA(pos_i, name, temp, temp, dt) * dt;
	pos_i = new_pos;
}

void CSimulation::RK(CSVector &pos_i, CSVector &vel_i, bool &delFlag,
		bool &succFlag, string name) {
	// http://doswa.com/2009/01/02/fourth-order-runge-kutta-numerical-integration.html

	float dt = m_par->dt;
	if (pos_i.magSquared() > this->m_SS.m_par->boundaryDistanceSquared) {
		cout << "OutofBounds for satellite ";
		delFlag = true;
	} else {
		CSVector pos[4], vel[4], acc[4];
		pos[0] = pos_i;
		vel[0] = vel_i;
		acc[0] = CalcA(pos[0], name, delFlag, succFlag, dt / 2.0);

		pos[1] = pos_i + vel[0] * (dt * 0.5);
		vel[1] = vel_i + acc[0] * (dt * 0.5);
		acc[1] = CalcA(pos[1], name, delFlag, succFlag, dt / 2.0);

		pos[2] = pos_i + vel[1] * (dt * 0.5);
		vel[2] = vel_i + acc[1] * (dt * 0.5);
		acc[2] = CalcA(pos[2], name, delFlag, succFlag, dt / 2.0);

		pos[3] = pos_i + vel[2] * dt;
		vel[3] = vel_i + acc[2] * dt;
		acc[3] = CalcA(pos[3], name, delFlag, succFlag, dt);

		pos_i += (vel[0] + vel[1] * 2.0 + vel[2] * 2.0 + vel[3]) * (dt / 6.0);
		vel_i += (acc[0] + acc[1] * 2.0 + acc[2] * 2.0 + acc[3]) * (dt / 6.0);
	}
	return;
}

CSVector CSimulation::CalcA(CSVector pos, string name, bool &delFlag,
		bool &succFlag, float dt) {
	CSVector sum, r_vec;
	double dist;

	for (vector<CPlanet>::iterator it = m_SS.m_planets.begin();
			it < m_SS.m_planets.end(); it++) {
		dist = pos.distance(it->getDynamics().m_position);

		if (it->m_name != name) {
			if (dist < it->m_radius) {
				cout << "Collision for satellite ";
				cout.flush();
				delFlag = true;
				if (m_par->destinationPlanetName == it->m_name) {
					succFlag = true;
				}
			} else {
				double dist_cubed = pow(dist, 3);
				if (dist_cubed == 0)
					dist_cubed = pow(10, -200);

				if (it->m_name != name) {
					r_vec = it->getDynamics().m_position + pos * -1.0;
					sum += (r_vec / dist_cubed) * (it->m_mass);
				}
			}
			if (m_par->destinationPlanetName == it->m_name
					&& dist < it->m_radius * m_par->goalDistanceMultiplier) {
				cout << pos << " : " << it->getDynamics().m_position << "\n";
				cout << "Noncollision success for satellite ";
				succFlag = true;
				delFlag = true;
				cout.flush();
			}
		}

		//else if (name == "Pluto")cout<<"v: "<<it->getDynamics().m_velocity << endl;
	}

	return (sum * this->m_SS.m_par->G);
}

void CSimulation::print_pos(ostream &pos_stream) {
	for (vector<CPlanet>::iterator it = m_SS.m_planets.begin();
			it < m_SS.m_planets.end(); it++) {
		pos_stream << (*(it->getDynamicsPtr())).m_position << endl;
	}
	for (vector<CSatellite>::iterator it = m_SS.m_sats.begin();
			it < m_SS.m_sats.end(); it++) {
		pos_stream << (*(it->getDynamicsPtr())).m_position << endl;
	}
	pos_stream << "-\n";
	pos_stream.flush();

}

void CSimulation::print_en(ostream &en_stream, double init_en) {

	en_stream << (m_SS.totalEnergy() - init_en) / init_en;

}

void CSimulation::update(unsigned long currStep) {
	//cout << "Using integration method " << m_par->integration_method << endl;
	CSolarSystem newSS(m_SS);

	for (vector<CSatellite>::iterator it = newSS.m_sats.begin();
			it < newSS.m_sats.end(); it++) {
		if (it->m_thrusts.size() > 0
				&& it->m_thrusts[it->m_thrusts.size() - 1].m_t[0] - currStep
						< .5) {
			CSVector a(it->m_thrusts[it->m_thrusts.size() - 1].m_dvx[0],
					it->m_thrusts[it->m_thrusts.size() - 1].m_dvy[0],
					it->m_thrusts[it->m_thrusts.size() - 1].m_dvz[0]);
			it->getDynamicsPtr()->m_velocity += a;
			it->m_fuel += a.magSquared();
			it->m_thrusts.pop_back();
		}
	}
	bool temp = false;
	for (vector<CPlanet>::iterator it = newSS.m_planets.begin();
			it < newSS.m_planets.end(); it++) {
		if (m_par->integration_method == 4) {
			RK(it->getDynamicsPtr()->m_position,
					it->getDynamicsPtr()->m_velocity, it->m_delFlag, temp,
					it->m_name);
		} else if (m_par->integration_method == 1) {
			Euler(it->getDynamicsPtr()->m_position,
					it->getDynamicsPtr()->m_velocity, it->m_name);
		}

		else
			cout << "ERROR, integration method set to "
					<< m_par->integration_method << endl;
	}
	for (vector<CSatellite>::iterator it = newSS.m_sats.begin();
			it < newSS.m_sats.end(); it++) {
		if (m_par->integration_method == 4) {
			RK(it->getDynamicsPtr()->m_position,
					it->getDynamicsPtr()->m_velocity, it->m_delFlag,
					it->m_succFlag, "Satellite");
		} else if (m_par->integration_method == 1) {
			Euler(it->getDynamicsPtr()->m_position,
					it->getDynamicsPtr()->m_velocity, "Satellite");
		}

		else
			cout << "ERROR, integration method set to "
					<< m_par->integration_method << endl;
	}
	unsigned int i = 0;

	//remove all flagged sats
	while (newSS.m_sats.begin() + i != newSS.m_sats.end()) {
		if (newSS.m_sats[i].m_delFlag) {
			if (newSS.m_sats[i].m_succFlag) {
				m_succList.push_back(newSS.m_sats[i]);
			} else {
				m_deadList.push_back(newSS.m_sats[i]);
			}
			cout << newSS.m_sats[i].m_name.selfID << " at timestep " << currStep
					<< endl;
			newSS.m_sats.erase(newSS.m_sats.begin() + i);
		} else {
			i++;
		}
	}

	m_SS = newSS;
}

void CSimulation::successPrint(ostream& output) {

	for (vector<CSatellite>::iterator it = m_succList.begin();
			it < m_succList.end(); it++) {
			output << "   Mom" << " : " << it->m_name.momID << endl;
			output << "   Dad" << " : " << it->m_name.dadID << endl;
			output << "   SelfID: " << it->m_name.selfID << endl;
			output << "   Generation: " << it->m_name.genNum << endl;
			output << "   Core: " << it->m_name.coreNum << endl;
			output << "   Sat: " << it->m_name.selfID << endl;
			it->printThrust(output);
	}
	output.flush();
}
void CSimulation::genePrint(string simName) {
	cout << "Simulation Genetic Summary for sim : " << simName << endl;
	unsigned int n = 0;
	for (vector<CSatellite>::iterator it = m_SS.m_sats.begin();
			it < m_SS.m_sats.end(); it++) {
		cout << "   Mom" << " : " << it->m_name.momID << endl;
		cout << "   Dad" << " : " << it->m_name.dadID << endl;
		cout << "   Sat" << n << ": " << it->m_name.selfID << endl;
		it->printThrust(cout);
		n++;
	}
}
