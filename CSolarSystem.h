/*
 * CSolarSystem.h
 *
 *  Created on: Apr 8, 2015
 *      Author: Mark Wittekind and Drew Murray
 */

#ifndef CSOLARSYSTEM_H_
#define CSOLARSYSTEM_H_

#include "CPlanet.h"
#include "Parameters.cpp"
#include <vector>
#include "CSatellite.h"

using namespace std;

class CSolarSystem {
public:
	CSolarSystem();
	CSolarSystem(Par* par);
	virtual ~CSolarSystem();
	vector<CPlanet> m_planets;
	vector<CSatellite> m_sats;
	Par* m_par;
	vector<double> indPotential();
	vector<double> indKinetic();
	vector<double> totalIndEnergy(vector<double> kinetic, vector<double> potential);
	double potential();
	double kinetic();
	double totalEnergy();
	CSVector calcCOM();
	void adjustMomentum();
	void setSatellites(vector<CSatellite> sats){m_sats = sats;};
	CCoordSet dynamicsFromName(string name);
	float getPlanetRadius(string name);
	CCoordSet getPlanetDynamics(string name);
	CSolarSystem operator=(const CSolarSystem& rhs);
	CSolarSystem(const CSolarSystem& rhs);

	friend ostream &operator<<(ostream &output, const CSolarSystem &D) {
		for (vector<CPlanet>::const_iterator it = D.m_planets.begin();it < D.m_planets.end(); it++) {
			output << *it;
		}
		output << "-\n";
		return output;
	}

};

#endif /* CSOLARSYSTEM_H_ */
