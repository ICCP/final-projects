/*
 * CSatellite.h
 *
 *  Created on: Apr 20, 2015
 *      Author: Mark Wittekind and Drew Murray
 */

#ifndef CSATELLITE_H_
#define CSATELLITE_H_

#include "CPlanet.h"
#include "CChromosome.h"
#include <vector>
#include <algorithm>
#include <string>

class CSolarSystem;

using namespace std;

class CSatellite : public CPlanet{
public:
	CSatellite(Par* par, CSolarSystem* ss, SatName name);
	CSatellite(CSatellite rhs, vector<bool> flag);
	CSatellite(const CSatellite &rhs);
	CSatellite();
	virtual ~CSatellite();
	double m_fuel;
	int m_thrustQuant;
	string m_homePlanetName;
	vector<CChromosome> m_thrusts;
	vector<CChromosome> m_thrustsBackup;
	void printThrust(ostream &output);
	void setDynamics(CCoordSet dynamics);
	CSatellite operator*(CSatellite& rhs);
	CSatellite operator=(const CSatellite& rhs);
	SatName m_name;
	unsigned long long unsignedLongLongRand();
	bool m_succFlag;
	CCoordSet m_startDynamics;

private:
	bool m_dynamicsSet;

};

#endif /* CSATELLITE_H_ */
