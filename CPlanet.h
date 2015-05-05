/*
 * CPlanet.h
 *
 *  Created on: Apr 8, 2015
 *      Author: Mark Wittekind and Drew Murray
 */

#ifndef CPLANET_H_
#define CPLANET_H_

#include "CCoordSet.h"
#include "CSVector.h"
#include <string>
#include "Parameters.cpp"
#include <vector>

class CSolarSystem;

class CPlanet {
public:
	CPlanet();
	CPlanet(const CPlanet& rhs);
	CPlanet(float angle, float mass, float peri, float ap, float rad, CCoordSet ccs, string name,
			Par* par, CSolarSystem* ss);
	CPlanet(vector<bool> flagforsatelliteuse);
	virtual ~CPlanet();

	CPlanet operator=(const CPlanet& rhs);
	bool operator==(const CPlanet& rhs);

	Par* m_par;
	CSolarSystem* m_ss;
	float m_mass;
	string m_name;
	float m_angle;
	float m_radius;
	bool m_delFlag;

	float m_aphelion; // furthest from sun
	float m_perihelion;

	float calcTheorVel();
	float calcEscapeVel(float r);
	CSVector calcForce(CSolarSystem* ss);
	float calcEccentricity();

	CCoordSet getDynamics() {
		return m_dynamics;
	}
	CCoordSet* getDynamicsPtr() {
		return &m_dynamics;
	}
	void setDynamics(CCoordSet newSet);

	friend ostream &operator<<(ostream &output, const CPlanet &D) {
		output << D.m_name << ": " << D.m_mass << "kg\n" << D.m_dynamics << "\n";
		return output;
	}
protected:
	CCoordSet m_dynamics;

};

#endif /* CPLANET_H_ */
