/*
 * CSimulation.h
 *
 *  Created on: Apr 8, 2015
 *      Author: Mark Wittekind and Drew Murray
 */

#ifndef CSIMULATION_H_
#define CSIMULATION_H_

#include "CPlanet.h"
#include "Parameters.cpp"
#include "CSolarSystem.h"
#include <fstream>
#include <vector>
#include <fstream>

class CSimulation {
public:
	CSimulation(CSolarSystem solarSystem, Par *par, string name, unsigned short coreNum, unsigned int genNum);
	CSimulation(CSolarSystem solarSystem, Par* par, string name, vector<CSatellite> sats, unsigned short coreNum, unsigned int genNum);
	CSimulation(string name, Par* par, const CSolarSystem solarSystem, CSimulation b, unsigned short coreNum);
	CSimulation (const CSimulation &rhs);

	virtual ~CSimulation();
	void RK(CSVector &pos_i, CSVector &vel_i, bool &delFlag, bool &succFlag, string name);
	void Euler(CSVector &pos_i, CSVector &vel_i, string name);
	CSVector CalcDv(float dt, vector<CPlanet> planets, unsigned short i);
	CSVector CalcA(CSVector pos, string name, bool &delFlag, bool &succFlag, float dt);
	void print_pos(ostream &pos_stream);
	void print_en(ostream &en_stream, double init_en);
	void update(unsigned long currStep);
	string m_name;
	unsigned short m_coreNum;
	unsigned int m_genNum;
	Par* m_par;
	void genePrint(string simName);
	void successPrint(ostream& output);
	vector<CSatellite> m_succList;
	vector<CSatellite> m_breedList;
	vector<CSatellite> m_deadList;
	CSimulation operator=(const CSimulation &rhs);
	CSimulation operator*(const CSimulation &rhs);

	CSolarSystem m_SS;
};

#endif /* CSIMULATION_H_ */
