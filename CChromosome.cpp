/*
 * CChromosome.cpp
 *
 *  Created on: Apr 24, 2015
 *      Author: Mark Wittekind and Drew Murray
 */

#include "CChromosome.h"
CChromosome::CChromosome(){
	cout << "Default constructor of CChrom used\n";
	m_par = NULL;
	m_coreNum = 0;
}


CChromosome::CChromosome(Par* par, unsigned short coreNum) {
	for(unsigned short i = 0; i<3;i++){
		m_dvx.push_back(0);
		m_dvy.push_back(0);
		m_dvz.push_back(0);
		m_t.push_back(0);
	}
	m_coreNum = coreNum;


	unsigned seed = chrono::high_resolution_clock::now().time_since_epoch().count()+m_coreNum;
	default_random_engine gen(seed);
	uniform_real_distribution<float> distThrust(-1.0*(par->maxThrusterVel),par->maxThrusterVel);
	uniform_int_distribution<int> distTime(0,par->maxTimeSteps);

	m_par = par;
	float x;
	m_dvx[1] = -1.0*(par->maxThrusterVel);
	m_dvx[2] = (par->maxThrusterVel);
	m_dvx[0] = x = distThrust(gen);

	m_dvy[1] = -1.0*(par->maxThrusterVel);
	m_dvy[2] = (par->maxThrusterVel);
	m_dvy[0] = distThrust(gen);

	m_dvz[1] = -1.0*(par->maxThrusterVel);
	m_dvz[2] = (par->maxThrusterVel);
	m_dvz[0] = distThrust(gen);

	m_t[1] = 0;
	m_t[2] = float((par->maxTimeSteps));
	m_t[0] = float(distTime(gen));
}

CChromosome::CChromosome(const CChromosome &rhs) {

	m_par = rhs.m_par;
	m_coreNum = rhs.m_coreNum;

	for(unsigned short i = 0; i<3;i++){
		m_dvx.push_back(0);
		m_dvy.push_back(0);
		m_dvz.push_back(0);
		m_t.push_back(0);
	}

	m_dvx[1] = rhs.m_dvx[1];
	m_dvx[2] = rhs.m_dvx[2];
	m_dvx[0] = rhs.m_dvx[0];

	m_dvy[1] = rhs.m_dvy[1];
	m_dvy[2] = rhs.m_dvy[2];
	m_dvy[0] = rhs.m_dvy[0];

	m_dvz[1] = rhs.m_dvz[1];
	m_dvz[2] = rhs.m_dvz[2];
	m_dvz[0] = rhs.m_dvz[0];

	m_t[1] = rhs.m_t[1];
	m_t[2] = rhs.m_t[2];
	m_t[0] = rhs.m_t[0];
}


CChromosome::~CChromosome() {
}


CChromosome CChromosome::operator=(const CChromosome& rhs){

	m_par->G = rhs.m_par->G;
	m_par = rhs.m_par;

	m_dvx[1] = rhs.m_dvx[1];
	m_dvx[2] = rhs.m_dvx[2];
	m_dvx[0] = rhs.m_dvx[0];

	m_dvy[1] = rhs.m_dvy[1];
	m_dvy[2] = rhs.m_dvy[2];
	m_dvy[0] = rhs.m_dvy[0];

	m_dvz[1] = rhs.m_dvz[1];
	m_dvz[2] = rhs.m_dvz[2];
	m_dvz[0] = rhs.m_dvz[0];

	m_t[1] = rhs.m_t[1];
	m_t[2] = rhs.m_t[2];
	m_t[0] = rhs.m_t[0];


	return rhs;
}

void CChromosome::printChrom(ostream &output){
	output << m_dvx[0] << " " << m_dvy[0] << " " << m_dvz[0] << " " << m_t[0] << "\n";
}

CChromosome CChromosome::operator*(CChromosome& rhs)
{
	CChromosome newChrom(*this);
	uniform_real_distribution<float> randOne(0,1);

	for(unsigned short i = 0;i<newChrom.chromSize();i++){
		unsigned seed = chrono::high_resolution_clock::now().time_since_epoch().count()+i+m_coreNum;
		default_random_engine gen(seed);
		uniform_real_distribution<float> dist((*newChrom[i])[1],(*newChrom[i])[2]);
		float r = randOne(gen);
		if (r < m_par->mutateChance){
			//cout << "Mutating\n";
			(*newChrom[i])[0] = dist(gen);
		}
		else if (r < (m_par->mutateChance+1.0)/2.0){
			(*newChrom[i])[0] = (*(*this)[i])[0];
			//cout << "Using lhs\n";
		}
		else {
			(*newChrom[i])[0] = (*rhs[i])[0];
			//cout << "Using rhs\n";
		}
	}
	return newChrom;
}

vector<float> *CChromosome::operator[](int nIndex)
{
	if (nIndex==0)return &m_dvx;
	else if (nIndex==1)return &m_dvy;
	else if (nIndex==2)return &m_dvz;
	else if (nIndex==3)return &m_t;
	else{
		cout << "ERROR CChromosome::operator[] has invalid index " << nIndex << "\n";\
		return  NULL;
	}

}


bool CChromosome::operator==(const CChromosome& rhs)
{
	if(m_dvx[0]==rhs.m_dvx[0] && m_dvy[0]==rhs.m_dvy[0] && m_dvz[0]==rhs.m_dvz[0] && m_t[0]==rhs.m_t[0])
		return true;
	return false;
}

bool CChromosome::operator!=(const CChromosome& rhs)
{
	if(m_dvx[0]!=rhs.m_dvx[0] || m_dvy[0]!=rhs.m_dvy[0] || m_dvz[0]!=rhs.m_dvz[0] || m_t[0]!=rhs.m_t[0])
		return true;
	return false;
}

bool CChromosome::operator<(const CChromosome& rhs){
	return(this->m_t[0]<rhs.m_t[0]);
}
bool CChromosome::operator>(const CChromosome& rhs){
	return(this->m_t[0]>rhs.m_t[0]);
}
bool CChromosome::operator>=(const CChromosome& rhs){
	return(this->m_t[0]>=rhs.m_t[0]);
}
bool CChromosome::operator<=(const CChromosome& rhs){
	return(this->m_t[0]<=rhs.m_t[0]);
}
