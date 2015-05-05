/*
 * CSolarSystem.cpp
 *
 *  Created on: Apr 8, 2015
 *      Author: Mark Wittekind and Drew Murray
 */

#include "CSolarSystem.h"
#include <math.h>

CSolarSystem::CSolarSystem() {
	m_par = NULL;
}


CSolarSystem::CSolarSystem(const CSolarSystem& rhs) {

	m_par = rhs.m_par;
	m_planets = rhs.m_planets;
	m_sats = rhs.m_sats;

}

CSolarSystem CSolarSystem::operator=(const CSolarSystem& rhs)
{
	if(this == &rhs)
	        return *this;
	m_par = rhs.m_par;
	m_planets = rhs.m_planets;
	m_sats = rhs.m_sats;
	return *this;
}

CSolarSystem::CSolarSystem(Par* par) {


	m_par = par;
	float sm = m_par->SolarMass;
	CSVector zero;

	{
	//Sun
	CCoordSet ccs_Sun(zero,zero,zero);
	CPlanet Sun(0,sm,0,0,696342000, ccs_Sun, string("Sun"),m_par,this);
	m_planets.push_back(Sun);
	}{
	//Mercury
	CCoordSet ccs_Mercury(zero,zero,zero);
	CPlanet Mercury(0.110653875,sm*1.66012*pow(10,-7),46001200000,69816900000,4879000/2.0, ccs_Mercury, string("Mercury"),m_par,this);
	m_planets.push_back(Mercury);
	}{
	//Venus
	CCoordSet ccs_Venus(zero,zero,zero);;
	CPlanet Venus(0.0382227106,sm*2.44781*pow(10,-6),107477000000,108939000000,12104000/2.0, ccs_Venus, string("Venus"),m_par,this);
	m_planets.push_back(Venus);
	}{
	//Earth
	CCoordSet ccs_Earth(zero,zero,zero);;
	CPlanet Earth(0.0274016693,sm*3.003467*pow(10,-6),147095000000,151930000000, 12742000/2.0, ccs_Earth, string("Earth"),m_par,this);
	m_planets.push_back(Earth);
	}{
	//Mars
	CCoordSet ccs_Mars(zero,zero,zero);;
	CPlanet Mars(0.0291469985,sm*3.22713*pow(10,-7),206700000000,249200000000,6779000/2.0, ccs_Mars, string("Mars"),m_par,this);
	m_planets.push_back(Mars);
	}{
	//Jupiter
	CCoordSet ccs_Jupiter(zero,zero,zero);;
	CPlanet Jupiter(0.00558505361,sm*9.5458*pow(10,-4),740573600000,816520800000, 69911000, ccs_Jupiter, string("Jupiter"),m_par,this);
	m_planets.push_back(Jupiter);
	}{
	CCoordSet ccs_Saturn(zero,zero,zero);;
	CPlanet Saturn(0.016231562,sm*2.85812*pow(10,-4),1353572956000,1513325783000,116464000/2.0, ccs_Saturn, string("Saturn"),m_par,this);
	m_planets.push_back(Saturn);
	}{
	CCoordSet ccs_Uranus(zero,zero,zero);;
	CPlanet Uranus(0.0178023584,sm*4.36576*pow(10,-5),2735118100000,3006224700000,50724000/2.0, ccs_Uranus, string("Uranus"),m_par,this);
	m_planets.push_back(Uranus);
	}{
	CCoordSet ccs_Neptune(zero,zero,zero);;
	CPlanet Neptune(0.0125663706,sm*5.15028*pow(10,-5),4459504400000,4537580900000,49244000/2.0, ccs_Neptune, string("Neptune"),m_par,this);
	m_planets.push_back(Neptune);
	}{
	CCoordSet ccs_Pluto(zero,zero,zero);;
	CPlanet Pluto(0.296705973,sm*6.583*pow(10,-9),4437000000000,7311000000000, 2390000/2.0, ccs_Pluto, string("Pluto"),m_par,this);
	m_planets.push_back(Pluto);
	}
	//cout << "m_aphelion*tan(0.296705973) = " << Pluto.m_aphelion*tan(0.296705973) << "\n";
	for (vector<CPlanet>::iterator it = m_planets.begin();it < m_planets.end(); it++) {
		it->getDynamicsPtr()->m_position.m_x=it->m_aphelion*cos(it->m_angle);
		it->getDynamicsPtr()->m_position.m_z=it->m_aphelion*sin(it->m_angle);
		it->getDynamicsPtr()->m_velocity.m_y=it->calcTheorVel();
	}
	adjustMomentum(); // Set total momentum to zero

}

CSolarSystem::~CSolarSystem() {
}

vector<double> CSolarSystem::indPotential(){
	short n = m_planets.size();
	short i = 0;
	vector<double> energies (n,1);
	for(i = 0; i < n; i++){
		for(short j = i + 1; j < n; j++){
			float r = m_planets[i].getDynamics().m_position.distance(m_planets[j].getDynamics().m_position);
			double total = -m_par->G*double(m_planets[i].m_mass)*double(m_planets[j].m_mass)/double(r);
			energies[i] += total;
			energies[j] += total;
		}
	}
	return energies;
}

vector<double> CSolarSystem::indKinetic(){
	vector<double> energies;
	for (vector<CPlanet>::iterator it = m_planets.begin();it < m_planets.end(); it++) {
		energies.push_back(double(.5)*it->m_mass*it->getDynamics().m_velocity.magSquared());
	}
	return energies;
}

vector<double> CSolarSystem::totalIndEnergy(vector<double> kinetic, vector<double> potential){
	vector<double> energies;
	for (unsigned int i = 0; i < kinetic.size(); i++){
		energies.push_back(kinetic[i] + potential[i]);
	}
	return energies;
}

double CSolarSystem::potential(){
	short n = m_planets.size();
	short i = 0;
	double total=0;
	for(i = 0; i < n; i++){
		for(short j = i + 1; j < n; j++){
			float r = m_planets[i].getDynamics().m_position.distance(m_planets[j].getDynamics().m_position);
			double subtotal = -m_par->G*double(m_planets[i].m_mass)*double(m_planets[j].m_mass)/double(r);
			total += subtotal;
		}
	}
	return total;
}

double CSolarSystem::kinetic(){
	double total=0;
	for (vector<CPlanet>::iterator it = m_planets.begin();it < m_planets.end(); it++) {
		total+=(double(.5)*it->m_mass*it->getDynamics().m_velocity.magSquared());
	}
	return total;
}

double CSolarSystem::totalEnergy(){
	return kinetic()+potential();
}

CSVector CSolarSystem::calcCOM(){
	CSVector com;
	double x = 0;
	double y = 0;
	double z = 0;
	double totalMass = 0;

	for (vector<CPlanet>::iterator it = m_planets.begin();it < m_planets.end(); it++) {
		x+=double((it->getDynamics().m_position.m_x))*(it->m_mass);
		y+=double((it->getDynamics().m_position.m_y))*(it->m_mass);
		z+=double((it->getDynamics().m_position.m_z))*(it->m_mass);
		totalMass+=it->m_mass;
	}
//	cout << "Total Mass: " << totalMass << "\n";
//	cout << "pre div x: " << x << "\n";
//	cout << "pre div y: " << y << "\n";
//	cout << "pre div z: " << z << "\n";
	x /= totalMass;
	y /= totalMass;
	z /= totalMass;
//	cout << "x: " << x << "\n";
//	cout << "y: " << y << "\n";
//	cout << "z: " << z << "\n";
	com.m_x = x;
	com.m_y = y;
	com.m_z = z;

	return (com);
}

void CSolarSystem::adjustMomentum(){
	double x=0;
	double y=0;
	double z=0;
	double M=0;
	for (vector<CPlanet>::iterator it = m_planets.begin();it < m_planets.end(); it++) {
		x+=(it->getDynamics().m_velocity.m_x*it->m_mass);
		y+=(it->getDynamics().m_velocity.m_y*it->m_mass);
		z+=(it->getDynamics().m_velocity.m_z*it->m_mass);
		M+=it->m_mass;
	}
	//cout << "momentum: " << x << " " << y << " " << z << "\n";


	for (vector<CPlanet>::iterator it = m_planets.begin();it < m_planets.end(); it++) {
			it->getDynamicsPtr()->m_velocity.m_x-=x/M;
			it->getDynamicsPtr()->m_velocity.m_y-=y/M;
			it->getDynamicsPtr()->m_velocity.m_z-=z/M;
		}


/*
	x = y = z = 0;
	for (vector<CPlanet>::iterator it = m_planets.begin();it < m_planets.end(); it++) {
			x+=(it->getDynamics().m_velocity.m_x*it->m_mass);
			y+=(it->getDynamics().m_velocity.m_y*it->m_mass);
			z+=(it->getDynamics().m_velocity.m_z*it->m_mass);
		}


		cout << "momentum: " << x << " " << y << " " << z << "\n";
		*/

	return;
}

CCoordSet CSolarSystem::getPlanetDynamics(string name){
	unsigned short i = 0;
	while(m_planets[i].m_name != name){
		i++;
	}
	return m_planets[i].getDynamics();
}

float CSolarSystem::getPlanetRadius(string name){
	unsigned short i = 0;
	while(m_planets[i].m_name != name){
		i++;
	}
	return m_planets[i].m_radius;
}
