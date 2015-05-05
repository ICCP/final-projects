/*
 * CSVector.cpp
 *
 *  Created on: Apr 8, 2015
 *      Author: Mark Wittekind and Drew Murray
 */

#include "CSVector.h"
#include <math.h>

using namespace std;

CSVector::CSVector(float x, float y, float z) {
	m_x = x;
	m_y = y;
	m_z = z;
}
CSVector::CSVector() {
	m_x = 0;
	m_y = 0;
	m_z = 0;
}

CSVector::CSVector(const CSVector& other) {
	m_x = other.m_x;
	m_y = other.m_y;
	m_z = other.m_z;

}

CSVector::~CSVector() {
}

float CSVector::distance(const CSVector other)
{
	double dist = pow(double(m_x)-double(other.m_x),2);
	dist += pow(double(m_y)-double(other.m_y),2);
	dist += pow(double(m_z)-double(other.m_z),2);
	return sqrt(dist);
}

CSVector CSVector::operator=(const CSVector& rhs)
{
	if(this == &rhs)
	        return *this;
	m_x = rhs.m_x;
	m_y = rhs.m_y;
	m_z = rhs.m_z;
	return *this;
}

CSVector CSVector::operator+=(const CSVector& rhs)
{
	if(this == &rhs)
	        return *this;
	m_x += rhs.m_x;
	m_y += rhs.m_y;
	m_z += rhs.m_z;
	return *this;
}

CSVector CSVector::operator-=(const CSVector& rhs)
{
	if(this == &rhs)
	        return *this;
	m_x -= rhs.m_x;
	m_y -= rhs.m_y;
	m_z -= rhs.m_z;
	return *this;
}

CSVector CSVector::operator+(const CSVector& rhs)
{
	float newX = rhs.m_x + m_x;
	float newY = rhs.m_y + m_y;
	float newZ = rhs.m_z + m_z;
	return CSVector(newX,newY,newZ);
}

CSVector CSVector::operator-(CSVector& rhs)
{
	return (*this + rhs * -1.0);
}

CSVector CSVector::operator*(const float rhs)
{
	float newX = m_x * rhs;
	float newY = m_y * rhs;
	float newZ = m_z * rhs;
	return CSVector(newX,newY,newZ);
}
CSVector CSVector::operator*(const double rhs)
{
	double newX = double(m_x) * rhs;
	double newY = double(m_y) * rhs;
	double newZ = double(m_z) * rhs;
	return CSVector(newX,newY,newZ);
}

CSVector CSVector::operator/(const float rhs)
{
	float newX = m_x / rhs;
	float newY = m_y / rhs;
	float newZ = m_z / rhs;
	return CSVector(newX,newY,newZ);
}

CSVector CSVector::operator/(const double rhs)
{
	double newX = double(m_x) / rhs;
	double newY = double(m_y) / rhs;
	double newZ = double(m_z) / rhs;
	return CSVector(newX,newY,newZ);
}

float CSVector::operator*(const CSVector& other)
{//dot product
	float dist = m_x+other.m_x;
	dist += m_y+other.m_y;
	dist += m_z+other.m_z;
	return dist;
}

bool CSVector::operator==(const CSVector& other)
{
	if(m_x==other.m_x && m_y==other.m_y && m_z==other.m_z)
		return true;
	return false;
}

bool CSVector::operator!=(const CSVector& other)
{
	if(m_x!=other.m_x || m_y!=other.m_y || m_z!=other.m_z)
		return true;
	return false;
}

bool CSVector::operator<(CSVector rhs){
	CSVector zero;
	return(this->distance(zero)<rhs.distance(zero));
}
bool CSVector::operator>(CSVector rhs){
	CSVector zero;
	return(this->distance(zero)>rhs.distance(zero));
}
bool CSVector::operator>=(CSVector rhs){
	CSVector zero;
	return(this->distance(zero)>=rhs.distance(zero));
}
bool CSVector::operator<=(CSVector rhs){
	CSVector zero;
	return(this->distance(zero)<=rhs.distance(zero));
}


float CSVector::mag(){
	CSVector zero;
	return distance(zero);
}
float CSVector::magSquared(){
	float dist = pow(m_x,2);
		dist += pow(m_y,2);
		dist += pow(m_z,2);
		return dist;
}

/*
ostream& CSVector::operator<<(ostream &out, const CSVector& c)     //output
{
        out<<"("<<c.m_x<<","<<c.m_y<<","<<c.m_z<<")";
        return out;
}
*/
/*
istream CSVector::&operator>>(istream &in, CSVector &c)     //input
{
        in>>c.m_x;
        in>>c.m_y;
        in>>c.m_z;
        return in;
}
*/
