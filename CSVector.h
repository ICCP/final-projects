/*
 * CSVector.h
 *
 *  Created on: Apr 8, 2015
 *      Author: Mark Wittekind and Drew Murray
 */

#ifndef CSVECTOR_H_
#define CSVECTOR_H_

#include <iostream>
using namespace std;

class CSVector {
public:
	CSVector(float x, float y, float z);
	CSVector(const CSVector& other);
	CSVector();
	virtual ~CSVector();
	float distance(const CSVector other);
	float mag();
	float magSquared();
	float m_x;
	float m_y;
	float m_z;

	friend ostream &operator<<( ostream &output, const CSVector &D ){
		output<<D.m_x<<","<<D.m_y<<","<<D.m_z;
		return output;
	}

	CSVector operator+(const CSVector& rhs);
	CSVector operator-(CSVector& rhs);
	CSVector operator=(const CSVector& rhs);
	CSVector operator-=(const CSVector& rhs);
	CSVector operator+=(const CSVector& rhs);
	CSVector operator*(const float rhs);
	CSVector operator*(const double rhs);
	CSVector operator/(const float rhs);
	CSVector operator/(const double rhs);
	bool operator==(const CSVector& other);
	bool operator!=(const CSVector& other);
	float operator*(const CSVector& rhs);
	bool operator<(CSVector rhs);
	bool operator>(CSVector rhs);
	bool operator>=(CSVector rhs);
	bool operator<=(CSVector rhs);


};





#endif /* CSVECTOR_H_ */
