/*
 * CCoordSet.h
 *
 *  Created on: Apr 8, 2015
 *      Author: Mark Wittekind and Drew Murray
 */

#ifndef CCOORDSET_H_
#define CCOORDSET_H_

#include "CSVector.h"

class CCoordSet {
public:
	CCoordSet();
	CCoordSet(const CCoordSet& other);
	CCoordSet(CSVector pos, CSVector vel, CSVector acc);
	CCoordSet operator=(const CCoordSet& rhs);
	virtual ~CCoordSet();
	CSVector m_position;
	CSVector m_velocity;
	CSVector m_acceleration;

	friend ostream &operator<<( ostream &output, const CCoordSet &D ){
		output << "Pos:" << D.m_position << "\nVel:" << D.m_velocity << "\nAcc:"
				<< D.m_acceleration << "\n";
		return output;
		}
};

#endif /* CCOORDSET_H_ */
