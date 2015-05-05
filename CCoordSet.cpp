/*
 * CCoordSet.cpp
 *
 *  Created on: Apr 8, 2015
 *      Author: Mark Wittekind and Drew Murray
 */

#include "CCoordSet.h"

CCoordSet::CCoordSet() {
	CSVector zero;
	m_position =zero;
	m_velocity =zero;
	m_acceleration =zero;

}

CCoordSet::CCoordSet(CSVector pos, CSVector vel, CSVector acc) {
	m_position = pos;
	m_velocity = vel;
	m_acceleration = acc;
}

CCoordSet CCoordSet::operator=(const CCoordSet& rhs)
{
	if(this == &rhs)
	        return *this;
	m_position = rhs.m_position;
	m_velocity = rhs.m_velocity;
	m_acceleration = rhs.m_acceleration;
	return *this;
}

CCoordSet::CCoordSet(const CCoordSet& other) {
	m_position = other.m_position;
	m_velocity = other.m_velocity;
	m_acceleration = other.m_acceleration;

}

CCoordSet::~CCoordSet() {
}

