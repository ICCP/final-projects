/*
 * CChromosome.h
 *
 *  Created on: Apr 24, 2015
 *      Author: Mark Wittekind and Drew Murray
 */

#ifndef CTRAIT_H_
#define CTRAIT_H_

#include "Parameters.cpp"
#include <random>
#include <chrono>

class CChromosome {
public:
	CChromosome(Par* par, unsigned short coreNum);
	CChromosome(const CChromosome &rhs);
	CChromosome();
	virtual ~CChromosome();
	Par* m_par;

	vector<float> m_dvx;
	vector<float> m_dvy;
	vector<float> m_dvz;
	vector<float> m_t;
	unsigned short m_coreNum;
	int chromSize() const {
		return 4;
	}
	void printChrom(ostream &output);

	bool operator==(const CChromosome& rhs);
	bool operator!=(const CChromosome& rhs);
	CChromosome operator*(CChromosome& rhs);
	bool operator<(const CChromosome& rhs);
	bool operator>(const CChromosome& rhs);
	bool operator>=(const CChromosome& rhs);
	bool operator<=(const CChromosome& rhs);
	CChromosome operator=(const CChromosome& rhs);
	vector<float> *operator[](int nIndex);


//	friend ostream &operator<<(ostream &output, CChromosome &chrom) {
//		const int size = chrom.chromSize();
//		cout << "ERROR: "
//		for (unsigned short i = 0; i < size; i++) {
//			//cout << chrom[i][0] << "\n";
//			//output << chrom[i][0] << "\n";
//		}
//		output << "\n";
//		return output;
//	}
};

#endif /* CTRAIT_H_ */
