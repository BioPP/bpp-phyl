//
// File: T92.h
// Created by:  <@bogdanof>
// Created on: Mon May 26 14:41:24 2003
//

#ifndef _T92_H_
#define _T92_H_

#include "NucleotideSubstitutionModel.h"

// From NumCalc:
#include <NumCalc/Constraints.h>

// From SeqLib:
#include <Seq/SequenceContainer.h>

class T92 : public NucleotideSubstitutionModel
{
	protected:
		Constraint * thetaConstraint;
		void fillMatrices();

	public:
		T92(const Alphabet * alpha, double kappa = 1., double theta = 0.5);
		~T92();

		double Pij_t    (int i, int j, double d) const;
		double dPij_dt  (int i, int j, double d) const;
		double d2Pij_dt2(int i, int j, double d) const;
		Matrix getPij_t(double d) const;
		Matrix getdPij_dt(double d) const;
		Matrix getd2Pij_dt2(double d) const;

		string getName() const;
	
	public:
		//specific method:
		void setThetaFromData(const SequenceContainer & data);
};

#endif	//_T92_H_
