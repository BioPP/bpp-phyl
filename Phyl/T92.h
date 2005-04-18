//
// File: T92.h
// Created by:  <Julien.Dutheil@univ-montp2.fr>
// Created on: Mon May 26 14:41:24 2003
//

#ifndef _T92_H_
#define _T92_H_

#include "NucleotideSubstitutionModel.h"

// From NumCalc:
#include <NumCalc/Constraints.h>

// From SeqLib:
#include <Seq/NucleicAlphabet.h>
#include <Seq/SequenceContainer.h>

class T92 : public NucleotideSubstitutionModel
{
	protected:
		Constraint * thetaConstraint;
		void updateMatrices();

	public:
		T92(const NucleicAlphabet * alpha, double kappa = 1., double theta = 0.5);
		~T92();

		double Pij_t    (int i, int j, double d) const;
		double dPij_dt  (int i, int j, double d) const;
		double d2Pij_dt2(int i, int j, double d) const;
		Mat getPij_t(double d) const;
		Mat getdPij_dt(double d) const;
		Mat getd2Pij_dt2(double d) const;

		string getName() const;
	
		/**
		 * @brief This method is over-defined to actualize the 'theta' parameter too.
		 */
		void setFreqFromData(const SequenceContainer & data);
};

#endif	//_T92_H_
