//
// File: TN93.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Thu Jan 22 10:26:51 2004
//

#ifndef _TN93_H_
#define _TN93_H_

#include "NucleotideSubstitutionModel.h"

// From NumCalc:
#include <NumCalc/Constraints.h>

// From SeqLib:
#include <Seq/NucleicAlphabet.h>
#include <Seq/SequenceContainer.h>

class TN93 : public NucleotideSubstitutionModel
{
	protected:
		Constraint * piConstraint;
		void updateMatrices();

	public:
		TN93(
			const NucleicAlphabet * alpha,
			double kappa1 = 1.,
			double kappa2 = 1.,
			double piA = 0.25,
			double piC = 0.25,
			double piG = 0.25,
			double piT = 0.25);
	
		~TN93();

		double Pij_t    (int i, int j, double d) const;
		double dPij_dt  (int i, int j, double d) const;
		double d2Pij_dt2(int i, int j, double d) const;
		Mat getPij_t    (double d) const;
		Mat getdPij_dt  (double d) const;
		Mat getd2Pij_dt2(double d) const;

		string getName() const;
	
		/**
		 * @brief This method is over-defined to actualize the corresponding parameters piA, piT, piG and piC too.
		 */
		void setFreqFromData(const SequenceContainer & data);

};


#endif	//_TN93_H_
