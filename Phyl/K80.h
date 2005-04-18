//
// File: K80.h
// Created by:  <Julien.Dutheil@univ-montp2.fr>
// Created on: Tue May 27 15:24:30 2003
//

#ifndef _K2P_H_
#define _K2P_H_


#include "NucleotideSubstitutionModel.h"

// From SeqLib:
#include <Seq/NucleicAlphabet.h>

class K80 : public NucleotideSubstitutionModel
{
	public:
		K80(const NucleicAlphabet * alpha, double kappa = 1.);
		~K80();

		double Pij_t    (int i, int j, double d) const;
		double dPij_dt  (int i, int j, double d) const;
		double d2Pij_dt2(int i, int j, double d) const;
		Mat getPij_t    (double d) const;
		Mat getdPij_dt  (double d) const;
		Mat getd2Pij_dt2(double d) const;

		string getName() const;
	
	protected:
		void updateMatrices();

};


#endif	//_K80_H_
