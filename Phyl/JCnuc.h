//
// File: JCnuc.h
// Created by:  <@bogdanof>
// Created on: Tue May 27 16:04:36 2003
//

#ifndef _JCNUC_H_
#define _JCNUC_H_

#include "NucleotideSubstitutionModel.h"

class JCnuc : public NucleotideSubstitutionModel
{
	public:
		JCnuc(const Alphabet * alpha);
		~JCnuc();
	
		double Pij_t    (int i, int j, double d) const;
		double dPij_dt  (int i, int j, double d) const;
		double d2Pij_dt2(int i, int j, double d) const;
		Matrix getPij_t    (double d) const;
		Matrix getdPij_dt  (double d) const;
		Matrix getd2Pij_dt2(double d) const;

		string getName() const;
	
	protected:
		void fillMatrices();
};


#endif	//_JCNUC_H_
