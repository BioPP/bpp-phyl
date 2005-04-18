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
		JCnuc(const NucleicAlphabet * alpha);
		~JCnuc();
	
		double Pij_t    (int i, int j, double d) const;
		double dPij_dt  (int i, int j, double d) const;
		double d2Pij_dt2(int i, int j, double d) const;
		Mat getPij_t    (double d) const;
		Mat getdPij_dt  (double d) const;
		Mat getd2Pij_dt2(double d) const;

		string getName() const;
	
	protected:
		/**
		 * In the case of the model of Jukes & Cantor, this method is useless since
		 * the generator is fixed! No matrice can be changed... This method is only
		 * used in the constructor of the class.
		 */
		void updateMatrices();
};


#endif	//_JCNUC_H_
