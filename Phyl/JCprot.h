//
// File: JCprot.h
// Created by:  <Julien.Dutheil@univ-montp2.fr>
// Created on: Tue May 25 16:04:36 2003
//

#ifndef _JCPROT_H_
#define _JCPROT_H_

#include <Seq/ProteicAlphabet.h>
#include "ProteinSubstitutionModel.h"

class JCprot : public ProteinSubstitutionModel
{
	public:
		JCprot(const ProteicAlphabet * alpha);
		~JCprot();
	
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


#endif	//_JCPROT_H_
