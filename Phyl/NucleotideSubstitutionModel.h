//
// File: NucleicSubstitutionModel.h
// Created by:  <Julien.Dutheil@univ-montp2.fr>
// Created on: Tue May 27 11:03:53 2003
//

#ifndef _NUCLEICSUBSTITUTIONMODEL_H_
#define _NUCLEICSUBSTITUTIONMODEL_H_

#include "AbstractSubstitutionModel.h"

// From SeqLib:
#include <Seq/NucleicAlphabet.h>

class NucleotideSubstitutionModel : public AbstractSubstitutionModel
{
	public:
		NucleotideSubstitutionModel(const NucleicAlphabet * alpha);
		virtual ~NucleotideSubstitutionModel();
};


#endif	//_NUCLEICSUBSTITUTIONMODEL_H_
