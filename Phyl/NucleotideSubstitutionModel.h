//
// File: NucleicSubstitutionModel.h
// Created by:  <@bogdanof>
// Created on: Tue May 27 11:03:53 2003
//

#ifndef _NUCLEICSUBSTITUTIONMODEL_H_
#define _NUCLEICSUBSTITUTIONMODEL_H_

#include "AbstractSubstitutionModel.h"

class NucleotideSubstitutionModel : public AbstractSubstitutionModel
{
	public:
		NucleotideSubstitutionModel(const Alphabet * alpha);
		virtual ~NucleotideSubstitutionModel();
};


#endif	//_NUCLEICSUBSTITUTIONMODEL_H_
