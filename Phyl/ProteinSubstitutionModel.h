//
// File: ProteinSubstitutionModel.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Wed Jan 21 13:59:18 2004
//

#ifndef _PROTEINSUBSTITUTIONMODEL_H_
#define _PROTEINSUBSTITUTIONMODEL_H_

#include "AbstractSubstitutionModel.h"

// From SeqLib:
#include <Seq/Alphabet.h>

class ProteinSubstitutionModel : public AbstractSubstitutionModel
{
	public:
		ProteinSubstitutionModel(const Alphabet * alpha);

		virtual ~ProteinSubstitutionModel();
};


#endif	//_PROTEINSUBSTITUTIONMODEL_H_
