//
// File: DSO78.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Tue Oct 05 18:49:44 2004
//

#ifndef _DSO78_H_
#define _DSO78_H_

#include "ProteinSubstitutionModel.h"

// From SeqLib:
#include <Seq/ProteicAlphabet.h>

class DSO78 : public ProteinSubstitutionModel
{
	public:
		DSO78(const ProteicAlphabet * alpha);
		virtual ~DSO78() {}
		
	public:
		string getName() const;

};


#endif	//_DSO78_H_
