//
// File: JTT92.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Wed Jan 21 14:09:43 2004
//

#ifndef _JTT92_H_
#define _JTT92_H_

#include "ProteinSubstitutionModel.h"

// From SeqLib:
#include <Seq/ProteicAlphabet.h>

class JTT92: public ProteinSubstitutionModel
{
	public:
		JTT92(const ProteicAlphabet * alpha);
		virtual ~JTT92() {}
			
	public:
		string getName() const;

};


#endif	//_JTT92_H_
