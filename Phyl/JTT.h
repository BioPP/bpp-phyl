//
// File: JTT.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Wed Jan 21 14:09:43 2004
//

#ifndef _JTT_H_
#define _JTT_H_

#include "ProteinSubstitutionModel.h"

// From SeqLib:
#include <Seq/Alphabet.h>

class JTT : public ProteinSubstitutionModel
{
	public:
		JTT(const Alphabet * alpha);
		virtual ~JTT() {}
			
	protected:
		void fillMatrices();
	
};


#endif	//_JTT_H_
