//
// File: JCprotF.h
// Created by: Julien Dutheil
// Created on: Tue Sept 02 11:24 2008
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004)

This software is a computer program whose purpose is to provide classes
for phylogenetic data analysis.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef _JCPROTF_H_
#define _JCPROTF_H_

#include <Seq/ProteicAlphabet.h>
#include "ProteinSubstitutionModelWithFrequencies.h"

namespace bpp
{

/**
 * @brief The Jukes-Cantor substitution model for proteins with free equilibrium frequencies.
 *
 * This model is the same as JCprot, but allows the equilibrium frequencies to be estimated directly from the data.
 *
 * @see JCprot, ProteinSubstitutionModelWithFrequencies
 */
class JCprotF:
  public ProteinSubstitutionModelWithFrequencies
{
  protected:
    mutable double _exp;
    mutable RowMatrix<double> _p;

	public:
		JCprotF(const ProteicAlphabet * alpha);

		virtual ~JCprotF() {}
	
#ifndef NO_VIRTUAL_COV
    JCprotF*
#else
    Clonable*
#endif
    clone() const { return new JCprotF(*this); }

  public:
		string getName() const { return "JCprot+F"; }
	
};

} //end of namespace bpp.

#endif	//_JCPROTF_H_

