//
// File: DSO78F.h
// Created by: Julien Dutheil
// Created on: Tue Sep 02 12:49 2008
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

#ifndef _DSO78F_H_
#define _DSO78F_H_

#include "ProteinSubstitutionModelWithFrequencies.h"

// From SeqLib:
#include <Seq/ProteicAlphabet.h>

namespace bpp
{

/**
 * @brief The Dayhoff, Schwartz and Orcutt substitution model for proteins with free equilibrium frequencies
 *
 * Exchangeabilities have been computed using the DCMut method of Kosiol and Goldman.
 * The exchangability matrix is normalized so that \f$Q = S . \pi\f$ and \f$\sum_i Q_{i,i}\pi_i = -1\f$.
 * Eigen values and vectors are obtained numerically.
 * 
 * References:
 * - Dayhoff MO, Schwartz RM and Orcutt BC (1978), _A model of evolutionary change in proteins_, 5(3) 345-352, in _Atlas of Protein Sequence and Structure_. 
 * - Kosiol C and Goldman N (2005), _Molecular Biology And Evolution_ 22(2) 193-9. 
 */
class DSO78F:
  public ProteinSubstitutionModelWithFrequencies
{
	public:
		DSO78F(const ProteicAlphabet * alpha);

		virtual ~DSO78F() {}

#ifndef NO_VIRTUAL_COV
    DSO78F*
#else
    Clonable*
#endif
    clone() const { return new DSO78F(*this); }
    
	public:
		string getName() const { return "DSO78+F"; }

};

} //end of namespace bpp.

#endif	//_DSO78F_H_

