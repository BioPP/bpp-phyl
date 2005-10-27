//
// File: GTR.h
// Created by: Julien Dutheil
// Created on: Tue Oct 25 10:17 2005
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

#ifndef _GTR_H_
#define _GTR_H_

#include "NucleotideSubstitutionModel.h"

// From NumCalc:
#include <NumCalc/Constraints.h>

// From SeqLib:
#include <Seq/NucleicAlphabet.h>

/**
 * @brief General Time-Reversible model for nucleotides.
 *
 * Parametrization:
 *
 * - Exchangeability matrix:
 *
 *   \f[
 *   	 \left(
 *     \begin{array}{cccc}
 *     --- & a   & b   & c   \\
 *     a   & --- & d   & e   \\
 *     b   & d   & --- & f   \\
 *     c   & e   & f   & --- \\
 *     \end{array}
 *     \right)
 *   \f]
 *
 * - Frequencies:
 *   \f[
 *     \left(
 *   	 \begin{array}{c}
 *   	 \pi_A \\
 *   	 \pi_C \\
 *   	 \pi_G \\
 *   	 \pi_T \\
 *   	 \end{array}
 *   	 \right)
 *   \f]
 */
class GTR : public virtual NucleotideSubstitutionModel
{
	protected:
		Constraint * piConstraint;
		double _a, _b, _c, _d, _e, _f;

	public:
		GTR(
			const NucleicAlphabet * alpha,
			double a = 1.,
			double b = 1.,
			double c = 1.,
			double d = 1.,
			double e = 1.,
			double f = 1.,
			double piA = 0.25,
			double piC = 0.25,
			double piG = 0.25,
			double piT = 0.25);
	
		virtual ~GTR();

		string getName() const;

		/**
		 * @brief This method is over-defined to actualize the corresponding parameters piA, piT, piG and piC too.
		 */
		void setFreqFromData(const SequenceContainer & data);
};

#endif	//_GTR_H_

