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
 * @brief The General Time-Reversible substitution model for nucleotides.
 *
 * This model is sometimes called the REV model, following Yang (1994), see references.
 * It was used in Lanave et al (1984), described in Tavare et al. 1886 and Rodriguez et al. 1990.
 * It is the most general reversible one, it has 6 substitution rates and 4 frequency
 * parameters.
 * We used the parametrization proposed by Yang (1994):
 * \f[
 * S = \begin{pmatrix}
 * \cdots & d & f & b \\ 
 * d & \cdots & e & a \\ 
 * f & e & \cdots & c \\ 
 * b & a & c & \cdots \\ 
 * \end{pmatrix}
 * \f]
 * \f[
 * \pi = \left(\pi_A, \pi_C, \pi_G, \pi_T\right)
 * \f]
 * Normalization: we set \f$f\f$ to 1, and scale the matrix so that \f$\sum_i Q_{i,i}\pi_i = -1\f$.
 * Parameters \f$a,b,c,d,e\f$ are hence relative rates.
 * \f[
 * S = \frac{1}{P}\begin{pmatrix}
 * \frac{-b\pi_T-\pi_G-d\pi_C}{\pi_A} & d & 1 & b \\ 
 * d & \frac{-a\pi_T-e\pi_G-d\pi_A}{\pi_C} & e & a \\ 
 * 1 & e & \frac{-c\pi_T-e\pi_C-\pi_A}{\pi_G} & c \\ 
 * b & a & c & \frac{-c\pi_G-a\pi_C-b\pi_A}{\pi_T} \\ 
 * \end{pmatrix}
 * \f]
 * with \f{eqnarray*}
 * P &=& \pi_A(d\pi_C+ \pi_G+b\pi_T)\\
 *	 &+& \pi_C(d\pi_A+e\pi_G+a\pi_T)\\
 *   &+& \pi_G( \pi_A+e\pi_C+c\pi_T)\\
 *	 &+& \pi_T(b\pi_A+a\pi_C+c\pi_G)\\
 *   &=& 2*(a*\pi_C*\pi_T+b*\pi_A*\pi_T+c*\pi_G*\pi_T+d*\pi_A*\pi_C+e*\pi_C*\pi_G+\pi_A*\pi_G)
 * \f}
 *
 * The normalized generator is obtained by taking the dot product of \f$S\f$ and \f$\pi\f$:
 * \f[
 * Q = S . \pi = \frac{1}{P}\begin{pmatrix}
 * -b\pi_T-\pi_G-d\pi_C & d\pi_C & \pi_G & b\pi_T \\ 
 * d\pi_A & -a\pi_T-e\pi_G-d\pi_A & e\pi_G & a\pi_T \\ 
 * \pi_A & e\pi_C & -c\pi_T-e\pi_C-\pi_A & c\pi_T \\ 
 * b\pi_A & a\pi_C & c\pi_G & -c\pi_G-a\pi_C-b\pi_A \\ 
 * \end{pmatrix}
 * \f]
 *
 * For now, the generator of this model is diagonalized numericaly.
 * See AbstractSubstitutionModel for details of how the porbabilities are computed.
 *
 * The parameters are named \c "a", \c "b", \c "c", \c "d", \c "e", \c "piA", \c "piC",
 * \c "piG" and \c "piT" and their values may be retrieve with the command 
 * \code
 * getParameterValue("a")
 * \endcode
 * for instance.
 * 
 * Reference:
 * - Yang Z (1994), _Journal Of Molecular Evolution_ 39(1) 105-11.
 * - Lanave C, Preparata G, Saccone C and Serio G (1984), _Journal Of Molecular Evolution_ 20 86-93.
 * - Tavaré S (1986), _Lect. Math. Life Sci._ 17 57-86.
 * - Rodriguez F (1990, _Journal Of Theoretical Biology_ 142(4) 485-501.
 */
class GTR:
  public NucleotideSubstitutionModel
{
	protected:
    double _a, _b, _c, _d, _e, _piA, _piC, _piG, _piT, _p;

	public:
		GTR(
			const NucleicAlphabet * alpha,
			double a = 1.,
			double b = 1.,
			double c = 1.,
			double d = 1.,
			double e = 1.,
			double piA = 0.25,
			double piC = 0.25,
			double piG = 0.25,
			double piT = 0.25);
	
		virtual ~GTR() {}

    GTR * clone() const { return new GTR(*this); }

  public:
		string getName() const;

    void updateMatrices();

		/**
		 * @brief This method is redefined to actualize the corresponding parameters piA, piT, piG and piC too.
		 */
		void setFreqFromData(const SequenceContainer & data);
};

#endif	//_GTR_H_

