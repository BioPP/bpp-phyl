//
// File: L95.h
// Created by: Julien Dutheil
// Created on: Tue Nov 4 11:46 2008
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

#ifndef _L95_H_
#define _L95_H_

#include "NucleotideSubstitutionModel.h"

// From NumCalc:
#include <NumCalc/Constraints.h>

// From SeqLib:
#include <Seq/NucleicAlphabet.h>

namespace bpp
{

/**
 * @brief The strand Symmetric substitution model for nucleotides.
 *
 * We used a parametrization derived from Hobolth et al 2007
 * This model contains 6 parameters:
 * \f[
 * S = \begin{pmatrix}
 * \cdots & \beta & 1 & \gamma \\ 
 * \beta & \cdots & \delta & 1 \\ 
 * 1 & \delta & \cdots & \beta \\ 
 * \gamma & 1 & \beta & \cdots \\ 
 * \end{pmatrix}
 * \f]
 * The quilibrium frequencies 
 * \f[
 * \pi = \left(1-\frac{\theta}{2}, \frac{\theta}{2}, \frac{\theta}{2}, 1-\frac{\theta}{2}\right)
 * \f]
 * This models hence includes four parameters, three relative rates \f$\beta, \gamma, \delta\f$ and the GC content \f$\theta\f$.
 *
 * Normalization: we set \f$f\f$ to 1, and scale the matrix so that \f$\sum_i Q_{i,i}\pi_i = -1\f$.
 * The normalized generator is obtained by taking the dot product of \f$S\f$ and \f$\pi\f$:
 * \f[
 * Q = S . \pi = \frac{1}{P}\begin{pmatrix}
 * -\gamma\pi_T-\pi_G-\beta\pi_C & \beta\pi_C & \pi_G & \gamma\pi_T \\ 
 * \beta\pi_A & -\pi_T-\delta\pi_G-\beta\pi_A & \delta\pi_G & \pi_T \\ 
 * \pi_A & \delta\pi_C & -\beta\pi_T-\delta\pi_C-\pi_A & \beta\pi_T \\ 
 * \gamma\pi_A & \pi_C & \beta\pi_G & -\beta\pi_G-\pi_C-\gamma\pi_A \\ 
 * \end{pmatrix}
 * \f]
 * where P is the normnalisation constant.
 * For now, the generator of this model is diagonalized numericaly.
 * See AbstractSubstitutionModel for details of how the probabilities are computed.
 *
 * The parameters are named \c "beta", \c "gamma", \c "delta", and \c "theta"
 * and their values may be retrieve with the command 
 * \code
 * getParameterValue("beta")
 * \endcode
 * for instance.
 * 
 * Reference:
 * - Hobolth A, Christensen O Fm Mailund T, Schierup M H (2007), _PLoS Genetics_ 3(2) e7.
 * - Lobry J R (1995), _Journal Of Molecular Evolution_ 40 326-330.
 */
class L95:
  public NucleotideSubstitutionModel
{
	protected:
    double _beta, _gamma, _delta, _theta, _piA, _piC, _piG, _piT;

	public:
		L95(
			const NucleicAlphabet * alpha,
			double beta = 1.,
			double gamma = 1.,
			double delta = 1.,
			double theta = 0.5);
	
		virtual ~L95() {}

#ifndef NO_VIRTUAL_COV
    L95*
#else
    Clonable*
#endif
    clone() const { return new L95(*this); }

  public:
		string getName() const { return "Lobry 1995"; }

    void updateMatrices();

		/**
		 * @brief This method is redefined to actualize the corresponding parameters theta too.
		 */
		void setFreqFromData(const SequenceContainer & data);
};

} //end of namespace bpp.

#endif	//_L95_H_

