//
// File: TN93.h
// Created by: Julien Dutheil
// Created on: Thu Jan 22 10:26:51 2004
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

#ifndef _TN93_H_
#define _TN93_H_

#include "NucleotideSubstitutionModel.h"

// From NumCalc:
#include <NumCalc/Constraints.h>

// From SeqLib:
#include <Seq/NucleicAlphabet.h>
#include <Seq/SequenceContainer.h>

/**
 * @brief The Tamura and Nei (1993) substitution model for nucleotides.
 *
 * This model has two rate of transitions and one rate of transversion.
 * It also allows distinct equilibrium frequencies between A, C, G and T.
 * This models hence includes six parameters, two transition / transversion
 * relative rates \f$\kappa_1\f$ and \f$\kappa_2\f$, and four frequencies \f$\pi_A, \pi_C, \pi_G, \pi_T\f$.
 * These four frequencies are not independent parameters, since they have the constraint to
 * sum to 1. Usually, these parameters are measured from the data and not optimized.
 * \f[
 * \begin{pmatrix}
 * \cdots & r & \kappa_1 r & r \\ 
 * r & \cdots & r & \kappa_2 r \\ 
 * \kappa_1 r & r & \cdots & r \\ 
 * r & \kappa_2 r & r & \cdots \\ 
 * \end{pmatrix}
 * \f]
 * \f[
 * \pi = \left(\pi_A, \pi_C, \pi_G, \pi_T\right)
 * \f]
 * Normalization: \f$r\f$ is set so that \f$\sum_i Q_{i,i}\pi_i = -1\f$:
 * \f[
 * S = \frac{1}{P}\begin{pmatrix}
 * \frac{-\pi_T-\kappa_1\pi_G-\pi_C}{\pi_A} & 1 & \kappa_1 & 1 \\ 
 * 1 & \frac{-\kappa_2\pi_T-\pi_G-\pi_A}{\pi_C} & 1 & \kappa_2 \\ 
 * \kappa_1 & 1 & \frac{-\pi_T-\pi_C-\kappa_1\pi_A}{\pi_G} & 1 \\ 
 * 1 & \kappa_2 & 1 & \frac{-\pi_G-\kappa_2\pi_C-\pi_A}{\pi_T} \\ 
 * \end{pmatrix}
 * \f]
 * with \f$P=2\left(\pi_A * \pi_C + \pi_C * \pi_G + \pi_A * \pi_T + \pi_G * \pi_T + kappa_2 * \pi_C * \pi_T + \kappa_1 * \pi_A * \pi_G\right)\f$.
 *
 * The normalized generator is obtained by taking the dot product of \f$S\f$ and \f$\pi\f$:
 * \f[
 * Q = S . \pi = \frac{1}{P}\begin{pmatrix}
 * -\pi_T-\kappa_1\pi_G-\pi_C & \pi_C & \kappa_1\pi_G & \pi_T \\ 
 * \pi_A & -\kappa_2\pi_T-\pi_G-\pi_A & \pi_G & \kappa_2\pi_T \\ 
 * \kappa_1\pi_A & \pi_C & -\pi_T-\pi_C-\kappa_1\pi_A & \pi_T \\ 
 * \pi_A & \kappa_2\pi_C & \pi_G & -\pi_G-\kappa_2\pi_C-\pi_A \\ 
 * \end{pmatrix}
 * \f]
 *
 * For now, the generator of this model is diagonalized numericaly.
 * See AbstractSubstitutionModel for details of how the porbabilities are computed.
 *
 * The parameters are named \c "kappa1", \c "kappa2", \c "piA", \c "piC", \c "piG" and \c "piT"
 * and their values may be retrieve with the command 
 * \code
 * getParameterValue("kappa1")
 * \endcode
 * for instance.
 *
 * Reference:
 * - Tamura N and Nei K (1993), _Molecular Biology And Evolution_ 10(3) 512-26. 
  */
class TN93 : public virtual NucleotideSubstitutionModel
{
	protected:
		Constraint * piConstraint;

	public:
		TN93(
			const NucleicAlphabet * alpha,
			double kappa1 = 1.,
			double kappa2 = 1.,
			double piA = 0.25,
			double piC = 0.25,
			double piG = 0.25,
			double piT = 0.25);
	
		virtual ~TN93();

		double Pij_t    (int i, int j, double d) const;
		double dPij_dt  (int i, int j, double d) const;
		double d2Pij_dt2(int i, int j, double d) const;
		RowMatrix<double> getPij_t    (double d) const;
		RowMatrix<double> getdPij_dt  (double d) const;
		RowMatrix<double> getd2Pij_dt2(double d) const;

		string getName() const;
	
		/**
		 * @brief This method is over-defined to actualize the corresponding parameters piA, piT, piG and piC too.
		 */
		void setFreqFromData(const SequenceContainer & data);

	protected:
		void updateMatrices();

};

#endif	//_TN93_H_

