//
// File: NewtonBrentMetaOptimizer.h
// Created by; Julien Dutheil
// Created on: ue Nov 17 17:22 2004
//

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

Julien.Dutheil@univ-montp2.fr
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

#ifndef _NEWTONBRENTMETAOPTIMIZER_H_
#define _NEWTONBRENTMETAOPTIMIZER_H_

#include "AbstractHomogeneousTreeLikelihood.h"
#include "PseudoNewtonOptimizer.h"

// From NumCalc:
#include <NumCalc/SimpleMultiDimensions.h>
#include <NumCalc/AbstractOptimizer.h>

// From the STL:
#include <vector>
using namespace std;

/**
 * @brief Phylogenetic optimizer.
 *
 * This optimizer uses a pseudo-newton algorithm for estimating branch lengths,
 * and Brent's one-dimensional estimation algorithm for rate distribution and 
 * substitution model parameters.
 *
 * One step of optimization concists of the following:
 * 1) Estimate branch lengths (pseudo-newton),
 * 2) Estimate rate distribution parameters (loop over each parameter using Brent's method),
 * 3) Estimate substitution model parameters (loop over each parameter using Brent's method).
 * 
 * A PseudoNewtonOptimizer object is used for branch lengths estimations, a SimpleMultiDimensions
 * optimizer for rate distribution, and another for substitution model parameters.
 *
 * Furthermore, it is possible to roughly estimate all parameters before starting to loop
 * over all parameters (see setRoughEstimationEnabled).
 * 
 * @see PseudoNewtonOptimizer, SimpleMultiDimensions
 */
class NewtonBrentMetaOptimizer: public virtual AbstractOptimizer
{
	public:
		static double BRANCH_LENGTHS_TOL;
		static double RATE_DISTRIBUTION_TOL;
		static double SUBSTITUTION_MODEL_TOL;
	
	protected:
		ParameterList _rateDistributionParameters;
		ParameterList _substitutionModelParameters;
		ParameterList _branchLengthsParameters;
		SimpleMultiDimensions * _rateDistributionOptimizer;
		SimpleMultiDimensions * _substitutionModelOptimizer;
		PseudoNewtonOptimizer * _branchLengthsOptimizer;
		unsigned int _nbRateDistParams;
		unsigned int _nbSubsModParams;
		unsigned int _nbBranchLengths;
		bool _rough;
		
	public:
		NewtonBrentMetaOptimizer(DiscreteRatesAcrossSitesTreeLikelihood * tl);
		virtual ~NewtonBrentMetaOptimizer();

	public:
		
		void init(const ParameterList & parameters) throw (Exception);
		double optimize() throw (Exception);
		double step() throw (Exception) { return 0.; }
//		double getFunctionValue() const;
		void setFunction(Function * function)
		{
			AbstractOptimizer::setFunction(dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(function));
		}

    /**
     * @brief Enable rough estimation prior to final estimation.
     *
     * This generally speeds up the global estimation.
     * Branch lengths estimation is generally time-consuming, and
     * estimates highly depends on the shape of the rate distribution.
     * Roughly estimating the rate distribution parameters hence
     * save some time during the branch length estimation.
     *
     * The tolerance of rough estimation is fixed by the
     * BRANCH_LENGTHS_TOL, RATE_DISTRIBUTION_TOL and SUBSTITUTION_MODEL_TOL
     * static variables.
     *
     * @param yn A boolean.
     */
		void setRoughEstimationEnabled(bool yn) { _rough = yn; }
    /**
     * @return True If rough estimation is enabled.
     */
		bool isRoughEstimationEnabled() const { return _rough; }
};

#endif //_NEWTONBRENTMETAOPTIMIZER_H_

