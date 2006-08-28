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
 * This optimizer is a combination of Newton-Raphson (NR) and Brent (B) algorithm.
 * Parameters are splitted in two kinds: some to be optimized with NS (type 1), the other with B (type 2).
 * 
 * One step of optimization consists of the following:
 * 1) Estimate type 1 parameters jointly with NR,
 * 2) Estimate separately each type 2 parameter with B.
 * 
 * A PseudoNewtonOptimizer object is used for type 1 parameters, a SimpleMultiDimensions
 * optimizer for type 2 parameters.
 *
 * For now, type 1 parameters are branch lengths, type 2 are other parameters (substitution model and rate distribution).
 * 
 * To decrease the optimization time, the precision of the two optimizers can be increased progressively:
 * if @f$\varepsilon@f$ is the final precision required, one may consider using a precision increment of @f$\sigma=\log_10(\varepsilon/n)@f$, where @f$n@f$ is the number of progressive steps.
 * During the first step optimization step, the precisions of type 1 and 2 optimizers are set to @f$10^{\sigma}@f$, @f$10^{2\sigma}@f$ for step 2, ... until precision @f$10^{n\sigma}=\varepsilon@f$ at step @f$n@f$ and later.
 * This saves some time spending in the first steps of the estimation.
 * The number of steps @f$n@f$ is set in the constructor of the optimizer.
 * 
 * @see PseudoNewtonOptimizer, SimpleMultiDimensions
 */
class NewtonBrentMetaOptimizer: public virtual AbstractOptimizer
{
	
	protected:
		ParameterList _newtonParameters;
		ParameterList _brentParameters;
		PseudoNewtonOptimizer * _newtonOptimizer;
		SimpleMultiDimensions * _brentOptimizer;
		unsigned int _nbNewtonParameters;
		unsigned int _nbBrentParameters;
    unsigned int _n;
    double _precisionStep;
    unsigned int _stepCount;
		
	public:
    /**
     * @brief Build a new NewtonBrentMetaOptimizer object.
     *
     * @param tl A likelihood function.
     * @param n  The number of progressive steps to use in optimization.
     */
		NewtonBrentMetaOptimizer(DiscreteRatesAcrossSitesTreeLikelihood * tl, unsigned int n = 1);
		virtual ~NewtonBrentMetaOptimizer();

	public:
		
		void init(const ParameterList & parameters) throw (Exception);
		double optimize() throw (Exception);
    double step() throw (Exception);
		void setFunction(Function * function)
		{
			AbstractOptimizer::setFunction(dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(function));
		}

};

#endif //_NEWTONBRENTMETAOPTIMIZER_H_

