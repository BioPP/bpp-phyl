//
// File: PseudoNewtonOptimizer.h
// Created by; Julien Dutheil
// Created on: Tue Nov 16 11:52 2004
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

#ifndef _PSEUDONEWTONOPTIMIZER_H_
#define _PSEUDONEWTONOPTIMIZER_H_

// From NumCalc:
#include <NumCalc/AbstractOptimizer.h>

namespace bpp
{

/**
 * @brief This Optimizer implements Newton's algorithm for finding a minimum of a function.
 * This is in fact a modified version of the algorithm, as suggested by Nicolas Galtier, for
 * the purpose of optimizing phylogenetic likelihoods.
 *
 * Only second simple order derivative are computed, no cross derivative, following Galtier's
 * algorithm.
 * Felsenstein and Churchill's (1996) correction is applied when new trial as a likelihood
 * lower than the starting point.
 */
class PseudoNewtonOptimizer:
  public AbstractOptimizer
{
  public:
    class PNStopCondition:
      public AbstractOptimizationStopCondition
    {
      public:
        PNStopCondition(PseudoNewtonOptimizer * pno):
          AbstractOptimizationStopCondition(pno) {}
        virtual ~PNStopCondition() {}

#ifndef NO_VIRTUAL_COV
        PNStopCondition*
#else
        Clonable*
#endif
        clone() const { return new PNStopCondition(*this); }
          
      public:
				void init() {}
        bool isToleranceReached() const;
    };
   
  friend class PNStopCondition;
         
	protected:

		ParameterList _previousPoint; // Current point is in _parameters

		double _previousValue;

		unsigned int _n; // Number of parameters

		vector<string> _params; // All parameter names

    string _mode;

    double _maxCorrection;

	public:

		PseudoNewtonOptimizer(DerivableSecondOrder * function);

		virtual ~PseudoNewtonOptimizer() {}

#ifndef NO_VIRTUAL_COV
    PseudoNewtonOptimizer*
#else
    Clonable*
#endif
    clone() const { return new PseudoNewtonOptimizer(*this); }

	public:
		/**
		 * @name The Optimizer interface.
		 *
		 * @{
		 */
		double getFunctionValue() const throw (NullPointerException) { return _currentValue; }
		/** @} */

		void doInit(const ParameterList & params) throw (Exception);

		double doStep() throw (Exception);

    void setMaximumNumberOfCorrections(unsigned int mx) { _maxCorrection = mx; }

};

} //end of namespace bpp.

#endif //_PSEUDONEWTONOPTIMIZER_H_

