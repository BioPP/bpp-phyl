//
// File: G2001.h
// Created by: Julien Dutheil
// Created on: Mon Aug 07 18:31 2006
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

#ifndef _G2001_H_
#define _G2001_H_

#include "MarkovModulatedSubstitutionModel.h"

//From NumCalc:
#include <NumCalc/DiscreteDistribution.h>
#include <NumCalc/Parameter.h>

namespace bpp
{

/**
 * @brief Galtier's 2001 covarion model.
 *
 * This model is a subclass of the so-called Markov-modulated substitution models,
 * with a Jukes-Cantor rate matrix, of parameter @f$\nu@f$.
 * the original version uses a discrete @f$\Gamma@f$ distribution for rates, but you can
 * use it with virtually any rate distribution.
 *
 * @see MarkovModulatedSubstitutionModel
 *
 * Galtier N., Maximum-likelihood phylogenetic analysis under a covarion-like model (2001).
 * _Molecular Biology and Evolution_, 18:866-73.
 */
class G2001:
  public MarkovModulatedSubstitutionModel
{
  protected:
    DiscreteDistribution * _rDist;

  public:
    /**
     * @brief Build a new G2001 substitution model.
     *
     * @param model The substitution model to use. May be of any alphabet type.
     * @param rDist The discrete distribution for rates. The class will own the DiscreteDistribution object,
     * which will be deleted together with this instance.
     * @param nu    The rate matrix parameter.
     * @param normalizeRateChanges Tell if the rate transition matrix should be normalized.
     */
    G2001(ReversibleSubstitutionModel * model, DiscreteDistribution * rDist, double nu = 1., bool normalizeRateChanges = false):
      MarkovModulatedSubstitutionModel(model, normalizeRateChanges, "G01."), _rDist(rDist)
    {
      _nbRates = _rDist->getNumberOfCategories();
      _ratesExchangeability.resize(_nbRates, _nbRates);
      _rates.resize(_nbRates, _nbRates);
      _ratesFreq = vector<double>(_nbRates, 1./(double)_nbRates);
      _rDist->setNamespace("G01.rdist");
      addParameters_(_rDist->getIndependentParameters());
      Parameter p("G01.nu", nu, &Parameter::R_PLUS);
      addParameter_(p);
      updateRatesModel();
      updateMatrices();
    }

    G2001(const G2001 & model): MarkovModulatedSubstitutionModel(model)
    {
      _rDist = dynamic_cast<DiscreteDistribution *>(model._rDist->clone());
    }

    G2001 & operator=(const G2001 & model)
    {
      MarkovModulatedSubstitutionModel::operator=(model);
      _rDist = dynamic_cast<DiscreteDistribution *>(model._rDist->clone());
      return *this;
    }

    virtual ~G2001() { delete _rDist; }

#ifndef NO_VIRTUAL_COV
    G2001*
#else
    Clonable*
#endif
    clone() const { return new G2001(*this); }

  public:
    string getName() const { return _model->getName() + "+G2001"; }

    /**
     * @brief Re-definition of the super-class method to update the rate distribution too.
     *
     * @param parameters The parameters that have been modified.
     */
		void fireParameterChanged(const ParameterList & parameters)
    {
      _rDist->matchParametersValues(parameters);
      MarkovModulatedSubstitutionModel::fireParameterChanged(parameters);
    }
    
  protected:
    void updateRatesModel()
    {
      double nu = getParameterValue("nu");
      for(unsigned int i = 0; i < _nbRates; i++)
      {
         _rates(i,i) = _rDist->getCategory(i);
         for(unsigned int j = 0; j < _nbRates; j++)
         {
            if(i==j)
            {
              _ratesExchangeability(i,j) = -(double)_nbRates*nu;
            }
            else
            {
              _ratesExchangeability(i,j) = _nbRates*nu/(_nbRates-1);
            }
         }
      }
    }
	
};

} //end of namespace bpp.

#endif // _G2001_H_

