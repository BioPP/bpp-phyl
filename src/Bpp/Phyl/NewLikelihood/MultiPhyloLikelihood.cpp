//
// File: MultiPhyloLikelihood.cpp
// Created by: Laurent Guéguen
// Created on: vendredi 12 juillet 2013, à 00h 32
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

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

#include "MultiPhyloLikelihood.h"

#include "SingleRecursiveTreeLikelihoodCalculation.h"

using namespace std;
using namespace bpp;
using namespace newlik;

/******************************************************************************/

MultiPhyloLikelihood::MultiPhyloLikelihood(
  SubstitutionProcessCollection* processColl,
  bool verbose,
  bool patterns) :
  AbstractParametrizable(""),
  data_(0),
  processColl_(processColl),
  recursivity_('S'),
  computeFirstOrderDerivatives_(true),
  computeSecondOrderDerivatives_(true),
  initialized_(false),
  verbose_(verbose),
  nbSites_(0),
  nbStates_(0),
  minusLogLik_(0),
  vpTreelik_()
{
  // initialize parameters:
  
  addParameters_(processColl_->getIndependentParameters());

  if (recursivity_ != 'S')
    throw Exception("MultiPhyloLikelihood::MultiPhyloLikelihood : Non simple recursivity not implemented yet!");

  size_t nbP = processColl_->getNumberOfSubstitutionProcess();

  for (size_t i = 0; i < nbP; i++)
  {
    SingleRecursiveTreeLikelihoodCalculation* rt = new SingleRecursiveTreeLikelihoodCalculation(
        processColl_->getSubstitutionProcess(i), patterns, i == 0);
    vpTreelik_.push_back(rt);
  }
}

/******************************************************************************/

MultiPhyloLikelihood::MultiPhyloLikelihood(
  const SiteContainer& data,
  SubstitutionProcessCollection* processColl,
  bool verbose,
  bool patterns) :
  AbstractParametrizable(""),
  data_(&data),
  processColl_(processColl),
  recursivity_('S'),
  computeFirstOrderDerivatives_(true),
  computeSecondOrderDerivatives_(true),
  initialized_(false),
  verbose_(verbose),
  nbSites_(0),
  nbStates_(0),
  minusLogLik_(0),
  vpTreelik_()
{
  // initialize parameters:
  addParameters_(processColl_->getIndependentParameters());

  if (recursivity_ != 'S')
    throw Exception("MultiPhyloLikelihood::MultiPhyloLikelihood : Non simple recursivity not implemented yet!");

  size_t nbP = processColl_->getNumberOfSubstitutionProcess();

  for (size_t i = 0; i < nbP; i++)
  {
    SingleRecursiveTreeLikelihoodCalculation* rt = new SingleRecursiveTreeLikelihoodCalculation(*data_->clone(),
        processColl_->getSubstitutionProcess(i), patterns, i == 0);
    vpTreelik_.push_back(rt);
  }

  setData(data);
}

/******************************************************************************/

void MultiPhyloLikelihood::setData(const SiteContainer& sites)
{
  for (size_t i = 0; i < vpTreelik_.size(); i++)
  {
    vpTreelik_[i]->setData(sites);
  }

  nbSites_ = sites.getNumberOfSites();
  nbStates_ = sites.getAlphabet()->getSize();

  initialized_ = true;
}

/******************************************************************************/

void MultiPhyloLikelihood::fireParameterChanged(const ParameterList& parameters)
{
  processColl_->matchParametersValues(parameters);

  if (parameters.size() > 0)
    for (size_t i = 0; i < vpTreelik_.size(); i++)
    {
      vpTreelik_[i]->computeTreeLikelihood();
    }
}

/******************************************************************************/

ParameterList MultiPhyloLikelihood::getSubstitutionProcessParameters() const
{
  return processColl_->getSubstitutionProcessParameters();
}

ParameterList MultiPhyloLikelihood::getSubstitutionModelParameters() const
{
  return processColl_->getSubstitutionModelParameters();
}

ParameterList MultiPhyloLikelihood::getRateDistributionParameters() const
{
  return processColl_->getRateDistributionParameters();
}

ParameterList MultiPhyloLikelihood::getRootFrequenciesParameters() const
{
  return processColl_->getRootFrequenciesParameters();
}

ParameterList MultiPhyloLikelihood::getBranchLengthsParameters() const
{
  return processColl_->getBranchLengthsParameters();
}

std::vector<const TreeTemplate<Node>* > MultiPhyloLikelihood::getTrees() const
{
  std::vector<const TreeTemplate<Node>* > vT;
  size_t i=1;
  while (processColl_->hasTreeNumber(i)){
    vT.push_back(&(processColl_->getTree(i)->getTree()));
    i++;
  }
  return vT;
}

/******************************************************************************/

Vdouble MultiPhyloLikelihood::getLikelihoodForEachSite() const
{
  Vdouble l(getNumberOfSites());
  for (unsigned int i = 0; i < l.size(); ++i)
  {
    l[i] = getLikelihoodForASite(i);
  }
  return l;
}

/******************************************************************************/

VVdouble MultiPhyloLikelihood::getLikelihoodForEachSiteForEachProcess() const
{
  VVdouble l(getNumberOfSites());
  for (size_t i = 0; i < l.size(); ++i)
    {
      Vdouble* l_i = &l[i];
      l_i->resize(getNumberOfSubstitutionProcess());
      for (size_t c = 0; c < l_i->size(); ++c)
        {
          (*l_i)[c] = getLikelihoodForASiteForAProcess(i, c);
        }
    }
  return l;
}

/******************************************************************************/

double MultiPhyloLikelihood::getDLogLikelihoodForASite(size_t site) const
{
  return getDLikelihoodForASite(site) / getLikelihoodForASite(site);
}

/******************************************************************************/

double MultiPhyloLikelihood::getD2LogLikelihoodForASite(size_t site) const
{
  return getD2LikelihoodForASite(site) / getLikelihoodForASite(site)
         - pow( getDLikelihoodForASite(site) / getLikelihoodForASite(site), 2);
}

/******************************************************************************/


void MultiPhyloLikelihood::setParameters(const ParameterList& parameters)
throw (ParameterNotFoundException, ConstraintException)
{
  setParametersValues(parameters);
}

/******************************************************************************/

double MultiPhyloLikelihood::getValue() const
throw (Exception)
{
  if (!isInitialized())
    throw Exception("MultiPhyloLikelihood::getValue(). Instance is not initialized.");
  return minusLogLik_;
}

/******************************************************************************/

double MultiPhyloLikelihood::getFirstOrderDerivative(const string& variable) const
throw (Exception)
{
  if (!hasParameter(variable))
    throw ParameterNotFoundException("MultiPhyloLikelihood::getFirstOrderDerivative().", variable);
  if (!processColl_->hasBranchLengthsParameter(variable))
  {
    throw Exception("Derivatives are only implemented for branch length parameters.");
  }
  computeDLikelihood_(variable);
  return -getDLogLikelihood();
}

/******************************************************************************/

double MultiPhyloLikelihood::getSecondOrderDerivative(const string& variable) const
throw (Exception)
{
  if (!hasParameter(variable))
    throw ParameterNotFoundException("MultiPhyloLikelihood::getSecondOrderDerivative().", variable);
  if (!processColl_->hasBranchLengthsParameter(variable))
  {
    throw Exception("Derivatives are only implemented for branch length parameters.");
  }
  computeD2Likelihood_(variable);
  return -getD2LogLikelihood();
}

