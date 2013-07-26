//
// File: MixturePhyloLikelihood.cpp
// Created by: Laurent Guéguen
// Created on: vendredi 12 juillet 2013, à 14h 55
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

#include "MixturePhyloLikelihood.h"

using namespace std;
using namespace bpp;
using namespace newlik;

/******************************************************************************/

MixturePhyloLikelihood::MixturePhyloLikelihood(
  SubstitutionProcessCollection* processColl,
  bool verbose,
  bool patterns) :
  MultiPhyloLikelihood(processColl, verbose, patterns),
  simplex_(processColl_->getNumberOfSubstitutionProcess(), 1)
{
  // initialize parameters:
  addParameters_(simplex_.getIndependentParameters());
}

MixturePhyloLikelihood::MixturePhyloLikelihood(
  const SiteContainer& data,
  SubstitutionProcessCollection* processColl,
  bool verbose,
  bool patterns) :
  MultiPhyloLikelihood(data, processColl, verbose, patterns),
  simplex_(processColl_->getNumberOfSubstitutionProcess(), 1)
{
  // initialize parameters:
  addParameters_(simplex_.getIndependentParameters());
  minusLogLik_ = -getLogLikelihood();
}

void MixturePhyloLikelihood::fireParameterChanged(const ParameterList& parameters)
{
  MultiPhyloLikelihood::fireParameterChanged(parameters);
  simplex_.matchParametersValues(parameters);

  minusLogLik_ = -getLogLikelihood();
}

ParameterList MixturePhyloLikelihood::getNonDerivableParameters() const
{
  ParameterList pl = processColl_->getNonDerivableParameters();
  pl.addParameters(simplex_.getParameters());
  return pl;
}

double MixturePhyloLikelihood::getLikelihoodForASite(size_t site) const
{
  double x = 0;
  for (size_t i = 0; i < vpTreelik_.size(); i++)
  {
    x += vpTreelik_[i]->getLikelihoodForASite(site) * simplex_.prob(i);
  }

  return x;
}

double MixturePhyloLikelihood::getDLikelihoodForASite(size_t site) const
{
  double x = 0;
  for (size_t i = 0; i < vpTreelik_.size(); i++)
  {
    x += vpTreelik_[i]->getDLikelihoodForASite(site) * simplex_.prob(i);
  }

  return x;
}

double MixturePhyloLikelihood::getD2LikelihoodForASite(size_t site) const
{
  double x = 0;
  for (size_t i = 0; i < vpTreelik_.size(); i++)
  {
    x += vpTreelik_[i]->getD2LikelihoodForASite(site) * simplex_.prob(i);
  }

  return x;
}


