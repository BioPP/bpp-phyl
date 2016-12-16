//
// File: HmmOfAlignedPhyloLikelihood.cpp
// Created by: Laurent Guéguen
// Created on: mercredi 28 octobre 2015, à 23h 39
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

#include "HmmOfAlignedPhyloLikelihood.h"

#include <Bpp/Numeric/Hmm/LogsumHmmLikelihood.h>
#include <Bpp/Numeric/Hmm/FullHmmTransitionMatrix.h>

using namespace std;
using namespace bpp;

/******************************************************************************/

HmmOfAlignedPhyloLikelihood::HmmOfAlignedPhyloLikelihood(
  PhyloLikelihoodContainer* pC, const std::vector<size_t>& nPhylo) :
  AbstractPhyloLikelihood(),
  AbstractAlignedPhyloLikelihood(0),
  SetOfAlignedPhyloLikelihood(pC, nPhylo),
  hma_(),
  htm_(),
  hpep_(),
  hmm_()
{
  hma_ = unique_ptr<HmmPhyloAlphabet>(new HmmPhyloAlphabet(*this));

  htm_ = unique_ptr<FullHmmTransitionMatrix>(new FullHmmTransitionMatrix(hma_.get(), "HMM."));

  hpep_ = unique_ptr<HmmPhyloEmissionProbabilities>(new HmmPhyloEmissionProbabilities(hma_.get()));

  hmm_ = unique_ptr<LogsumHmmLikelihood>(new LogsumHmmLikelihood(hma_.get(), htm_.get(), hpep_.get(), false));

  addParameters_(htm_->getParameters());
}

void HmmOfAlignedPhyloLikelihood::setNamespace(const std::string& nameSpace)
{
  SetOfAlignedPhyloLikelihood::setNamespace(nameSpace);
  hmm_->setNamespace(nameSpace);
}

void HmmOfAlignedPhyloLikelihood::fireParameterChanged(const ParameterList& parameters)
{
  SetOfAlignedPhyloLikelihood::fireParameterChanged(parameters);
  hpep_->update();
  hmm_->matchParametersValues(parameters);
}

ParameterList HmmOfAlignedPhyloLikelihood::getNonDerivableParameters() const
{
  ParameterList pl=HmmOfAlignedPhyloLikelihood::getNonDerivableParameters();
  
  pl.includeParameters(htm_->getParameters());

  return pl;
}
