//
// File: HmmPhyloLikelihood.cpp
// Created by: Laurent Guéguen
// Created on: lundi 23 septembre 2013, à 22h 56
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

#include "HmmPhyloLikelihood.h"

#include "HmmPhyloEmissionProbabilities.h"

#include <Bpp/Numeric/Hmm/LogsumHmmLikelihood.h>
#include <Bpp/Numeric/Hmm/FullHmmTransitionMatrix.h>

using namespace std;
using namespace bpp;
using namespace newlik;

/******************************************************************************/

HmmPhyloLikelihood::HmmPhyloLikelihood(
  SubstitutionProcessCollection* processColl,
  char recursivity,
  bool verbose,
  bool patterns) :
  MultiPhyloLikelihood(processColl, recursivity, verbose, patterns),
  Hmm_(0)
{
  HmmProcessAlphabet* hpa=new HmmProcessAlphabet(processColl_.get());

  HmmTransitionMatrix* hptm=new FullHmmTransitionMatrix(hpa);

  HmmPhyloEmissionProbabilities* hpep=new HmmPhyloEmissionProbabilities(hpa, this);
    
  Hmm_ = auto_ptr<LogsumHmmLikelihood>(new LogsumHmmLikelihood(hpa, hptm, hpep, ""));
    
  // initialize parameters:
  addParameters_(Hmm_->getHmmTransitionMatrix().getParameters());
  addParameters_(Hmm_->getHmmStateAlphabet().getParameters());
}

HmmPhyloLikelihood::HmmPhyloLikelihood(
  const SiteContainer& data,
  SubstitutionProcessCollection* processColl,
  char recursivity,
  bool verbose,
  bool patterns) :
  MultiPhyloLikelihood(data, processColl, recursivity, verbose, patterns),
  Hmm_(0)
{
  HmmProcessAlphabet* hpa=new HmmProcessAlphabet(processColl_.get());

  HmmTransitionMatrix* hptm=new FullHmmTransitionMatrix(hpa);

  HmmPhyloEmissionProbabilities* hpep=new HmmPhyloEmissionProbabilities(hpa, this);

  Hmm_ = auto_ptr<LogsumHmmLikelihood>(new LogsumHmmLikelihood(hpa, hptm, hpep, ""));
  // initialize parameters:

  addParameters_(Hmm_->getHmmTransitionMatrix().getParameters());
  addParameters_(Hmm_->getHmmStateAlphabet().getParameters());

  minusLogLik_ = -getLogLikelihood();
}

void HmmPhyloLikelihood::fireParameterChanged(const ParameterList& parameters)
{
  MultiPhyloLikelihood::fireParameterChanged(parameters);

  Hmm_->matchParametersValues(parameters);
  
  minusLogLik_ = -Hmm_->getLogLikelihood();
}

ParameterList HmmPhyloLikelihood::getNonDerivableParameters() const
{
  ParameterList pl = processColl_->getNonDerivableParameters();
  pl.addParameters(Hmm_->getHmmTransitionMatrix().getParameters());
  pl.addParameters(Hmm_->getHmmStateAlphabet().getParameters());

  return pl;
}

