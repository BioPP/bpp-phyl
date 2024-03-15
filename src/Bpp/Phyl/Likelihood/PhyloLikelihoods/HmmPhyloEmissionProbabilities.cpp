// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "HmmPhyloEmissionProbabilities.h"

using namespace bpp;
using namespace std;

HmmPhyloEmissionProbabilities::HmmPhyloEmissionProbabilities(std::shared_ptr<HmmPhyloAlphabet> alphabet) :
  AbstractParametrizable(""),
  context_(alphabet->getContext()),
  phylAlph_(alphabet),
  emProb_(),
  nbSites_(alphabet->getNumberOfSites())
{
  setHmmStateAlphabet(alphabet);
}


void HmmPhyloEmissionProbabilities::setHmmStateAlphabet(std::shared_ptr<HmmStateAlphabet> stateAlphabet)
{
  if (!stateAlphabet)
    throw HmmUnvalidAlphabetException("Null alphabet in HmmPhyloEmissionProbabilities::setHmmStateAlphabet");
  if (!dynamic_cast<const HmmPhyloAlphabet*>(stateAlphabet.get()))
    throw HmmUnvalidAlphabetException("Non PhyloLikelihood alphabet in HmmPhyloEmissionProbabilities::setHmmStateAlphabet");

  phylAlph_ = dynamic_pointer_cast<HmmPhyloAlphabet>(stateAlphabet);

  std::vector<std::shared_ptr<Node_DF> > vEM;

  auto nbStates = phylAlph_->getNumberOfStates();

  for (size_t i = 0; i < nbStates; i++)
  {
    auto tmp = phylAlph_->alignedPhyloLikelihood(i).alignedLikelihoodCalculation().getSiteLikelihoods(false);
    vEM.push_back(tmp);
  }

  // Compound to put site log lik of different processes in a matrix

  emProb_ = EmissionLogk::create(context_, std::move(vEM), Dimension<MatrixLik>(Eigen::Index(nbStates), Eigen::Index(nbSites_)));
}
