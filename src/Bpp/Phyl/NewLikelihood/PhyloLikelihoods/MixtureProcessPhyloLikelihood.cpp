//
// File: MixtureProcessPhyloLikelihood.cpp
// Authors:
//   Laurent GuÃ©guen
// Created: vendredi 12 juillet 2013, Ã  14h 55
//

/*
  Copyright or Â© or Copr. Bio++ Development Team, (November 16, 2004)
  
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#include <Bpp/Numeric/VectorTools.h>

#include "MixtureProcessPhyloLikelihood.h"

using namespace std;
using namespace bpp;
using namespace numeric;

/******************************************************************************/

MixtureProcessPhyloLikelihood::MixtureProcessPhyloLikelihood(
  const AlignedValuesContainer& data,
  MixtureSequenceEvolution& processSeqEvol,
  CollectionNodes& collNodes,
  size_t nSeqEvol,
  size_t nData,
  bool verbose,
  bool patterns) :
  AbstractPhyloLikelihood(collNodes.getContext()),
  AbstractAlignedPhyloLikelihood(collNodes.getContext(), data.getNumberOfSites()),
  MultiProcessSequencePhyloLikelihood(data, processSeqEvol, collNodes, nSeqEvol, nData, verbose, patterns),
  mSeqEvol_(processSeqEvol),
  likCal_(new AlignedLikelihoodCalculation(collNodes.getContext()))
{
  if (vLikCal_.size() == 0)
    throw Exception("MixtureProcessPhyloLikelihood::MixtureProcessPhyloLikelihood : empty singleprocesslikelihoods set.");

  auto& simplex = mSeqEvol_.getSimplex();

  // parameters of the simplex
  const auto& param = simplex.getParameters();
  ParameterList paramList;

  for (size_t i = 0; i < param.size(); i++)
  {
    paramList.shareParameter(ConfiguredParameter::create(getContext(), param[i]));
  }

  shareParameters_(paramList);

  // make Simplex DF & Frequencies from it
  simplex_ = ConfiguredParametrizable::createConfigured<Simplex, ConfiguredSimplex>(getContext(), simplex, paramList, "");

  // for derivates
  auto deltaNode = NumericMutable<double>::create(getContext(), 0.001);
  auto config = NumericalDerivativeType::ThreePoints;

  simplex_->config.delta = deltaNode;
  simplex_->config.type = config;

  auto fsf = ConfiguredParametrizable::createRowVector<ConfiguredSimplex, FrequenciesFromSimplex, Eigen::RowVectorXd>(getContext(), {simplex_}, RowVectorDimension (Eigen::Index(simplex.dimension())));

  // get RowVectorXd for each single Calculation
  std::vector<std::shared_ptr<Node_DF> > vSL;

  for (auto& lik: vLikCal_)
  {
    vSL.push_back(lik->getSiteLikelihoods(true));
  }

  // put probabilities of the simplex
  vSL.push_back(fsf);

  auto single0 = vLikCal_[0];
  auto nbSite = single0->getNumberOfDistinctSites();

  auto sL = CWiseMean<RowLik, ReductionOf<RowLik>, Eigen::RowVectorXd>::create(getContext(), std::move(vSL), RowVectorDimension ((int)nbSite));

  likCal_->setSiteLikelihoods(sL, true);

  // likelihoods per site
  likCal_->setSiteLikelihoods(single0->expandVector(sL), false);

  auto su = SumOfLogarithms<RowLik>::create (getContext(), {sL, single0->getRootWeights()}, RowVectorDimension ((int)nbSite));

  likCal_->setLikelihoodNode(su);
}


/******************************************************************************/

VVdouble MixtureProcessPhyloLikelihood::getPosteriorProbabilitiesPerSitePerProcess() const
{
  size_t nbProcess = getNumberOfSubstitutionProcess();

  auto pb = getLikelihoodPerSitePerProcess();
  auto l = getLikelihoodPerSite();

  const auto& freq = simplex_->getTargetValue()->getFrequencies();

  for (size_t i = 0; i < nbSites_; ++i)
  {
    for (size_t j = 0; j < nbProcess; ++j)
    {
      pb[i][j] = convert(pb[i][j] * freq[j] / l[i]);
    }
  }
  return pb;
}
