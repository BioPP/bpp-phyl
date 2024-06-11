// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/VectorTools.h>

#include "MixtureProcessPhyloLikelihood.h"

using namespace std;
using namespace bpp;
using namespace numeric;

/******************************************************************************/

MixtureProcessPhyloLikelihood::MixtureProcessPhyloLikelihood(
    shared_ptr<const AlignmentDataInterface> data,
    shared_ptr<MixtureSequenceEvolution> processSeqEvol,
    shared_ptr<CollectionNodes> collNodes,
    size_t nSeqEvol,
    size_t nData) :
  AbstractPhyloLikelihood(collNodes->context()),
  AbstractAlignedPhyloLikelihood(collNodes->context(), data->getNumberOfSites()),
  AbstractSingleDataPhyloLikelihood(
      collNodes->context(),
      data->getNumberOfSites(),
      (processSeqEvol->getSubstitutionProcessNumbers().size() != 0)
          ? processSeqEvol->substitutionProcess(processSeqEvol->getSubstitutionProcessNumbers()[0]).getNumberOfStates()
    : 0,
      nData),
  AbstractSequencePhyloLikelihood(collNodes->context(), processSeqEvol, nData),
  AbstractParametrizable(""),
  MultiProcessSequencePhyloLikelihood(data, processSeqEvol, collNodes, nSeqEvol, nData),
  mSeqEvol_(processSeqEvol),
  likCal_(make_shared<AlignedLikelihoodCalculation>(collNodes->context()))
{
  if (vLikCal_.size() == 0)
    throw Exception("MixtureProcessPhyloLikelihood::MixtureProcessPhyloLikelihood : empty singleprocesslikelihoods set.");

  auto& simplex = mSeqEvol_->simplex();

  // parameters of the simplex
  const auto& param = simplex.getParameters();
  ParameterList paramList;

  for (size_t i = 0; i < param.size(); i++)
  {
    paramList.shareParameter(ConfiguredParameter::create(context(), param[i]));
  }

  shareParameters_(paramList);

  // make Simplex DF & Frequencies from it
  simplex_ = ConfiguredParametrizable::createConfigured<Simplex, ConfiguredSimplex>(context(), simplex, paramList, "");

  // for derivates
  auto deltaNode = NumericMutable<double>::create(context(), 0.001);
  auto config = NumericalDerivativeType::ThreePoints;

  simplex_->config.delta = deltaNode;
  simplex_->config.type = config;

  auto fsf = ConfiguredParametrizable::createRowVector<ConfiguredSimplex, FrequenciesFromSimplex, Eigen::RowVectorXd>(context(), {simplex_}, RowVectorDimension (Eigen::Index(simplex.dimension())));

  // get RowVectorXd for each single Calculation
  std::vector<std::shared_ptr<Node_DF>> vSL;

  for (auto& lik : vLikCal_)
  {
    vSL.push_back(lik->getSiteLikelihoods(true));
  }

  // put probabilities of the simplex
  vSL.push_back(fsf);

  auto single0 = vLikCal_[0];
  auto nbSite = single0->getNumberOfDistinctSites();

  auto sL = CWiseMean<RowLik, ReductionOf<RowLik>, Eigen::RowVectorXd>::create(context(), std::move(vSL), RowVectorDimension ((int)nbSite));

  likCal_->setSiteLikelihoods(sL, true);

  // likelihoods per site
  likCal_->setSiteLikelihoods(single0->expandVector(sL), false);

  auto su = SumOfLogarithms<RowLik>::create (context(), {sL, single0->getRootWeights()}, RowVectorDimension ((int)nbSite));

  likCal_->setLikelihoodNode(su);
}


/******************************************************************************/

VVdouble MixtureProcessPhyloLikelihood::getPosteriorProbabilitiesPerSitePerProcess() const
{
  size_t nbProcess = getNumberOfSubstitutionProcess();

  auto pb = getLikelihoodPerSitePerProcess();
  auto l = getLikelihoodPerSite();

  const auto& freq = simplex_->targetValue()->getFrequencies();

  for (size_t i = 0; i < nbSites_; ++i)
  {
    for (size_t j = 0; j < nbProcess; ++j)
    {
      pb[i][j] = pb[i][j] * freq[j] / convert(l[i]);
    }
  }
  return pb;
}
