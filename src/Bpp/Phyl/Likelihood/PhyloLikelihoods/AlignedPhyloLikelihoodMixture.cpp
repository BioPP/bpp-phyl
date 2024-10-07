// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "AlignedPhyloLikelihoodMixture.h"

using namespace bpp;
using namespace std;

AlignedPhyloLikelihoodMixture::AlignedPhyloLikelihoodMixture(
    Context& context,
    std::shared_ptr<PhyloLikelihoodContainer> pC,
    const std::vector<size_t>& nPhylo,
    bool inCollection) :
  AbstractPhyloLikelihood(context),
  AbstractParametrizable(""),
  AbstractPhyloLikelihoodSet(context, pC, nPhylo, inCollection),
  AbstractAlignedPhyloLikelihood(context, 0),
  AbstractAlignedPhyloLikelihoodSet(context, pC, nPhylo, inCollection, ""),
  likCal_(make_shared<AlignedLikelihoodCalculation>(context))
{
  Simplex simplex(getNumbersOfPhyloLikelihoods().size(), 1, false, "Mixture.");

  // parameters of the simplex
  const auto& param = simplex.getParameters();
  ParameterList paramList;

  for (size_t i = 0; i < param.size(); ++i)
  {
    paramList.shareParameter(ConfiguredParameter::create(this->context(), param[i]));
  }

  shareParameters_(paramList);

  // make Simplex DF & Frequencies from it

  simplex_ = ConfiguredParametrizable::createConfigured<Simplex, ConfiguredSimplex>(this->context(), simplex, paramList, "");

  // for derivates
  auto deltaNode = NumericMutable<double>::create(this->context(), 0.001);
  auto config = NumericalDerivativeType::ThreePoints;

  simplex_->config.delta = deltaNode;
  simplex_->config.type = config;

  auto fsf = ConfiguredParametrizable::createRowVector<ConfiguredSimplex, FrequenciesFromSimplex, Eigen::RowVectorXd>(this->context(), {simplex_}, RowVectorDimension (Eigen::Index(simplex.dimension())));

  // get RowVectorXd for each single Calculation
  std::vector<std::shared_ptr<Node_DF>> vSL;

  for (auto np : nPhylo)
  {
    vSL.push_back(alignedPhyloLikelihood(np).alignedLikelihoodCalculation().getSiteLikelihoods(false));
  }

  // put probabilities of the simplex

  vSL.push_back(fsf);

  auto sL = CWiseMean<RowLik, ReductionOf<RowLik>, Eigen::RowVectorXd>::create(this->context(), std::move(vSL), RowVectorDimension (Eigen::Index(nbSites_)));

  likCal_->setSiteLikelihoods(sL);

  auto su = SumOfLogarithms<RowLik>::create (this->context(), {sL}, RowVectorDimension (Eigen::Index (nbSites_)));

  likCal_->setLikelihoodNode(su);
}

Vdouble AlignedPhyloLikelihoodMixture::getPhyloProbabilities() const
{
  return accessValueConstCast<const Simplex*>(*simplex_)->getFrequencies();
}

double AlignedPhyloLikelihoodMixture::getPhyloProb(size_t index) const
{
  return accessValueConstCast<const Simplex*>(*simplex_)->prob(index);
}

/******************************************************************************/

void AlignedPhyloLikelihoodMixture::fireParameterChanged(const ParameterList& parameters)
{
  // simplex_->matchParametersValues(parameters);
  // SetOfAbstractPhyloLikelihood::fireParameterChanged(parameters);
}

void AlignedPhyloLikelihoodMixture::setPhyloProb(const Simplex& si)
{
  simplex_->matchParametersValues(si.getParameters());
  //  matchParametersValues(simplex_.getParameters());
}
