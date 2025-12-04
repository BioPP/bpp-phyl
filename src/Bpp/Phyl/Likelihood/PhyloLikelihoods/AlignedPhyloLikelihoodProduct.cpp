// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "AlignedPhyloLikelihoodProduct.h"

using namespace bpp;
using namespace std;

AlignedPhyloLikelihoodProduct::AlignedPhyloLikelihoodProduct(
    Context& context,
    std::shared_ptr<PhyloLikelihoodContainer> pC,
    bool inCollection) :
  AbstractPhyloLikelihood(context),
  AbstractParametrizable(""),
  AbstractPhyloLikelihoodSet(context, pC, {}, inCollection),
  AbstractAlignedPhyloLikelihood(context, 0),
  AbstractAlignedPhyloLikelihoodSet(context, pC, inCollection),
  likCal_()
{
  auto nPhylo = pC->getNumbersOfPhyloLikelihoods();

  // get RowVectorXd for each single Calculation
  std::vector<std::shared_ptr<Node_DF>> vSL;

  for (auto np : nPhylo)
  {
    vSL.push_back(alignedPhyloLikelihood(np).alignedLikelihoodCalculation().getSiteLikelihoods(false));
  }

  auto sL = CWiseMul<RowLik, ReductionOf<RowLik>>::create(this->context(), std::move(vSL), RowVectorDimension (nbSites_));

  likCal_->setSiteLikelihoods(sL);

  auto su = SumOfLogarithms<RowLik>::create (this->context(), {sL}, RowVectorDimension (Eigen::Index (nbSites_)));

  likCal_->setLikelihoodNode(su);
}

AlignedPhyloLikelihoodProduct::AlignedPhyloLikelihoodProduct(
    Context& context,
    std::shared_ptr<PhyloLikelihoodContainer> pC,
    const std::vector<size_t>& nPhylo,
    bool inCollection) :
  AbstractPhyloLikelihood(context),
  AbstractParametrizable(""),
  AbstractPhyloLikelihoodSet(context, pC, {}, inCollection),
  AbstractAlignedPhyloLikelihood(context, 0),
  AbstractAlignedPhyloLikelihoodSet(context, pC, nPhylo, inCollection),
  likCal_(make_shared<AlignedLikelihoodCalculation>(context))
{
  // get RowVectorXd for each single Calculation
  std::vector<std::shared_ptr<Node_DF>> vSL;

  for (auto np : nPhylo)
  {
    vSL.push_back(alignedPhyloLikelihood(np).alignedLikelihoodCalculation().getSiteLikelihoods(false));
  }

  auto sL = CWiseMul<RowLik, ReductionOf<RowLik>>::create(this->context(), std::move(vSL), RowVectorDimension (nbSites_));

  likCal_->setSiteLikelihoods(sL);

  auto su = SumOfLogarithms<Eigen::RowVectorXd>::create (this->context(), {sL}, RowVectorDimension (nbSites_));

  likCal_->setLikelihoodNode(su);
}
