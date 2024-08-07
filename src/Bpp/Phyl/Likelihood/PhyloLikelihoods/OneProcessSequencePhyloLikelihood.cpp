// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "OneProcessSequencePhyloLikelihood.h"

using namespace std;
using namespace bpp;
using namespace numeric;

/******************************************************************************/

OneProcessSequencePhyloLikelihood::OneProcessSequencePhyloLikelihood(
    Context& context,
    std::shared_ptr<OneProcessSequenceEvolution> evol,
    size_t nSeqEvol) :
  AbstractPhyloLikelihood(context),
  AbstractAlignedPhyloLikelihood(context, 0),
  AbstractSingleDataPhyloLikelihood(context, 0, (evol->getSubstitutionProcessNumbers().size() != 0) ? evol->substitutionProcess(evol->getSubstitutionProcessNumbers()[0]).getNumberOfStates() : 0, 0),
  AbstractSequencePhyloLikelihood(context, evol, nSeqEvol),
  AbstractParametrizable(""),
  AbstractParametrizableSequencePhyloLikelihood(context, evol, nSeqEvol),
  mSeqEvol_(evol),
  likCal_()
{
  resetParameters_();

  auto sp = evol->getSubstitutionProcess();
  likCal_ = std::make_shared<LikelihoodCalculationSingleProcess>(context, sp);

  shareParameters_(likCal_->getParameters());
}

/******************************************************************************/

OneProcessSequencePhyloLikelihood::OneProcessSequencePhyloLikelihood(
    Context& context,
    std::shared_ptr<const AlignmentDataInterface> data,
    std::shared_ptr<OneProcessSequenceEvolution> evol,
    size_t nSeqEvol,
    size_t nData) :
  AbstractPhyloLikelihood(context),
  AbstractAlignedPhyloLikelihood(context, data->getNumberOfSites()),
  AbstractSingleDataPhyloLikelihood(context, data->getNumberOfSites(), (evol->getSubstitutionProcessNumbers().size() != 0) ? evol->substitutionProcess(evol->getSubstitutionProcessNumbers()[0]).getNumberOfStates() : 0, nData),
  AbstractSequencePhyloLikelihood(context, evol, nData),
  AbstractParametrizable(""),
  AbstractParametrizableSequencePhyloLikelihood(context, evol, nSeqEvol),
  mSeqEvol_(evol),
  likCal_()
{
  resetParameters_();

  const auto& sp = evol->getSubstitutionProcess();
  likCal_ = std::make_shared<LikelihoodCalculationSingleProcess>(context, data, sp);
  shareParameters_(likCal_->getParameters());
}

/******************************************************************************/

OneProcessSequencePhyloLikelihood::OneProcessSequencePhyloLikelihood(
    std::shared_ptr<const AlignmentDataInterface> data,
    std::shared_ptr<OneProcessSequenceEvolution> evol,
    std::shared_ptr<CollectionNodes> collNodes,
    size_t nSeqEvol,
    size_t nData) :
  AbstractPhyloLikelihood(collNodes->context()),
  AbstractAlignedPhyloLikelihood(collNodes->context(), data->getNumberOfSites()),
  AbstractSingleDataPhyloLikelihood(collNodes->context(), data->getNumberOfSites(), (evol->getSubstitutionProcessNumbers().size() != 0) ? evol->substitutionProcess(evol->getSubstitutionProcessNumbers()[0]).getNumberOfStates() : 0, nData),
  AbstractSequencePhyloLikelihood(collNodes->context(), evol, nData),
  AbstractParametrizable(""),
  AbstractParametrizableSequencePhyloLikelihood(collNodes->context(), evol, nSeqEvol),
  mSeqEvol_(evol),
  likCal_()
{
  resetParameters_();

  likCal_ = std::make_shared<LikelihoodCalculationSingleProcess>(collNodes, data, nSeqEvol);

  shareParameters_(likCal_->getParameters());
}


/******************************************************************************/

VVdouble OneProcessSequencePhyloLikelihood::getLikelihoodPerSitePerClass() const
{
  VVdouble vd;
  copyEigenToBpp(getLikelihoodCalculationSingleProcess()->getSiteLikelihoodsForAllClasses(), vd);
  return vd;
}


/******************************************************************************/

Vdouble OneProcessSequencePhyloLikelihood::getPosteriorProbabilitiesForSitePerClass(size_t pos) const
{
  auto rates = getLikelihoodCalculationSingleProcess()->substitutionProcess().getRateDistribution();

  if (!rates || rates->getNumberOfCategories() == 1)
    return Vdouble(1, 1);
  else
  {
    auto probas = rates->getProbabilities();
    std::vector<DataLik> vv(rates->getNumberOfCategories());
    for (size_t i = 0; i < vv.size(); i++)
    {
      vv[i] = probas[i] * (getLikelihoodCalculationSingleProcess()->getSiteLikelihoodsForAClass(i))(Eigen::Index(pos));
    }

    auto sv = VectorTools::sum(vv);
    Vdouble res(rates->getNumberOfCategories());

    for (size_t i = 0; i < vv.size(); i++)
    {
      res[i] = convert(vv[i] / sv);
    }
    return res;
  }
}

/******************************************************************************/

VVdouble OneProcessSequencePhyloLikelihood::getPosteriorProbabilitiesPerSitePerClass() const
{
  auto rates = getLikelihoodCalculationSingleProcess()->substitutionProcess().getRateDistribution();

  auto nbS = getLikelihoodCalculationSingleProcess()->getNumberOfSites();
  VVdouble vv(nbS);

  if (!rates || rates->getNumberOfCategories() == 1)
  {
    for (auto& v:vv)
    {
      v.resize(1, 1);
    }
  }
  else
  {
    Eigen::VectorXd probas;
    copyBppToEigen(rates->getProbabilities(), probas);

    auto vvLik = getLikelihoodCalculationSingleProcess()->getSiteLikelihoodsForAllClasses();
    for (size_t i = 0; i < nbS; i++)
    {
      vv[i].resize(size_t(vvLik.rows()));
      VectorLik sv(numeric::cwise(vvLik.col(Eigen::Index(i))) * probas.array());
      sv /= sv.sum();
      copyEigenToBpp(sv, vv[i]);
    }
  }
  return vv;
}

/******************************************************************************/

vector<size_t> OneProcessSequencePhyloLikelihood::getClassWithMaxPostProbPerSite() const
{
  size_t nbSites = getNumberOfSites();
  auto l(getLikelihoodPerSitePerClass());
  vector<size_t> classes(nbSites);
  for (size_t i = 0; i < nbSites; ++i)
  {
    classes[i] = VectorTools::whichMax<double>(l[i]);
  }
  return classes;
}


/******************************************************************************/

Vdouble OneProcessSequencePhyloLikelihood::getPosteriorRatePerSite() const
{
  auto probas = likelihoodCalculationSingleProcess().substitutionProcess().rateDistribution().getProbabilities();
  auto rates = likelihoodCalculationSingleProcess().substitutionProcess().rateDistribution().getCategories();

  size_t nbSites   = getNumberOfSites();
  size_t nbClasses = getNumberOfClasses();
  auto pb = getLikelihoodPerSitePerClass();
  auto l  = getLikelihoodPerSite();
  Vdouble prates(nbSites, 0.);
  for (size_t i = 0; i < nbSites; i++)
  {
    for (size_t j = 0; j < nbClasses; j++)
    {
      prates[i] += convert((pb[i][j] / l[i]) * probas[j] *  rates[j]);
    }
  }
  return prates;
}

/******************************************************************************/


Vdouble OneProcessSequencePhyloLikelihood::getPosteriorStateFrequencies(uint nodeId)
{
  auto vv = getLikelihoodCalculationSingleProcess()->getLikelihoodsAtNode(nodeId)->targetValue();

  size_t nbSites   = getNumberOfSites();
  VVdouble pp;
  pp.resize(nbSites);

  for (uint i = 0; i < (uint)nbSites; i++)
  {
    copyEigenToBpp(vv.col(i) / vv.col(i).sum(), pp[size_t(i)]);
  }

  Vdouble v(nbStates_);
  for (uint st = 0; st < (uint)nbStates_; st++)
  {
    auto s = 0.0;
    for (uint i = 0; i < (uint)nbSites; i++)
    {
      s += pp[(size_t)i][size_t(st)];
    }

    v[size_t(st)] = s / (double)nbSites;
  }
  return v;
}
