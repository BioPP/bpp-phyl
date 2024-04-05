// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "SingleProcessPhyloLikelihood.h"

using namespace bpp;
using namespace std;
using namespace numeric;


Vdouble SingleProcessPhyloLikelihood::getPosteriorProbabilitiesForSitePerClass(size_t pos) const
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

VVdouble SingleProcessPhyloLikelihood::getPosteriorProbabilitiesPerSitePerClass() const
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

VVDataLik SingleProcessPhyloLikelihood::getLikelihoodPerSitePerClass() const
{
  VVDataLik vd;
  auto eg = getLikelihoodCalculationSingleProcess()->getSiteLikelihoodsForAllClasses();
  copyEigenToBpp(eg.transpose(), vd);
  return vd;
}

/******************************************************************************/

vector<size_t> SingleProcessPhyloLikelihood::getClassWithMaxPostProbPerSite() const
{
  size_t nbSites = getNumberOfSites();
  auto l = getLikelihoodPerSitePerClass();
  vector<size_t> classes(nbSites);
  for (size_t i = 0; i < nbSites; ++i)
  {
    classes[i] = VectorTools::whichMax<DataLik>(l[i]);
  }
  return classes;
}

/******************************************************************************/

Vdouble SingleProcessPhyloLikelihood::getPosteriorRatePerSite() const
{
  auto probas = getLikelihoodCalculationSingleProcess()->substitutionProcess().getRateDistribution()->getProbabilities();
  auto rates = getLikelihoodCalculationSingleProcess()->substitutionProcess().getRateDistribution()->getCategories();

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

Vdouble SingleProcessPhyloLikelihood::getPosteriorStateFrequencies(uint nodeId)
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
    for (size_t i = 0; i < (size_t)nbSites; i++)
    {
      s += pp[i][size_t(st)];
    }

    v[size_t(st)] = s / (double)nbSites;
  }
  return v;
}
