//
// File: SingleProcessPhyloLikelihood.cpp
// Authors:
//   FranÃÂ§ois Gindraud, Laurent GuÃÂ©guen
//   Creation: lundi 27 mai 2019, Ã 06h 35
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  
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


#include "SingleProcessPhyloLikelihood.h"

using namespace bpp;
using namespace std;
using namespace numeric;


Vdouble SingleProcessPhyloLikelihood::getPosteriorProbabilitiesForSitePerClass(size_t pos) const
{
  auto rates = getLikelihoodCalculationSingleProcess()->getSubstitutionProcess().getRateDistribution();

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
  auto rates = getLikelihoodCalculationSingleProcess()->getSubstitutionProcess().getRateDistribution();

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
  auto probas = getLikelihoodCalculationSingleProcess()->getSubstitutionProcess().getRateDistribution()->getProbabilities();
  auto rates = getLikelihoodCalculationSingleProcess()->getSubstitutionProcess().getRateDistribution()->getCategories();

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
  auto vv = getLikelihoodCalculationSingleProcess()->getLikelihoodsAtNode(nodeId)->getTargetValue();

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