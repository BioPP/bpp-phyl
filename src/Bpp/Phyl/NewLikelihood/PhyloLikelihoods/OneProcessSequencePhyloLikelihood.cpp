//
// File: OneProcessSequencePhyloLikelihood.cpp
// Created by: Julien Dutheil
// Created on: mardi 28 avril 2015, à 13h 11
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

#include "OneProcessSequencePhyloLikelihood.h"

using namespace std;
using namespace bpp;

/******************************************************************************/

OneProcessSequencePhyloLikelihood::OneProcessSequencePhyloLikelihood(
  Context& context,
  OneProcessSequenceEvolution& evol,
  size_t nSeqEvol,
  bool verbose,
  bool patterns) :
  AbstractPhyloLikelihood(context),
  AbstractAlignedPhyloLikelihood(context, 0),
  AbstractSequencePhyloLikelihood(context, evol, nSeqEvol),
  mSeqEvol_(evol),
  likCal_()
{
  resetParameters_();
  const auto& sp = evol.getSubstitutionProcess();
  likCal_ = std::make_shared<LikelihoodCalculationSingleProcess>(context, sp);

  shareParameters_(likCal_->getParameters());
}

/******************************************************************************/

OneProcessSequencePhyloLikelihood::OneProcessSequencePhyloLikelihood(
  Context& context,
  const AlignedValuesContainer& data,
  OneProcessSequenceEvolution& evol,
  size_t nSeqEvol,
  size_t nData,
  bool verbose,
  bool patterns) :
  AbstractPhyloLikelihood(context),
  AbstractAlignedPhyloLikelihood(context, data.getNumberOfSites()),
  AbstractSequencePhyloLikelihood(context, evol, nData),
  mSeqEvol_(evol),
  likCal_()
{
  resetParameters_();
  
  const auto& sp = evol.getSubstitutionProcess();
  likCal_ = std::make_shared<LikelihoodCalculationSingleProcess>(context, data, sp);  
  shareParameters_(likCal_->getParameters());
}

/******************************************************************************/

OneProcessSequencePhyloLikelihood::OneProcessSequencePhyloLikelihood(
  const AlignedValuesContainer& data,
  OneProcessSequenceEvolution& evol,
  CollectionNodes& collNodes,
  size_t nSeqEvol,
  size_t nData,
  bool verbose,
  bool patterns) :
  AbstractPhyloLikelihood(collNodes.getContext()),
  AbstractAlignedPhyloLikelihood(collNodes.getContext(), data.getNumberOfSites()),
  AbstractSequencePhyloLikelihood(collNodes.getContext(), evol, nData),
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
  copyEigenToBpp(getLikelihoodCalculationSingleProcess()->getSiteLikelihoodsForAllClasses(),vd);
  return vd;
}


/******************************************************************************/

Vdouble OneProcessSequencePhyloLikelihood::getPosteriorProbabilitiesForSitePerClass(size_t pos) const
{
  auto rates=getLikelihoodCalculationSingleProcess()->getSubstitutionProcess().getRateDistribution();
  
  if (!rates || rates->getNumberOfCategories()==1)
    return Vdouble(1,1);
  else
  {
    auto probas = rates->getProbabilities();
    Vdouble vv(rates->getNumberOfCategories());
    for (size_t i=0;i<vv.size();i++)
      vv[i]=probas[i] * (getLikelihoodCalculationSingleProcess()->getSiteLikelihoodsForAClass(i))(Eigen::Index(pos));
      
    vv/=VectorTools::sum(vv);
    return vv;
  }
}

/******************************************************************************/

VVdouble OneProcessSequencePhyloLikelihood::getPosteriorProbabilitiesPerSitePerClass() const
{
  auto rates=getLikelihoodCalculationSingleProcess()->getSubstitutionProcess().getRateDistribution();

  auto nbS=getLikelihoodCalculationSingleProcess()->getNumberOfSites();
  VVdouble vv(nbS);

  if (!rates || rates->getNumberOfCategories()==1)
  {
    for (auto& v:vv)
      v.resize(1,1);
  }
  else
  {
    Eigen::VectorXd probas;
    copyBppToEigen(rates->getProbabilities(),probas);
    
    auto vvLik=getLikelihoodCalculationSingleProcess()->getSiteLikelihoodsForAllClasses();
    for (size_t i=0;i<nbS;i++)
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
  auto probas = getLikelihoodCalculationSingleProcess()->getSubstitutionProcess().getRateDistribution()->getProbabilities();
  auto rates = getLikelihoodCalculationSingleProcess()->getSubstitutionProcess().getRateDistribution()->getCategories();
  
  size_t nbSites   = getNumberOfSites();
  size_t nbClasses = getNumberOfClasses();
  VVdouble pb = getLikelihoodPerSitePerClass();
  Vdouble l  = getLikelihoodPerSite();
  Vdouble prates(nbSites, 0.);
  for (size_t i = 0; i < nbSites; i++)
  {
    for (size_t j = 0; j < nbClasses; j++)
    {
      prates[i] += (pb[i][j] / l[i]) * probas[j] *  rates[j];
    }
  }
  return prates;
}
  
/******************************************************************************/


Vdouble OneProcessSequencePhyloLikelihood::getPosteriorStateFrequencies(uint nodeId)
{
  auto vv = getLikelihoodCalculationSingleProcess()->getLikelihoodsAtNode(nodeId)->getTargetValue();
  
  for (auto i=0;i<vv.cols();i++)
    vv.col(i)/=vv.col(i).sum();

  VectorLik vs(vv.rowwise().mean());
    
  Vdouble v;
  copyEigenToBpp(vs, v);
  return v;
}
