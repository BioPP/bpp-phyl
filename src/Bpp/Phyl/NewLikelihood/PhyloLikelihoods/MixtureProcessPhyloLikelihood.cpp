//
// File: MixtureProcessPhyloLikelihood.cpp
// Created by: Laurent Guéguen
// Created on: vendredi 12 juillet 2013, à 14h 55
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

#include "MixtureProcessPhyloLikelihood.h"

#include <Bpp/Numeric/VectorTools.h>

using namespace std;
using namespace bpp;

/******************************************************************************/

MixtureProcessPhyloLikelihood::MixtureProcessPhyloLikelihood(
  Context& context,
  const AlignedValuesContainer& data,
  MixtureSequenceEvolution& processSeqEvol,
  size_t nSeqEvol,
  size_t nData,
  bool verbose,
  bool patterns) :
  AbstractPhyloLikelihood(context),
  AbstractAlignedPhyloLikelihood(context, data.getNumberOfSites()),
  MultiProcessSequencePhyloLikelihood(context, data, processSeqEvol, nSeqEvol, nData, verbose, patterns),
  mSeqEvol_(processSeqEvol)
{

  if (likCal_.getNumberOfSingleProcess()==0)
    throw Exception("MixtureProcessPhyloLikelihood::MixtureProcessPhyloLikelihood : empty singleprocesslikelihoods set.");

  auto& simplex = mSeqEvol_.getSimplex();
  
  //parameters of the simplex
  const auto& param=simplex.getParameters();
  ParameterList paramList;

  for (size_t i=0;i<param.size();i++)
    paramList.shareParameter(ConfiguredParameter::create(getContext(), param[i]));
  
  shareParameters_(paramList);

  // make Simplex DF & Frequencies from it
  simplex_ = ConfiguredParametrizable::createConfigured<Simplex, ConfiguredSimplex>(getContext(), simplex, paramList, "");

  auto fsf = ConfiguredParametrizable::createVector<ConfiguredSimplex, FrequenciesFromSimplex>(getContext(), {simplex_}, rowVectorDimension (Eigen::Index(simplex.dimension())));

  // get RowVectorXd for each single Calculation
  std::vector<std::shared_ptr<Node_DF>> vSL;
  
  for (size_t i=0; i< likCal_.getNumberOfSingleProcess() ; i++)
    vSL.push_back(likCal_.getSingleLikelihood(i)->getSiteLikelihoods(true));
  
  // put probabilities of the simplex
  vSL.push_back(fsf);

  auto single0 = likCal_.getSingleLikelihood(0);
  auto nbSite = single0->getNumberOfDistinctSites();

  siteLikelihoods_ = CWiseMean<Eigen::RowVectorXd, ReductionOf<Eigen::RowVectorXd>, Eigen::RowVectorXd>::create(getContext(), std::move(vSL), rowVectorDimension (Eigen::Index(nbSite)));

  // likelihoods per site
  patternedSiteLikelihoods_ = single0->expandVector(siteLikelihoods_);

  auto su = SumOfLogarithms<Eigen::RowVectorXd>::create (getContext(), {siteLikelihoods_, single0->getRootWeights()}, rowVectorDimension (Eigen::Index (nbSite)));

  likCal_.setLikelihoodNode(su);
}

/******************************************************************************/

double MixtureProcessPhyloLikelihood::getLikelihoodForASite(size_t site) const
{
  return patternedSiteLikelihoods_->getTargetValue()(site);
}

/******************************************************************************/

VVdouble MixtureProcessPhyloLikelihood::getPosteriorProbabilitiesPerSitePerProcess() const
{
  size_t nbProcess = getNumberOfSubstitutionProcess();

  VVdouble pb = getLikelihoodPerSitePerProcess();
  Vdouble l = getLikelihoodPerSite();

  const auto& freq = simplex_->getTargetValue()->getFrequencies();
  
  for (size_t i = 0; i < nbSites_; ++i)
  {
    for (size_t j = 0; j < nbProcess; ++j)
    {
      pb[i][j] = pb[i][j] * freq[j] / l[i];
    }
  }
  return pb;
}


