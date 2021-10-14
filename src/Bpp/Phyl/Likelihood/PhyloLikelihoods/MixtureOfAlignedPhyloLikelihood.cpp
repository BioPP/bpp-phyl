//
// File: MixtureOfAlignedPhyloLikelihood.cpp
// Authors:
//   Laurent GuÃÂ©guen
// Created: lundi 26 octobre 2015, ÃÂ  21h 56
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


#include "MixtureOfAlignedPhyloLikelihood.h"

using namespace bpp;
using namespace std;

MixtureOfAlignedPhyloLikelihood::MixtureOfAlignedPhyloLikelihood(Context& context, std::shared_ptr<PhyloLikelihoodContainer> pC, const std::vector<size_t>& nPhylo, bool inCollection) :
  AbstractPhyloLikelihood(context),
  AbstractAlignedPhyloLikelihood(context, 0),
  SetOfAlignedPhyloLikelihood(context, pC, nPhylo, inCollection, ""),
  likCal_(new AlignedLikelihoodCalculation(context))
{
  Simplex simplex(getNumbersOfPhyloLikelihoods().size(), 1, false, "Mixture.");

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

  for (auto np:nPhylo)
  {
    vSL.push_back(getPhyloLikelihood(np)->getAlignedLikelihoodCalculation()->getSiteLikelihoods(false));
  }

  // put probabilities of the simplex

  vSL.push_back(fsf);

  auto sL = CWiseMean<RowLik, ReductionOf<RowLik>, Eigen::RowVectorXd>::create(getContext(), std::move(vSL), RowVectorDimension (Eigen::Index(nbSites_)));

  likCal_->setSiteLikelihoods(sL);

  auto su = SumOfLogarithms<RowLik>::create (getContext(), {sL}, RowVectorDimension (Eigen::Index (nbSites_)));

  likCal_->setLikelihoodNode(su);
}

MixtureOfAlignedPhyloLikelihood::MixtureOfAlignedPhyloLikelihood(const MixtureOfAlignedPhyloLikelihood& sd) :
  AbstractPhyloLikelihood(sd),
  AbstractAlignedPhyloLikelihood(sd),
  SetOfAlignedPhyloLikelihood(sd),
  simplex_(sd.simplex_),
  likCal_(sd.likCal_)
{}

Vdouble MixtureOfAlignedPhyloLikelihood::getPhyloProbabilities() const
{
  return accessValueConstCast<const Simplex*>(*simplex_)->getFrequencies();
}

double MixtureOfAlignedPhyloLikelihood::getPhyloProb(size_t index) const
{
  return accessValueConstCast<const Simplex*>(*simplex_)->prob(index);
}

/******************************************************************************/

void MixtureOfAlignedPhyloLikelihood::fireParameterChanged(const ParameterList& parameters)
{
  // simplex_->matchParametersValues(parameters);
  // SetOfAbstractPhyloLikelihood::fireParameterChanged(parameters);
}

void MixtureOfAlignedPhyloLikelihood::setPhyloProb(const Simplex& si)
{
  simplex_->matchParametersValues(si.getParameters());
//  matchParametersValues(simplex_.getParameters());
}
