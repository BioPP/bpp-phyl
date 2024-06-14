// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "PhyloLikelihoodSet.h"

using namespace bpp;
using namespace std;

AbstractPhyloLikelihoodSet::AbstractPhyloLikelihoodSet(
    Context& context,
    std::shared_ptr<PhyloLikelihoodContainer> pC,
    bool inCollection,
    const std::string& prefix) :
  AbstractPhyloLikelihood(context),
  AbstractParametrizable(prefix),
  pPhyloCont_(pC),
  nPhylo_(),
  vLikCal_()
{
  for (auto np:getPhyloContainer()->getNumbersOfPhyloLikelihoods())
  {
    addPhyloLikelihood(np, inCollection ? "" : "_" + TextTools::toString(np));
  }
}

AbstractPhyloLikelihoodSet::AbstractPhyloLikelihoodSet(
    Context& context,
    std::shared_ptr<PhyloLikelihoodContainer> pC,
    const std::vector<size_t>& nPhylo,
    bool inCollection,
    const std::string& prefix) :
  AbstractPhyloLikelihood(context),
  AbstractParametrizable(prefix),
  pPhyloCont_(pC),
  nPhylo_(),
  vLikCal_()
{
  for (auto np:nPhylo)
  {
    addPhyloLikelihood(np, inCollection ? "" : "_" + TextTools::toString(np));
  }
}

bool AbstractPhyloLikelihoodSet::addPhyloLikelihood(size_t nPhyl, const std::string& suff)
{
  auto aPL = getPhyloLikelihood(nPhyl);

  if (aPL)
  {
    nPhylo_.push_back(nPhyl);
    if (suff != "") // Use specific parameters names
    {
      const auto& pl = aPL->getParameters();

      for (size_t i = 0; i < pl.size(); i++)
      {
        auto confP = dynamic_pointer_cast<ConfiguredParameter>(pl.getParameter(i));
        if (confP == 0)
          throw Exception("SetOfAbstractPhyloLikelihood::addPhyloLikelihood: Parameter " + pl[i].getName() + " is not configured in PhyloLikelihood " + TextTools::toString(nPhyl));

        auto name = confP->getName() + suff;

        if (!hasParameter(name))
        {
          Parameter par(pl[i]);
          par.setName(name);

          auto cfPar = ConfiguredParameter::create(context_, {confP->dependency(0)}, par);
          shareParameter_(cfPar);
        }
      }
    }
    else
      shareParameters_(aPL->getParameters());


    vLikCal_.push_back(aPL->getLikelihoodCalculation());
    return true;
  }
  return false;
}

ParameterList AbstractPhyloLikelihoodSet::getNonDerivableParameters() const
{
  ParameterList pl;
  for (size_t i = 0; i < nPhylo_.size(); i++)
  {
    pl.includeParameters(getPhyloLikelihood(nPhylo_[i])->getNonDerivableParameters());
  }

  return pl;
}

ParameterList AbstractPhyloLikelihoodSet::getDerivableParameters() const
{
  ParameterList pl;
  for (size_t i = 0; i < nPhylo_.size(); i++)
  {
    pl.includeParameters(getPhyloLikelihood(nPhylo_[i])->getDerivableParameters());
  }

  return pl;
}


ParameterList AbstractPhyloLikelihoodSet::getBranchLengthParameters() const
{
  ParameterList pl;
  for (size_t i = 0; i < nPhylo_.size(); i++)
  {
    pl.includeParameters(getPhyloLikelihood(nPhylo_[i])->getBranchLengthParameters());
  }

  return pl;
}

ParameterList AbstractPhyloLikelihoodSet::getSubstitutionModelParameters() const
{
  ParameterList pl;
  for (size_t i = 0; i < nPhylo_.size(); i++)
  {
    pl.includeParameters(getPhyloLikelihood(nPhylo_[i])->getSubstitutionModelParameters());
  }

  return pl;
}

ParameterList AbstractPhyloLikelihoodSet::getRateDistributionParameters() const
{
  ParameterList pl;
  for (size_t i = 0; i < nPhylo_.size(); i++)
  {
    pl.includeParameters(getPhyloLikelihood(nPhylo_[i])->getRateDistributionParameters());
  }

  return pl;
}

ParameterList AbstractPhyloLikelihoodSet::getRootFrequenciesParameters() const
{
  ParameterList pl;
  for (size_t i = 0; i < nPhylo_.size(); i++)
  {
    pl.includeParameters(getPhyloLikelihood(nPhylo_[i])->getRootFrequenciesParameters());
  }

  return pl;
}
