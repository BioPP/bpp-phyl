//
// File: SetOfAbstractPhyloLikelihood.cpp
// Authors:
//   Laurent GuÃ©guen
// Created: jeudi 14 mai 2015, Ã  17h 07
//

/*
  Copyright or Â© or Copr. Bio++ Development Team, (November 16, 2004)
  
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


#include "SetOfAbstractPhyloLikelihood.h"

using namespace bpp;
using namespace std;

SetOfAbstractPhyloLikelihood::SetOfAbstractPhyloLikelihood(Context& context, std::shared_ptr<PhyloLikelihoodContainer> pC, bool inCollection, const std::string& prefix) :
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

SetOfAbstractPhyloLikelihood::SetOfAbstractPhyloLikelihood(Context& context, std::shared_ptr<PhyloLikelihoodContainer> pC, const std::vector<size_t>& nPhylo, bool inCollection, const std::string& prefix) :
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


SetOfAbstractPhyloLikelihood::SetOfAbstractPhyloLikelihood(const SetOfAbstractPhyloLikelihood& sd) :
  AbstractPhyloLikelihood(sd),
  AbstractParametrizable(sd),
  pPhyloCont_(sd.pPhyloCont_),
  nPhylo_(sd.nPhylo_),
  vLikCal_(sd.vLikCal_)
{}

bool SetOfAbstractPhyloLikelihood::addPhyloLikelihood(size_t nPhyl, const std::string& suff)
{
  AbstractPhyloLikelihood* aPL = getAbstractPhyloLikelihood(nPhyl);

  if (aPL != NULL)
  {
    nPhylo_.push_back(nPhyl);
    if (suff != "") // Use specific parameters names
    {
      const auto& pl = aPL->getParameters();

      for (size_t i = 0; i < pl.size(); i++)
      {
        auto confP = dynamic_pointer_cast<ConfiguredParameter>(pl.getSharedParameter(i));
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

ParameterList SetOfAbstractPhyloLikelihood::getNonDerivableParameters() const
{
  ParameterList pl;
  for (size_t i = 0; i < nPhylo_.size(); i++)
  {
    pl.includeParameters(getPhyloLikelihood(nPhylo_[i])->getNonDerivableParameters());
  }

  return pl;
}

ParameterList SetOfAbstractPhyloLikelihood::getDerivableParameters() const
{
  ParameterList pl;
  for (size_t i = 0; i < nPhylo_.size(); i++)
  {
    pl.includeParameters(getPhyloLikelihood(nPhylo_[i])->getDerivableParameters());
  }

  return pl;
}


ParameterList SetOfAbstractPhyloLikelihood::getBranchLengthParameters() const
{
  ParameterList pl;
  for (size_t i = 0; i < nPhylo_.size(); i++)
  {
    pl.includeParameters(getPhyloLikelihood(nPhylo_[i])->getBranchLengthParameters());
  }

  return pl;
}

ParameterList SetOfAbstractPhyloLikelihood::getSubstitutionModelParameters() const
{
  ParameterList pl;
  for (size_t i = 0; i < nPhylo_.size(); i++)
  {
    pl.includeParameters(getPhyloLikelihood(nPhylo_[i])->getSubstitutionModelParameters());
  }

  return pl;
}

ParameterList SetOfAbstractPhyloLikelihood::getRateDistributionParameters() const
{
  ParameterList pl;
  for (size_t i = 0; i < nPhylo_.size(); i++)
  {
    pl.includeParameters(getPhyloLikelihood(nPhylo_[i])->getRateDistributionParameters());
  }

  return pl;
}

ParameterList SetOfAbstractPhyloLikelihood::getRootFrequenciesParameters() const
{
  ParameterList pl;
  for (size_t i = 0; i < nPhylo_.size(); i++)
  {
    pl.includeParameters(getPhyloLikelihood(nPhylo_[i])->getRootFrequenciesParameters());
  }

  return pl;
}
