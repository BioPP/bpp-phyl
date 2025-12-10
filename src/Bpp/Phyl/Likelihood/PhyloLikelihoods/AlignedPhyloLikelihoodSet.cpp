// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "AlignedPhyloLikelihoodSet.h"

using namespace std;
using namespace bpp;


AbstractAlignedPhyloLikelihoodSet::AbstractAlignedPhyloLikelihoodSet(
    Context& context,
    std::shared_ptr<PhyloLikelihoodContainer> pC,
    bool inCollection,
    const std::string& prefix) :
  AbstractPhyloLikelihood(context),
  AbstractPhyloLikelihoodSet(context, pC, inCollection, prefix),
  AbstractAlignedPhyloLikelihood(context, 0)
{
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wuninitialized" // Remove warning (TODO check whether this is a real one!)
  for (auto np:nPhylo_)
  {
    auto aPL = getAlignedPhyloLikelihood(np);

    if (!aPL)
      throw Exception("AlignedPhyloLikelihoodSet::AlignedPhyloLikelihoodSet  :non aligned PhyloLikelihood: " + TextTools::toString(np));
    if (getNumberOfSites() == 0)
      setNumberOfSites(aPL->getNumberOfSites());
    else if (aPL->getNumberOfSites() != getNumberOfSites())
      throw BadSizeException("AlignedPhyloLikelihoodSet::AlignedPhyloLikelihoodSet: mismatch lengths between aligned PhyloLikelihood: ", aPL->getNumberOfSites(), getNumberOfSites());
  }
#pragma GCC diagnostic pop
}

/*************************************************************/

AbstractAlignedPhyloLikelihoodSet::AbstractAlignedPhyloLikelihoodSet(
    Context& context,
    std::shared_ptr<PhyloLikelihoodContainer> pC,
    const std::vector<size_t>& nPhylo,
    bool inCollection,
    const std::string& prefix) :
  AbstractPhyloLikelihood(context),
  AbstractPhyloLikelihoodSet(context, pC, nPhylo, inCollection, prefix),
  AbstractAlignedPhyloLikelihood(context, 0)
{
  for (auto np : nPhylo)
  {
    auto aPL = getAlignedPhyloLikelihood(np);

    if (!aPL)
      throw Exception("AlignedPhyloLikelihoodSet::AlignedPhyloLikelihoodSet  :non aligned PhyloLikelihood: " + TextTools::toString(np));
    if (getNumberOfSites() == 0)
      setNumberOfSites(aPL->getNumberOfSites());
    else if (aPL->getNumberOfSites() != getNumberOfSites())
      throw BadSizeException("AlignedPhyloLikelihoodSet::AlignedPhyloLikelihoodSet: mismatch lengths between aligned PhyloLikelihood: ", aPL->getNumberOfSites(), getNumberOfSites());
  }
}

/******************************************************************************/

bool AbstractAlignedPhyloLikelihoodSet::addPhyloLikelihood(size_t nPhyl, const std::string& suff)
{
  auto aPL = getAlignedPhyloLikelihood(nPhyl);

  if (aPL && (getNumberOfSites() == 0 || aPL->getNumberOfSites() == getNumberOfSites()))
  {
    if (AbstractPhyloLikelihoodSet::addPhyloLikelihood(nPhyl, suff))
    {
      if (getNumberOfSites() == 0)
        setNumberOfSites(aPL->getNumberOfSites());
      return true;
    }
    else
      return false;
  }

  return false;
}

/******************************************************************************/
