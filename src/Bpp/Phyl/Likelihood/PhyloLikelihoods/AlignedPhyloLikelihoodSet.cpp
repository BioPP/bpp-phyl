//
// File: AlignedPhyloLikelihoodSet.cpp
// Authors:
//   Laurent GuÃÂ©guen
// Created: mercredi 7 octobre 2015, ÃÂ  14h 21
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

