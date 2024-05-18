// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "PatternTools.h"

// From bpp-seq:
#include <Bpp/Seq/SiteTools.h>

using namespace bpp;

// From the STL:
#include <iostream>
#include <algorithm>

using namespace std;

/******************************************************************************/

std::unique_ptr<AlignmentDataInterface> PatternTools::getSequenceSubset(
    const AlignmentDataInterface& sequenceSet,
    const Node& node)
{
  try
  {
    const auto& siteContainer = dynamic_cast<const SiteContainerInterface&>(sequenceSet);
    return getSequenceSubset(siteContainer, node);
  }
  catch (std::bad_cast& e)
  {}

  try
  {
    const auto& siteContainer = dynamic_cast<const ProbabilisticSiteContainerInterface&>(sequenceSet);
    return getSequenceSubset(siteContainer, node);
  }
  catch (std::bad_cast& e)
  {}

  throw Exception("PatternTools::getSequenceSubset : unsupported sequence type.");
}

std::unique_ptr<SiteContainerInterface> PatternTools::getSequenceSubset(
    const SiteContainerInterface& sequenceSet,
    const Node& node)
{
  auto alphabet = sequenceSet.getAlphabet();
  size_t nbSites = sequenceSet.getNumberOfSites();
  auto sequenceSubset = std::make_unique<VectorSiteContainer>(alphabet);

  auto leaves = TreeTemplateTools::getLeaves(node);

  for (auto i : leaves)
  {
    if (i->hasName())
    {
      // Use sequence name as key.
      try
      {
        auto newSeq = std::make_unique<Sequence>(sequenceSet.sequence(i->getName()));
        sequenceSubset->addSequence(i->getName(), newSeq);
      }
      catch (std::exception& e)
      {
        ApplicationTools::displayWarning("PatternTools::getSequenceSubset : Leaf name not found in sequence file: " + i->getName() + " : Replaced with unknown sequence");

        auto seq = std::make_unique<Sequence>(i->getName(), "", alphabet);
        seq->setToSizeR(nbSites);
        SymbolListTools::changeGapsToUnknownCharacters(*seq);
        sequenceSubset->addSequence(i->getName(), seq);
      }
    }
  }
  sequenceSubset->setSiteCoordinates(sequenceSet.getSiteCoordinates());
  return sequenceSubset;
}

std::unique_ptr<ProbabilisticSiteContainerInterface> PatternTools::getSequenceSubset(
    const ProbabilisticSiteContainerInterface& sequenceSet,
    const Node& node)
{
  auto alphabet = sequenceSet.getAlphabet();
  size_t nbSites = sequenceSet.getNumberOfSites();
  auto sequenceSubset = std::make_unique<ProbabilisticVectorSiteContainer>(alphabet);

  auto leaves = TreeTemplateTools::getLeaves(node);

  for (auto i : leaves)
  {
    if (i->hasName())
    {
      // Use sequence name as key.
      try
      {
        auto newSeq = std::make_unique<ProbabilisticSequence>(sequenceSet.sequence(i->getName()));
        sequenceSubset->addSequence(i->getName(), newSeq);
      }
      catch (std::exception const& e)
      {
        ApplicationTools::displayWarning("PatternTools::getSequenceSubset : Leaf name not found in sequence file: " + i->getName() + " : Replaced with unknown sequence");

        auto newSeq = std::make_unique<ProbabilisticSequence>(i->getName(), Table<double>(alphabet->getSize(), 0), alphabet);
        newSeq->setToSizeR(nbSites);
        SymbolListTools::changeGapsToUnknownCharacters(*newSeq);
        sequenceSubset->addSequence(i->getName(), newSeq);
      }
    }
  }
  sequenceSubset->setSiteCoordinates(sequenceSet.getSiteCoordinates());
  return sequenceSubset;
}

/******************************************************************************/

std::unique_ptr<AlignmentDataInterface> PatternTools::getSequenceSubset(
    const AlignmentDataInterface& sequenceSet,
    const std::vector<std::string>& names)
{
  auto alphabet = sequenceSet.getAlphabet();

  try
  {
    const auto& siteContainer = dynamic_cast<const SiteContainerInterface&>(sequenceSet);
    return getSequenceSubset(siteContainer, names);
  }
  catch (std::bad_cast& e)
  {}

  try
  {
    const auto& siteContainer = dynamic_cast<const ProbabilisticSiteContainerInterface&>(sequenceSet);
    return getSequenceSubset(siteContainer, names);
  }
  catch (std::bad_cast& e)
  {}

  throw Exception("PatternTools::getSequenceSubset : unsupported sequence type.");
}


std::unique_ptr<SiteContainerInterface> PatternTools::getSequenceSubset(
    const SiteContainerInterface& sequenceSet,
    const std::vector<std::string>& names)
{
  auto alphabet = sequenceSet.getAlphabet();
  auto sequenceSubset = std::make_unique<VectorSiteContainer>(alphabet);

  for (auto& i : names)
  {
    if (sequenceSet.hasSequence(i))
    {
      auto newSeq = std::make_unique<Sequence>(sequenceSet.sequence(i));
      sequenceSubset->addSequence(i, newSeq);
    }
    else
      throw SequenceNotFoundException("PatternTools ERROR: name not found in sequence container: ", i);
  }
  sequenceSubset->setSiteCoordinates(sequenceSet.getSiteCoordinates());
  return std::move(sequenceSubset);
}

std::unique_ptr<ProbabilisticSiteContainerInterface> PatternTools::getSequenceSubset(
    const ProbabilisticSiteContainerInterface& sequenceSet,
    const std::vector<std::string>& names)
{
  auto alphabet = sequenceSet.getAlphabet();
  auto sequenceSubset = std::make_unique<ProbabilisticVectorSiteContainer>(alphabet);

  for (auto& i : names)
  {
    if (sequenceSet.hasSequence(i))
    {
      auto newSeq = std::make_unique<ProbabilisticSequence>(sequenceSet.sequence(i));
      sequenceSubset->addSequence(i, newSeq);
    }
    else
      throw SequenceNotFoundException("PatternTools ERROR: name not found in sequence container: ", i);
  }
  sequenceSubset->setSiteCoordinates(sequenceSet.getSiteCoordinates());
  return std::move(sequenceSubset);
}

/******************************************************************************/

std::unique_ptr<AlignmentDataInterface> PatternTools::shrinkSiteSet(
    const AlignmentDataInterface& siteSet)
{
  try
  {
    const auto& sc = dynamic_cast<const SiteContainerInterface&>(siteSet);
    return shrinkSiteSet(sc);
  }
  catch (std::bad_cast& e)
  {}

  try
  {
    const auto& psc = dynamic_cast<const ProbabilisticSiteContainerInterface&>(siteSet);
    return shrinkSiteSet(psc);
  }
  catch (std::bad_cast& e)
  {}

  throw Exception("PatternTools::shrinkSiteSet : unsupported sequence type.");
}

std::unique_ptr<SiteContainerInterface> PatternTools::shrinkSiteSet(
    const SiteContainerInterface& siteSet)
{
  auto alphabet = siteSet.getAlphabet();

  if (siteSet.getNumberOfSites() == 0)
    throw Exception("PatternTools::shrinkSiteSet siteSet is void.");

  vector<std::unique_ptr<Site>> sites;

  for (unsigned int i = 0; i < siteSet.getNumberOfSites(); ++i)
  {
    const auto& currentSite = siteSet.site(i);
    bool siteExists = false;
    for (unsigned int j = 0; !siteExists && j < sites.size(); ++j)
    {
      if (SiteTools::areSymbolListsIdentical(currentSite, *sites[j]))
        siteExists = true;
    }
    if (!siteExists)
      sites.push_back(make_unique<Site>(currentSite));
  }
  auto result = make_unique<VectorSiteContainer>(sites, alphabet, false);
  result->setSequenceNames(siteSet.getSequenceNames(), true); // Update keys too
  return result;
}

std::unique_ptr<ProbabilisticSiteContainerInterface> PatternTools::shrinkSiteSet(
    const ProbabilisticSiteContainerInterface& siteSet)
{
  auto alphabet = siteSet.getAlphabet();

  if (siteSet.getNumberOfSites() == 0)
    throw Exception("PatternTools::shrinkSiteSet siteSet is void.");

  vector<unique_ptr<ProbabilisticSite>> sites;

  for (unsigned int i = 0; i < siteSet.getNumberOfSites(); ++i)
  {
    const auto& currentSite = siteSet.site(i);
    bool siteExists = false;
    for (unsigned int j = 0; !siteExists && j < sites.size(); ++j)
    {
      if (SiteTools::areSymbolListsIdentical(currentSite, *sites[j]))
        siteExists = true;
    }
    if (!siteExists)
      sites.push_back(make_unique<ProbabilisticSite>(currentSite));
  }
  auto result = make_unique<ProbabilisticVectorSiteContainer>(sites, alphabet, false);
  result->setSequenceNames(siteSet.getSequenceNames(), false);
  return result;
}

/******************************************************************************/

Vint PatternTools::getIndexes(
    const AlignmentDataInterface& sequences1,
    const AlignmentDataInterface& sequences2)
{
  size_t nbSites = sequences1.getNumberOfSites();
  Vint indexes(nbSites);

  try
  {
    const auto& sc1 = dynamic_cast<const SiteContainerInterface&>(sequences1);
    const auto& sc2 = dynamic_cast<const SiteContainerInterface&>(sequences2);

    for (size_t i = 0; i < nbSites; ++i)
    {
      // For each site in sequences1,
      indexes[i] = -1;
      const auto& site1 = dynamic_cast<const Site&>(sc1.site(i));
      for (size_t j = 0; j < sequences2.getNumberOfSites(); ++j)
      {
        const auto& site2 = dynamic_cast<const Site&>(sc2.site(i));
        if (SiteTools::areSymbolListsIdentical(site1, site2))
        {
          indexes[i] = static_cast<int>(j);
          break;
        }
      }
    }
    return indexes;
  }
  catch (std::bad_cast& e)
  {}

  try
  {
    const auto& psc1 = dynamic_cast<const ProbabilisticSiteContainerInterface&>(sequences1);
    const auto& psc2 = dynamic_cast<const ProbabilisticSiteContainerInterface&>(sequences2);

    for (size_t i = 0; i < nbSites; ++i)
    {
      // For each site in sequences1,
      indexes[i] = -1;
      const auto& site1 = dynamic_cast<const ProbabilisticSite&>(psc1.site(i));
      for (size_t j = 0; j < sequences2.getNumberOfSites(); ++j)
      {
        const auto& site2 = dynamic_cast<const ProbabilisticSite&>(psc2.site(i));
        if (SiteTools::areSymbolListsIdentical(site1, site2))
        {
          indexes[i] = static_cast<int>(j);
          break;
        }
      }
    }
    return indexes;
  }
  catch (std::bad_cast& e)
  {}

  throw Exception("PatternTools::shrinkSiteSet : unsupported sequence type.");
}

/******************************************************************************/
