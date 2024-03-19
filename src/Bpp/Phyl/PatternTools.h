// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_PATTERNTOOLS_H
#define BPP_PHYL_PATTERNTOOLS_H

#include <Bpp/Numeric/VectorTools.h>

#include "Tree/PhyloTree.h"
#include "Tree/TreeTemplateTools.h"

// From SeqLib:
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/SymbolListTools.h>
#include <Bpp/Seq/Site.h>
#include <Bpp/Seq/Container/SiteContainer.h>
#include <Bpp/Seq/Container/AlignmentData.h>

// From the STL:
#include <map>

namespace bpp
{
/**
 * @brief Utilitary methods to compute site patterns.
 *
 * Theses methods are mainly designed to save computation in likelihood
 * and parsimony methods.
 */
class PatternTools
{
public:
  /**
   * @brief Extract the sequences corresponding to a given subtree.
   *
   * @param sequenceSet The container to look in.
   * @param node        The root node of the subtree to check.
   * @param tree        The tree owing the node.
   * @return A new site container with corresponding sequences.
   * @throw Exception if an error occured.
   */
  template<class N, class E, class I>
  static std::unique_ptr<AlignmentDataInterface> getSequenceSubset(
    const AlignmentDataInterface& sequenceSet,
    const std::shared_ptr<N> node,
    const AssociationTreeGraphImplObserver<N, E, I>& tree);

  /**
   * @brief Extract the sequences corresponding to a given subtree.
   *
   * @param sequenceSet The container to look in.
   * @param node        The root node of the subtree to check.
   * @return A new site container with corresponding sequences.
   * @throw Exception if an error occured.
   */
  static std::unique_ptr<AlignmentDataInterface> getSequenceSubset(
    const AlignmentDataInterface& sequenceSet,
    const Node& node);

  /**
   * @brief Extract the sequences corresponding to a given set of names.
   *
   * @param sequenceSet The container to look in.
   * @param names       The names of the sequences to look for.
   * @return A new site container with corresponding sequences.
   * @throw Exception if an error occured.
   */
  static std::unique_ptr<AlignmentDataInterface> getSequenceSubset(
    const AlignmentDataInterface& sequenceSet,
    const std::vector<std::string>& names);

  /**
   * @brief Compress a site container by removing duplicated sites.
   *
   * @param sequenceSet The container to look in.
   * @return A new site container with unique sites.
   * @throw Exception if an error occured.
   */
  static std::unique_ptr<AlignmentDataInterface> shrinkSiteSet(
    const AlignmentDataInterface& sequenceSet);

  /**
   * @brief Look for the occurence of each site in sequences1 in sequences2 and send the
   * position of the first occurence, or -1 if not found.
   *
   * @param sequences1 First container.
   * @param sequences2 Second container.
   * @return A vector of positions.
   */
  static Vint getIndexes(
    const AlignmentDataInterface& sequences1,
    const AlignmentDataInterface& sequences2);
};

template<class N, class E, class I>
std::unique_ptr<AlignmentDataInterface> PatternTools::getSequenceSubset(
    const AlignmentDataInterface& sequenceSet,
    const std::shared_ptr<N> node,
    const AssociationTreeGraphImplObserver<N, E, I>& tree)
{
  size_t nbSites = sequenceSet.getNumberOfSites();
  auto alphabet = sequenceSet.getAlphabet();

  try
  {
    const auto& sitecontainer = dynamic_cast<const SiteContainerInterface&>(sequenceSet);

    auto sequenceSubset = std::make_unique<VectorSiteContainer>(alphabet);

    std::vector<std::shared_ptr<N>> leaves = tree.getLeavesUnderNode(node);

    for (auto i : leaves)
    {
      if (i->hasName())
      {
        // Use sequence name as key.
        try
        {
          auto newSeq = std::make_unique<Sequence>(sitecontainer.sequence(i->getName()));
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
  catch (std::bad_cast& e)
  {}

  try
  {
    const ProbabilisticSiteContainerInterface& sitecontainer = dynamic_cast<const ProbabilisticSiteContainerInterface&>(sequenceSet);

    auto sequenceSubset = std::make_unique<ProbabilisticVectorSiteContainer>(alphabet);

    std::vector<std::shared_ptr<N>> leaves = tree.getLeavesUnderNode(node);

    for (auto i : leaves)
    {
      if (i->hasName())
      {
        // Use sequence name as key.
        try
        {
          auto newSeq = std::make_unique<ProbabilisticSequence>(sitecontainer.sequence(i->getName()));
          sequenceSubset->addSequence(newSeq->getName(), newSeq);
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
  catch (std::bad_cast& e)
  {}

  throw Exception("PatternTools::getSequenceSubset : unsupported sequence type.");
}
} // end of namespace bpp.
#endif // BPP_PHYL_PATTERNTOOLS_H
