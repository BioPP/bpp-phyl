// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/App/ApplicationTools.h>

#include "../PatternTools.h"
#include "../Tree/TreeTemplateTools.h"
#include "AbstractTreeParsimonyScore.h"

using namespace bpp;
using namespace std;

AbstractTreeParsimonyScore::AbstractTreeParsimonyScore(
  shared_ptr<TreeTemplate<Node>> tree,
  shared_ptr<const SiteContainerInterface> data,
  bool verbose,
  bool includeGaps) :
  treePtr_(std::move(tree)),
  data_(nullptr),
  alphabet_(data->getAlphabet()),
  statesMap_(nullptr),
  nbStates_(0)
{
  statesMap_ = make_shared<CanonicalStateMap>(alphabet_, includeGaps);
  nbStates_  = statesMap_->getNumberOfModelStates();
  init_(data, verbose);
}


AbstractTreeParsimonyScore::AbstractTreeParsimonyScore(
  shared_ptr<TreeTemplate<Node>> tree,
  shared_ptr<const SiteContainerInterface> data,
  std::shared_ptr<const StateMapInterface> statesMap,
  bool verbose) :
  treePtr_(std::move(tree)),
  data_(nullptr),
  alphabet_(data->getAlphabet()),
  statesMap_(statesMap),
  nbStates_(statesMap->getNumberOfModelStates())
{
  init_(data, verbose);
}

void AbstractTreeParsimonyScore::init_(shared_ptr<const SiteContainerInterface> data, bool verbose)
{
  if (treePtr_->isRooted())
  {
    if (verbose)
      ApplicationTools::displayWarning("Tree has been unrooted.");
    treePtr_->unroot();
  }
  TreeTemplateTools::deleteBranchLengths(*treePtr_->getRootNode());

  // Sequences will be in the same order than in the tree:
  shared_ptr<AlignmentDataInterface> tmp = PatternTools::getSequenceSubset(*data, *treePtr_->getRootNode());
  data_ = dynamic_pointer_cast<const SiteContainerInterface>(tmp);
  if (!data_)
    throw Exception("AbstractTreeParsimonyScore::init_ : Data must be plain alignments.");

  if (data_->getNumberOfSequences() == 1)
    throw Exception("Error, only 1 sequence!");
  if (data_->getNumberOfSequences() == 0)
    throw Exception("Error, no sequence!");
  if (data_->getAlphabet()->getSize() > 20)
    throw Exception("Error, only alphabet with size <= 20 are supported. See the source file of AbstractTreeParsimonyScore.");
}

std::vector<unsigned int> AbstractTreeParsimonyScore::getScorePerSite() const
{
  vector<unsigned int> scores(data_->getNumberOfSites());
  for (size_t i = 0; i < scores.size(); i++)
  {
    scores[i] = getScoreForSite(i);
  }
  return scores;
}
