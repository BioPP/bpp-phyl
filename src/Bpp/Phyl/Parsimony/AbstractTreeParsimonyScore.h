//
// File: AbstractTreeParsimonyScore.h
// Authors:
//   Julien Dutheil
// Created: 2005-07-28 17:25:00
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

#ifndef BPP_PHYL_PARSIMONY_ABSTRACTTREEPARSIMONYSCORE_H
#define BPP_PHYL_PARSIMONY_ABSTRACTTREEPARSIMONYSCORE_H


#include "../Model/StateMap.h"
#include "TreeParsimonyScore.h"

// From bpp-seq:
#include <Bpp/Seq/Container/SiteContainer.h>

namespace bpp
{
/**
 * @brief Partial implementation of the TreeParsimonyScore interface.
 */
class AbstractTreeParsimonyScore :
  public virtual TreeParsimonyScoreInterface
{
private:
  std::shared_ptr<TreeTemplate<Node>> treePtr_;
  std::shared_ptr<const SiteContainerInterface> data_;
  std::shared_ptr<const Alphabet> alphabet_;
  std::shared_ptr<const StateMapInterface> statesMap_;
  size_t nbStates_;

public:
  AbstractTreeParsimonyScore(
    std::shared_ptr<TreeTemplate<Node>> tree,
    std::shared_ptr<const SiteContainerInterface> data,
    bool verbose,
    bool includeGaps);

  AbstractTreeParsimonyScore(
    std::shared_ptr<TreeTemplate<Node>> tree,
    std::shared_ptr<const SiteContainerInterface> data,
    std::shared_ptr<const StateMapInterface> statesMap,
    bool verbose);

  virtual ~AbstractTreeParsimonyScore() {}

private:
  void init_(std::shared_ptr<const SiteContainerInterface> data, bool verbose);

public:
  const Tree& tree() const override { return *treePtr_; }
  virtual const TreeTemplate<Node>& treeTemplate() const { return *treePtr_; }
  virtual std::shared_ptr<const TreeTemplate<Node>> getTreeTemplate() const { return treePtr_; }
  std::vector<unsigned int> getScorePerSite() const override;
  const StateMapInterface& stateMap() const override { return *statesMap_; }
  std::shared_ptr<const StateMapInterface> getStateMap() const override { return statesMap_; }

protected:
  virtual Tree& tree_() { return *treePtr_; }
  virtual TreeTemplate<Node>& treeTemplate_() { return *treePtr_; }
  virtual std::shared_ptr<TreeTemplate<Node>> getTreeTemplate_() { return treePtr_; }
};
} // end of namespace bpp.
#endif // BPP_PHYL_PARSIMONY_ABSTRACTTREEPARSIMONYSCORE_H
