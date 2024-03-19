// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_IO_NEXUSIOTREE_H
#define BPP_PHYL_IO_NEXUSIOTREE_H


#include "../Tree/PhyloTree.h"
#include "../Tree/TreeTemplate.h"
#include "IoTree.h"

namespace bpp
{
/**
 * @brief a simple parser for reading trees from a Nexus file.
 *
 * This reader is not supposed to be a full parser of the Nexus files,
 * but only extract the tree data. Only a basic subset of the options
 * are and will be supported.
 *
 * This format is described in the following paper:
 * Maddison D, Swofford D, and Maddison W (1997), _Syst Biol_ 46(4):590-621
 *
 * @author Julien Dutheil
 */
class NexusIOTree :
  public virtual AbstractITree,
  public virtual AbstractOTree,
  public virtual AbstractIMultiTree,
  public virtual AbstractOMultiTree,
  public AbstractIPhyloTree,
  public AbstractOPhyloTree,
  public AbstractIMultiPhyloTree,
  public AbstractOMultiPhyloTree
{
public:
  /**
   * @brief Build a new Nexus tree parser.
   */
  NexusIOTree() {}

  virtual ~NexusIOTree() {}

public:
  /**
   * @name The IOTree interface
   *
   * @{
   */
  const std::string getFormatName() const override;
  const std::string getFormatDescription() const override;
  /* @} */

  /**
   * @name The ITree interface
   *
   * @{
   */
  using AbstractITree::readTreeTemplate;

  std::unique_ptr<TreeTemplate<Node>> readTreeTemplate(std::istream& in) const override;

  using AbstractIPhyloTree::readPhyloTree;

  std::unique_ptr<PhyloTree> readPhyloTree(std::istream& in) const override;

  /** @} */

  /**
   * @name The OTree interface
   *
   * @{
   */
  using AbstractOTree::writeTree;

  void writeTree(const Tree& tree, std::ostream& out) const override
  {
    write_(tree, out);
  }

  using AbstractOPhyloTree::writePhyloTree;

  void writePhyloTree(const PhyloTree& tree, std::ostream& out) const override
  {
    write_(tree, out);
  }
  /** @} */

  /**
   * @name The IMultiTree interface
   *
   * @{
   */
  using AbstractIMultiTree::readTrees;

  void readTrees(std::istream& in, std::vector<std::unique_ptr<Tree>>& trees) const override;

  using AbstractIMultiPhyloTree::readPhyloTrees;

  void readPhyloTrees(std::istream& in, std::vector<std::unique_ptr<PhyloTree>>& trees) const override;
  /**@}*/

  /**
   * @name The OMultiTree interface
   *
   * @{
   */
  using AbstractOMultiTree::writeTrees;

  void writeTrees(const std::vector<const Tree*>& trees, std::ostream& out) const override
  {
    write_(trees, out);
  }

  using AbstractOMultiPhyloTree::writePhyloTrees;

  void writePhyloTrees(const std::vector<const PhyloTree*>& trees, std::ostream& out) const override
  {
    write_(trees, out);
  }
  /** @} */

protected:
  void write_(const Tree& tree, std::ostream& out) const;

  void write_(const PhyloTree& tree, std::ostream& out) const;

  template<class N>
  void write_(const TreeTemplate<N>& tree, std::ostream& out) const;

  void write_(const std::vector<const Tree*>& trees, std::ostream& out) const;

  void write_(const std::vector<const PhyloTree*>& trees, std::ostream& out) const;

  template<class N>
  void write_(const std::vector<TreeTemplate<N>*>& trees, std::ostream& out) const;
};
} // end of namespace bpp.
#endif // BPP_PHYL_IO_NEXUSIOTREE_H
