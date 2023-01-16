//
// File: NexusIoTree.h
// Authors:
//   Julien Dutheil
// Created: 2009-05-27 19:06:00
//

/*
  Copyright or ÃÂ© or Copr. CNRS, (November 16, 2004)
  
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
