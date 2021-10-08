//
// File: PhyloTreeTools.h
// Authors:
//   Julien Dutheil
// Created: mercredi 31 aoÃ»t 2016, Ã  13h 13
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

#ifndef BPP_PHYL_TREE_PHYLOTREETOOLS_H
#define BPP_PHYL_TREE_PHYLOTREETOOLS_H

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/VectorTools.h>

#include "BipartitionList.h"
#include "PhyloBranch.h"
#include "PhyloNode.h"
#include "PhyloTree.h"
#include "PhyloTreeExceptions.h"

// From SeqLib:
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/DistanceMatrix.h>

#include <memory>

namespace bpp
{
/**
 * @brief Generic utilitary methods dealing with trees.
 *
 * These methods work with all Tree object.
 * However, depending on the tree implementation, they may not be the most efficient.
 *
 * @see PhyloTree
 */

class PhyloTreeTools
{
public:
  PhyloTreeTools() {}
  virtual ~PhyloTreeTools() {}

public:
  /**
   * @brief Get the height of the subtree defined by node 'node', i.e. the maximum
   * distance between leaves and the root of the subtree.
   *
   * The distance do not include the branch length of the subtree root node.
   * The height of a leaf is hence 0.
   *
   * @param tree The tree.
   * @param node The node defining the subtree.
   * @return The height of the subtree.
   */

  static double getHeight(const PhyloTree& tree, const std::shared_ptr<PhyloNode> node);

  /**
   * @brief Grafen's method to initialize branch lengths.
   *
   * Each height of the node (total distance from the leaves) is set equal to the number of
   * leaf nodes for the corresponding subtrees - 1 for inner nodes, 0 for leaves.
   *
   * If the tree already has branch lengths, they will be ignored.
   *
   * Reference:
   * Grafen A. The phylogenetic regression. Philos Trans R Soc Lond B Biol Sci. 1989; 326(1233):119-57
   *
   * @param tree The tree.
   */
  static void initBranchLengthsGrafen(PhyloTree& tree);

  /**
   * @brief Compute branch lengths using Grafen's method.
   *
   * The 'height' of each node is devided by the total height of the tree, and the ratio is raised at power 'rho'.
   * A value of rho=0 hence returns a star tree.
   *
   * Reference:
   * Grafen A. The phylogenetic regression. Philos Trans R Soc Lond B Biol Sci. 1989; 326(1233):119-57
   *
   * @param tree The tree to use.
   * @param power The rho parameter.
   * @param init Tell if the height must be initialized by calling the initBranchLengthsGrafen() method.
   *             Otherwise use branch lengths.
   */

  static void computeBranchLengthsGrafen(PhyloTree& tree, double power = 1, bool init = true);

private:
  static size_t initBranchLengthsGrafen(PhyloTree& tree, std::shared_ptr<PhyloNode> node);

  static void computeBranchLengthsGrafen(PhyloTree& tree, std::shared_ptr<PhyloNode> node, double power, double total, double& height, double& heightRaised);

public:
  /**
   * @brief Modify a tree's branch lengths to make a clock tree, by rebalancing branch lengths.
   *
   * The height of each node is set to the mean height of all son nodes.
   * This may however lead to negative branch lengths, since the mean heigth
   * may be inferior to one of the son heights, due to short branch lengths.
   * If the 'noneg' is set to yes, the mean height is checked against all son
   * heights. If it is inferior to one of the son heights, the maximum son
   * height is used instead. This results in a multifurcation.
   *
   * This method is recursive and will be applied on all sons nodes.
   *
   * @param tree The tree to use.
   * @param node The node defining the subtree.
   * @return The modified height of the node.
   */

  static double convertToClockTree(PhyloTree& tree, std::shared_ptr<PhyloNode> node);

  /**
   * @brief Modify a tree's branch lengths to make a clock tree, by rescaling subtrees.
   *
   * The height of each node is set to the mean height of all son nodes.
   * All branch lengths of the corresponding subtrees are updated proportionally.
   * This algorithm is smaller than the convertToClockTree method, but may be more accurate.
   *
   * This method is recursive and will be applied on all sons nodes.
   *
   * @param tree The tree to use.
   * @param node The node defining the subtree.
   * @return The modified height of the node.
   */

  static double convertToClockTree2(PhyloTree& tree, std::shared_ptr<PhyloNode> node);


  /**
   * @brief Determine the mid-point position of the root along the branch that already contains the root. Consequently, the topology of the rooted tree remains identical.
   *
   * This code uses two inner functions to compute the mid-point position: statFromNode_ and bestRootPosition_.
   * This code is inspired by a code performing a similar calculation in Seaview (Guindon et al., 2010, Mol. Biol. Evol. 27(2):221-4).
   *
   * @param tree The rooted tree for which the root has to be moved to its mid-point position, along the branch where it already stands.
   */
  static void constrainedMidPointRooting(PhyloTree& tree);

  /**
   * @name Some properties.
   *
   * @{
   */

private:
  struct Moments_
  {
    double N;
    double sum, squaredSum;
    Moments_() : N(0),
      sum(0),
      squaredSum(0) {}
  };

  static Moments_ statFromNode_(const PhyloTree& tree, const std::shared_ptr<PhyloNode> root);
  static double bestRootPosition_(const PhyloTree& tree, const std::shared_ptr<PhyloNode> node1, const std::shared_ptr<PhyloNode> node2, double length);


  /** @} */
};
} // end of namespace bpp.
#endif // BPP_PHYL_TREE_PHYLOTREETOOLS_H
