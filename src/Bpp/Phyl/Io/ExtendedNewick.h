// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_IO_EXTENDEDNEWICK_H
#define BPP_PHYL_IO_EXTENDEDNEWICK_H


#include "../Tree/PhyloDAG.h"
#include "IoDAG.h"

namespace bpp
{
/**
 * @brief The so-called 'ExtendedNewick' parenthetic format for
 * phylogenetic networks, where hybridization nodes are mandatory
 * tagged and appear several times.
 *
 * Branch values are supported:
 *
 * ex:
 * <code>
 * ((A:0.1, ((B:0.2, (C:0.3, (D:0.2)Y#H1:0.2):0.1):0.2, (((Y#H1:0.1, E:0.4)h, 6)f)X#H2)c)a, ((X#H2, 7)d, 8)b)r;
 * </code>
 *
 *
 * Those values can be lengths, or probabilities for hybridization
 * branches, depending on the biological meaning of the network.
 *
 * Reference : Cardona, G, ó, F, Valiente, G (2008). Extended Newick:
 * it is time for a standard representation of phylogenetic networks.
 * BMC Bioinformatics, 9:532.
 *
 */

class ExtendedNewick :
  public AbstractIPhyloDAG,
  public AbstractOPhyloDAG,
  public AbstractIMultiPhyloDAG,
  public AbstractOMultiPhyloDAG
{
protected:
  bool allowComments_;
  bool writeId_;
  bool verbose_;
  
public:
  /**
   * @brief Build a new ExtendedNewick reader/writer.
   *
   * Some extendednewick format allow comments between hooks ('[' ']').
   *
   * @param allowComments Tell if comments between [] are allowed in file.
   * @param writeId       If true, nodes ids will be written in place of bootstrap values.
   * @param verbose       If some info should be displayed, such as progress bar etc.
   */
  ExtendedNewick(bool allowComments = false, bool writeId = false, bool verbose = false) :
    allowComments_(allowComments),
    writeId_(writeId),
    verbose_(verbose)
  {}

  virtual ~ExtendedNewick() {}

public:
  /**
   * @name The IODAGraph interface
   *
   * @{
   */
  const std::string getFormatName() const override;
  const std::string getFormatDescription() const override;
  /* @} */

  /**
   * @name The IDAGraph interface
   *
   * @{
   */

  using AbstractIPhyloDAG::readPhyloDAG;

  std::unique_ptr<PhyloDAG> readPhyloDAG(std::istream& in) const override;

private:
  std::shared_ptr<PhyloNode> parenthesisToNode(
      PhyloDAG& dag,
      std::shared_ptr<PhyloNode> father,
      const std::string& description,
      unsigned int& nodeCounter,
      unsigned int& branchCounter,
      std::map<std::string, std::shared_ptr<PhyloNode> >& mapEvent,
      bool withId,
      bool verbose) const;

public:
  std::unique_ptr<PhyloDAG> parenthesisToPhyloDAG(
      const std::string& description,
      bool withId,
      bool verbose = false) const;

/** @} */

  /**
   * @name The ODAGraph interface
   *
   * @{
   */

public:
  using AbstractOPhyloDAG::writePhyloDAG;

  void writePhyloDAG(const PhyloDAG& dag, std::ostream& out) const override
  {
    write_(dag, out);
  }
  /** @} */

  /**
   * @name The IMultiDAGraph interface
   *
   * @{
   */

  using AbstractIMultiPhyloDAG::readPhyloDAGs;

  void readPhyloDAGs(
      std::istream& in,
      std::vector<std::unique_ptr<PhyloDAG>>& dags) const override;

  /**@}*/

  using AbstractOMultiPhyloDAG::writePhyloDAGs;

  void writePhyloDAGs(const std::vector<const PhyloDAG*>& dags, std::ostream& out) const override
  {
    write_(dags, out);
  }
  /** @} */

protected:
  void write_(const PhyloDAG& tree, std::ostream& out) const;

  void write_(const std::vector<const PhyloDAG*>& dags, std::ostream& out) const;

  Element getElement(const std::string& elt) const override;

/**
 * @brief Get the ExtendedNewick description of a subdag.
 *
 * @param dag The dag to convert.
 * @param node The top of the subdag to convert.
 * @param writeId Tells if node ids must be printed.
 *                Leaves id will be added to the leave names, separated by a '_' character.
 * @return A string in the parenthesis format.
 */

  std::string edgeToParenthesis(const PhyloDAG& dag, std::shared_ptr<PhyloBranch> edge, std::vector<std::shared_ptr<PhyloNode>>& writtenNodes, bool writeId = false) const;

/**
 * @brief Get the parenthesis description of a tree.
 *
 * @param dag The dag to convert.
 * @param writeId Tells if node ids must be printed.
 *                Leaves id will be added to the leave names, separated by a '_' character.
 * @return A string in the parenthesis format.
 */

  std::string dagToParenthesis(const PhyloDAG& dag, bool writeId = false) const;

};
} // end of namespace bpp.
#endif // BPP_PHYL_IO_EXTENDEDNEWICK_H
