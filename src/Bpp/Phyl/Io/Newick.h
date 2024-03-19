// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_IO_NEWICK_H
#define BPP_PHYL_IO_NEWICK_H


#include "../Tree/PhyloTree.h"
#include "../Tree/TreeTemplate.h"
#include "IoTree.h"

namespace bpp
{
/**
 * @brief The so-called 'newick' parenthetic format.
 *
 * Branch lengths and bootstraps are supported:
 *
 * ex:
 * <code>
 * ((A:0.1, B:0.15)90:0.2, C:0.27);
 * </code>
 *
 * Code example:
 * @code
 * #include <Phyl/Newick.h>
 * #include <Phyl/Tree/Tree.h>
 *
 * Newick * newickReader = new Newick(false); //No comment allowed!
 * try {
 *  Tree * tree = newickReader->read("MyTestTree.dnd"); // Tree in file MyTestTree.dnd
 *  cout << "Tree has " << tree->getNumberOfLeaves() << " leaves." << endl;
 * } catch (Exception e) {
 *	cout << "Error when reading tree." << endl;
 * }
 * delete tree;
 * delete newickReader;
 * @endcode
 *
 * Bootstrap values are stored as node properties, as Number<double> objects and with the tag TreeTools::BOOTSTRAP.
 *
 * This is also possible to read a "tagged" tree, where additional info is provided in place of bootstrap values:
 * ((A,B)N2,(C,D)N3)N1;
 * This is achieved by calling the enableExtendedBootstrapProperty method, and providing a property name to use.
 * The additional information will be stored at each node as a property, in a String object.
 * The disableExtendedBootstrapProperty method restores the default behavior.
 */
class Newick :
  public AbstractITree,
  public AbstractOTree,
  public AbstractIMultiTree,
  public AbstractOMultiTree,
  public AbstractIPhyloTree,
  public AbstractOPhyloTree,
  public AbstractIMultiPhyloTree,
  public AbstractOMultiPhyloTree
{
protected:
  bool allowComments_;
  bool writeId_;
  bool useBootstrap_;
  std::string bootstrapPropertyName_;
  bool verbose_;

public:
  /**
   * @brief Build a new Newick reader/writer.
   *
   * Some newick format allow comments between hooks ('[' ']').
   *
   * @param allowComments Tell if comments between [] are allowed in file.
   * @param writeId       If true, nodes ids will be written in place of bootstrap values.
   * @param verbose       If some info should be displayed, such as progress bar etc.
   */
  Newick(bool allowComments = false, bool writeId = false, bool verbose = false) :
    allowComments_(allowComments),
    writeId_(writeId),
    useBootstrap_(true),
    bootstrapPropertyName_("bootstrap"),
    verbose_(verbose) {}

  virtual ~Newick() {}

public:
  void enableExtendedBootstrapProperty(const std::string& propertyName)
  {
    useBootstrap_ = false;
    bootstrapPropertyName_ = propertyName;
  }

  void disableExtendedBootstrapProperty()
  {
    useBootstrap_ = true;
    bootstrapPropertyName_ = "bootstrap";
  }

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

private:
  std::shared_ptr<PhyloNode> parenthesisToNode(
      PhyloTree& tree,
      std::shared_ptr<PhyloNode> father,
      const std::string& description,
      unsigned int& nodeCounter,
      bool bootstrap,
      const std::string& propertyName,
      bool withId,
      bool verbose) const;

public:
  std::unique_ptr<PhyloTree> parenthesisToPhyloTree(
      const std::string& description,
      bool bootstrap = false,
      const std::string& propertyName = "",
      bool withId = false,
      bool verbose = false) const;

/** @} */

  /**
   * @name The OTree interface
   *
   * @{
   */

public:
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

  void readTrees(
      std::istream& in,
      std::vector<std::unique_ptr<Tree>>& trees) const override;

  using AbstractIMultiPhyloTree::readPhyloTrees;

  void readPhyloTrees(
      std::istream& in,
      std::vector<std::unique_ptr<PhyloTree>>& trees) const override;

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

  IOTree::Element getElement(const std::string& elt) const override;

/**
 * @brief Get the Newick description of a subtree.
 *
 * @param tree The tree to convert.
 * @param node The top of the subtree to convert.
 * @param writeId Tells if node ids must be printed.
 *                This will overwrite bootstrap values if there are ones.
 *                Leaves id will be added to the leave names, separated by a '_' character.
 * @return A string in the parenthesis format.
 */

  std::string nodeToParenthesis(const PhyloTree& tree, std::shared_ptr<PhyloNode> node, bool writeId = false) const;

/* @brief Get the parenthesis description of a subtree.
 *
 * @param tree The tree
 * @param node The node defining the subtree.
 * @param bootstrap Tell is bootstrap values must be writen.
 * If so, the content of the property with name "bootstrap" will be written as bootstrap value.
 * The property should be a Number<double> object.
 * Otherwise, the content of the property with name 'propertyName' will be written.
 * In this later case, the property should be a String object.
 * @param propertyName The name of the property to use. Only used if bootstrap = false.
 * @return A string in the parenthesis format.
 */

  std::string nodeToParenthesis(const PhyloTree& tree, std::shared_ptr<PhyloNode> node, bool bootstrap, const std::string& propertyName) const;

/**
 * @brief Get the parenthesis description of a tree.
 *
 * @param tree The tree to convert.
 * @param writeId Tells if node ids must be printed.
 *                This will overwrite bootstrap values if there are ones.
 *                Leaves id will be added to the leave names, separated by a '_' character.
 * @return A string in the parenthesis format.
 */

  std::string treeToParenthesis(const PhyloTree& tree, bool writeId = false) const;

/**
 * @brief Get the parenthesis description of a tree.
 *
 * @param tree The tree to convert.
 * @param bootstrap Tell is bootstrap values must be writen.
 * If so, the content of the property with name "bootstrap" will be written as bootstrap value.
 * The property should be a Number<double> object.
 * Otherwise, the content of the property with name 'propertyName' will be written.
 * In this later case, the property should be a String object.
 * @param propertyName The name of the property to use. Only used if bootstrap = false.
 * @return A string in the parenthesis format.
 */
  std::string treeToParenthesis(const PhyloTree& tree, bool bootstrap, const std::string& propertyName) const;
};
} // end of namespace bpp.
#endif // BPP_PHYL_IO_NEWICK_H
