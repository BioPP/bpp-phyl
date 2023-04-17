//
// File: Nhx.h
// Authors:
//   Bastien Boussau
// Created: Tu Oct 19 10:24:03 2010
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

#ifndef BPP_PHYL_IO_NHX_H
#define BPP_PHYL_IO_NHX_H


#include "../Tree/PhyloTree.h"
#include "../Tree/TreeTemplate.h"
#include "IoTree.h"

// From the STL:
#include <set>

namespace bpp
{
/**
 * @brief The so-called 'Nhx - New Hampshire eXtended' parenthetic format.
 *
 * See http://www.phylosoft.org/NHX/ for details.
 *
 * Branch lengths and node annotations are supported:
 *
 * ex:
 * <code>
 * ((A:0.1[&&NHX:S=human], B:0.15[&&NHX:S=chimp]):0.2[&&NHX:B=90:D=N:S=primates], C:0.27[&&NHX:S=mouse]);
 * </code>
 *
 * Where e.g. "S=human" means that the sequence A comes from species "human", "B=90" stands for a support value,
 * and "D=N" means that there was no duplication at the node. Other tags are allowed, see http://www.phylosoft.org/NHX/.
 * By default, comments between "[" and "]" are removed, unless the opening bracket is followed by "&&NHX".
 *
 * Code example:
 * @code
 * #include <Phyl/Nhx.h>
 * #include <Phyl/Tree/Tree.h>
 *
 * Nhx * NhxReader = new Nhx();
 * try {
 *   Tree * tree = nhxReader->read("MyTestTree.dnd"); // Tree in file MyTestTree.dnd
 *   cout << "Tree has " << tree->getNumberOfLeaves() << " leaves." << endl;
 * } catch (Exception e) {
 *  cout << "Error when reading tree." << endl;
 * }
 * delete tree;
 * delete nhxReader;
 * @endcode
 *
 * All node annotations are stored as node properties, with type bppString for all properties except for support values, where a Number is used.
 *
 */
class Nhx :
  public AbstractITree,
  public AbstractOTree,
  public AbstractIMultiTree,
  public AbstractOMultiTree,
  public AbstractIPhyloTree,
  public AbstractOPhyloTree,
  public AbstractIMultiPhyloTree,
  public AbstractOMultiPhyloTree
{
public:
  struct Property
  {
public:
    /**
     * @brief The name of the property, which will be used in parsed trees.
     */
    std::string name;
    /**
     * @brief The tag of the property, as it will be found in the tree file.
     */
    std::string tag;
    /**
     * @brief Tells if the property is a branch property instead of a node property.
     */
    bool onBranch;
    /**
     * @brief The type of the property. 0 is string, 1 is integer, 2 is double, 3 is boolean.
     */
    short type;

public:
    Property(const std::string& pptName, const std::string& pptTag, bool pptOnBranch = false, short pptType = 0) :
      name(pptName), tag(pptTag), onBranch(pptOnBranch), type(pptType) {}

    bool operator<(const Property& ppt) const
    {
      return name < ppt.name;
    }
  };

private:
  std::set<Property> supportedProperties_;
  bool useTagsAsPropertyNames_;
  mutable bool hasIds_;

public:
  /**
   * @brief Build a new Nhx reader/writer.
   *
   * Comments between hooks ('[' ']') are ignored.
   *
   * @param useTagsAsPptNames Tells if the NHX tag should be used as a property name in the parsed tree.
   */
  Nhx(bool useTagsAsPptNames = true);
  virtual ~Nhx() {}

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

  std::unique_ptr<TreeTemplate<Node>> parenthesisToTree(const std::string& description) const;

  std::unique_ptr<PhyloTree> parenthesisToPhyloTree(const std::string& description) const;

  std::string treeToParenthesis(const TreeTemplate<Node>& tree) const;

  std::string treeToParenthesis(const PhyloTree& tree) const;

  void registerProperty(const Property& property)
  {
    supportedProperties_.insert(property);
  }

  /**
   * @brief Convert property names from tag to names.
   *
   * If a tree has been parsed using useTagsAsPropertyNames=true,
   * this method allows to convert the tree as is it was parsed using
   * the option set to false.
   *
   * @param node The root node of the subtree to convert.
   */
  void changeTagsToNames(Node& node) const;
  void changeTagsToNames(PhyloTree& tree, std::shared_ptr<PhyloNode> node) const;

  /**
   * @brief Convert property names from names to tags.
   *
   * If a tree has been parsed using useTagsAsPropertyNames=false,
   * this method allows to convert the tree as is it was parsed using
   * the option set to true.
   *
   * @param node The root node of the subtree to convert.
   */
  void changeNamesToTags(Node& node) const;
  void changeNamesToTags(PhyloTree& tree, std::shared_ptr<PhyloNode> node) const;

  void useTagsAsPropertyNames(bool yn) { useTagsAsPropertyNames_ = yn; }
  bool useTagsAsPropertyNames() const { return useTagsAsPropertyNames_; }

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

private:
  Node* parenthesisToNode(const std::string& description) const;

  std::shared_ptr<PhyloNode> parenthesisToNode(PhyloTree& tree, std::shared_ptr<PhyloNode> father, const std::string& description) const;

public:
  std::string propertiesToParenthesis(const Node& node) const;

  std::string propertiesToParenthesis(const PhyloTree& tree, const std::shared_ptr<PhyloNode> node) const;

protected:
  std::string nodeToParenthesis(const Node& node) const;

  std::string nodeToParenthesis(const PhyloTree& tree, const std::shared_ptr<PhyloNode> node) const;

  bool setNodeProperties(Node& node, const std::string properties) const;

  bool setNodeProperties(PhyloTree& tree, std::shared_ptr<PhyloNode> node, const std::string properties) const;

protected:
  static std::string propertyToString_(const Clonable* pptObject, short type);
  static Clonable* stringToProperty_(const std::string& pptDesc, short type);

  /**
   * @brief check and fill all nodes ids.
   */
  void checkNodesId_(PhyloTree& tree) const;
};
} // end of namespace bpp.
#endif // BPP_PHYL_IO_NHX_H
