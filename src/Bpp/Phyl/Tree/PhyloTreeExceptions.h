// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_TREE_PHYLOTREEEXCEPTIONS_H
#define BPP_PHYL_TREE_PHYLOTREEEXCEPTIONS_H

#include <Bpp/Exceptions.h>
#include <Bpp/Text/TextTools.h>


// From the STL:
#include <string>
#include <memory>

namespace bpp
{
class PhyloNode;
class PhyloBranch;
class PhyloTree;

/**
 * @brief General exception thrown when something is wrong with a particular node.
 */
class PhyloNodeException :
  public Exception
{
protected:
  unsigned int nodeId_;

public:
  /**
   * @brief Build a new PhyloNodeException.
   * @param text A message to be passed to the exception hierarchy.
   * @param nodeId The id of the node that threw the exception.
   */
  PhyloNodeException(const std::string& text, unsigned int nodeId) :
    Exception("PhyloNodeException: " + text + "(id:" + TextTools::toString(nodeId) + ")"),
    nodeId_(nodeId) {}

  virtual ~PhyloNodeException() {}

public:
  /**
   * @brief Get the id of node that threw the exception.
   *
   * @return The id of the faulty node.
   */
  virtual unsigned int getNodeId() const { return nodeId_; }
};


/**
 * @brief General exception thrown when something is wrong with a particular node.
 */
class PhyloNodePException :
  public PhyloNodeException
{
private:
  const PhyloNode* node_;

public:
  /**
   * @brief Build a new PhyloNodePException.
   * @param text A message to be passed to the exception hierarchy.
   * @param tree the phyloTree owning the branch
   * @param node A const pointer toward the node that threw the exception.
   */
  PhyloNodePException(const std::string& text, const PhyloTree& tree, const std::shared_ptr<PhyloNode> node);

  /**
   * @brief Build a new PhyloNodePException.
   * @param text A message to be passed to the exception hierarchy.
   * @param node A const pointer toward the node that threw the exception.
   */
  PhyloNodePException(const std::string& text, const PhyloNode* node);

  /**
   * @brief Build a new PhyloNodePException.
   * @param text A message to be passed to the exception hierarchy.
   * @param nodeId The id of the node that threw the exception.
   */
  PhyloNodePException(const std::string& text, unsigned int nodeId) :
    PhyloNodeException(text, nodeId), node_(0) {}

  PhyloNodePException(const PhyloNodePException& nex) :
    PhyloNodeException(nex),
    node_(nex.node_)
  {}

  PhyloNodePException& operator=(const PhyloNodePException& nex)
  {
    PhyloNodeException::operator=(nex);
    node_ = nex.node_;
    return *this;
  }

  virtual ~PhyloNodePException() {}

public:
  /**
   * @brief Get the node that threw the exception.
   *
   * @return A pointer toward the faulty node.
   */
  virtual const PhyloNode* getNode() const { return node_; }
  /**
   * @brief Get the id of node that threw the exception.
   *
   * @return The id of the faulty node.
   */
  virtual unsigned int getNodeId() const { return nodeId_; }
};

/**
 * @brief General exception thrown if a property could not be found.
 */
class PhyloNodePropertyNotFoundException :
  public PhyloNodePException
{
private:
  std::string propertyName_;

public:
  /**
   * @brief Build a new PropertyNotFoundException (Node).
   * @param text A message to be passed to the exception hierarchy.
   * @param propertyName The name of the property.
   * @param tree the phyloTree owning the node
   * @param node A const pointer toward the node that threw the exception.
   */
  PhyloNodePropertyNotFoundException(const std::string& text, const std::string& propertyName, const PhyloTree& tree, const std::shared_ptr<PhyloNode> node) :
    PhyloNodePException("Property not found: " + propertyName + ". " + text, tree, node),
    propertyName_(propertyName) {}

  /**
   * @brief Build a new PropertyNotFoundException (Node).
   * @param text A message to be passed to the exception hierarchy.
   * @param propertyName The name of the property.
   * @param node A const pointer toward the node that threw the exception.
   */
  PhyloNodePropertyNotFoundException(const std::string& text, const std::string& propertyName, const PhyloNode* node) :
    PhyloNodePException("Property not found: " + propertyName + ". " + text, node),
    propertyName_(propertyName) {}

  /**
   * @brief Build a new PropertyNotFoundException.
   *
   * @param text A message to be passed to the exception hierarchy.
   * @param propertyName The name of the property.
   * @param nodeId The id of the node that threw the exception.
   */
  PhyloNodePropertyNotFoundException(const std::string& text, const std::string& propertyName, unsigned int nodeId) :
    PhyloNodePException("Property not found: " + propertyName + ". " + text, nodeId),
    propertyName_(propertyName) {}

  virtual ~PhyloNodePropertyNotFoundException() {}

public:
  /**
   * @brief Get the name of the property that could not be found.
   *
   * @return The name of the missing property.
   */
  virtual const std::string& getPropertyName() const { return propertyName_; }
};

/**
 * @brief Exception thrown when something is wrong with a particular node.
 */
class PhyloNodeNotFoundException :
  public Exception
{
private:
  std::string id_;

public:
  /**
   * @brief Build a new PhyloNodeNotFoundException.
   *
   * @param text A message to be passed to the exception hierarchy.
   * @param id   A string describing the node.
   */
  PhyloNodeNotFoundException(const std::string& text, const std::string& id);

  /**
   * @brief Build a new PhyloNodeNotFoundException.
   *
   * @param text A message to be passed to the exception hierarchy.
   * @param id   A node identifier.
   */
  PhyloNodeNotFoundException(const std::string& text, unsigned int id);

  virtual ~PhyloNodeNotFoundException() {}

public:
  /**
   * @brief Get the node id that threw the exception.
   *
   * @return The id of the node.
   */
  virtual std::string getId() const { return id_; }
};

/**
 * @brief General exception thrown when something is wrong with a particular branch.
 */
class PhyloBranchException :
  public Exception
{
protected:
  unsigned int branchId_;

public:
  /**
   * @brief Build a new PhyloBranchPException.
   * @param text A message to be passed to the exception hierarchy.
   * @param branchId The id of the branch that threw the exception.
   */
  PhyloBranchException(const std::string& text, unsigned int branchId) :
    Exception("PhyloBranchException: " + text + "(id:" + TextTools::toString(branchId) + ")"),
    branchId_(branchId) {}

  virtual ~PhyloBranchException() {}

public:
  /**
   * @brief Get the id of branch that threw the exception.
   *
   * @return The id of the faulty branch.
   */
  virtual unsigned int getBranchId() const { return branchId_; }
};


/**
 * @brief General exception thrown when something is wrong with a particular branch.
 */
class PhyloBranchPException :
  public PhyloBranchException
{
private:
  const PhyloBranch* branch_;

public:
  /**
   * @brief Build a new PhyloBranchPException.
   * @param text A message to be passed to the exception hierarchy.
   * @param tree the phyloTree owning the branch
   * @param branch A const pointer toward the branch that threw the exception.
   */
  PhyloBranchPException(const std::string& text, const PhyloTree& tree, const std::shared_ptr<PhyloBranch> branch);

  /**
   * @brief Build a new PhyloBranchPException.
   * @param text A message to be passed to the exception hierarchy.
   * @param branch A const pointer toward the branch that threw the exception.
   */
  PhyloBranchPException(const std::string& text, const PhyloBranch* branch);

  /**
   * @brief Build a new PhyloBranchPException.
   * @param text A message to be passed to the exception hierarchy.
   * @param branchId The id of the branch that threw the exception.
   */
  PhyloBranchPException(const std::string& text, unsigned int branchId) :
    PhyloBranchException(text, branchId), branch_(0) {}

  PhyloBranchPException(const PhyloBranchPException& nex) :
    PhyloBranchException(nex),
    branch_(nex.branch_)
  {}

  PhyloBranchPException& operator=(const PhyloBranchPException& nex)
  {
    PhyloBranchException::operator=(nex);
    branch_ = nex.branch_;
    return *this;
  }

  virtual ~PhyloBranchPException() {}

public:
  /**
   * @brief Get the branch that threw the exception.
   *
   * @return A pointer toward the faulty branch.
   */
  virtual const PhyloBranch* getBranch() const { return branch_; }
  /**
   * @brief Get the id of branch that threw the exception.
   *
   * @return The id of the faulty branch.
   */
  virtual unsigned int getBranchId() const { return branchId_; }
};

/**
 * @brief General exception thrown if a property could not be found.
 */
class PhyloBranchPropertyNotFoundException :
  public PhyloBranchPException
{
private:
  std::string propertyName_;

public:
  /**
   * @brief Build a new PropertyNotFoundException (Branch).
   * @param text A message to be passed to the exception hierarchy.
   * @param propertyName The name of the property.
   * @param tree the phyloTree owning the branch
   * @param branch A const pointer toward the branch that threw the exception.
   */

  PhyloBranchPropertyNotFoundException(const std::string& text, const std::string& propertyName, const PhyloTree& tree, const std::shared_ptr<PhyloBranch> branch) :
    PhyloBranchPException("Property not found: " + propertyName + ". " + text, tree, branch),
    propertyName_(propertyName) {}

  /**
   * @brief Build a new PropertyNotFoundException (Branch).
   * @param text A message to be passed to the exception hierarchy.
   * @param propertyName The name of the property.
   * @param branch A const pointer toward the branch that threw the exception.
   */

  PhyloBranchPropertyNotFoundException(const std::string& text, const std::string& propertyName, const PhyloBranch* branch) :
    PhyloBranchPException("Property not found: " + propertyName + ". " + text, branch),
    propertyName_(propertyName) {}

  /**
   * @brief Build a new PropertyNotFoundException.
   *
   * @param text A message to be passed to the exception hierarchy.
   * @param propertyName The name of the property.
   * @param branchId The id of the branch that threw the exception.
   */
  PhyloBranchPropertyNotFoundException(const std::string& text, const std::string& propertyName, unsigned int branchId) :
    PhyloBranchPException("Property not found: " + propertyName + ". " + text, branchId),
    propertyName_(propertyName) {}

  virtual ~PhyloBranchPropertyNotFoundException() {}

public:
  /**
   * @brief Get the name of the property that could not be found.
   *
   * @return The name of the missing property.
   */
  virtual const std::string& getPropertyName() const { return propertyName_; }
};

/**
 * @brief Exception thrown when something is wrong with a particular branch.
 */
class PhyloBranchNotFoundException :
  public Exception
{
private:
  std::string id_;

public:
  /**
   * @brief Build a new PhyloBranchNotFoundException.
   *
   * @param text A message to be passed to the exception hierarchy.
   * @param id   A string describing the branch.
   */
  PhyloBranchNotFoundException(const std::string& text, const std::string& id);

  /**
   * @brief Build a new PhyloBranchNotFoundException.
   *
   * @param text A message to be passed to the exception hierarchy.
   * @param id   A branch identifier.
   */
  PhyloBranchNotFoundException(const std::string& text, unsigned int id);

  virtual ~PhyloBranchNotFoundException() {}

public:
  /**
   * @brief Get the branch id that threw the exception.
   *
   * @return The id of the branch.
   */
  virtual std::string getId() const { return id_; }
};

/**
 * @brief General exception thrown when something wrong happened in a tree.
 */
class PhyloTreeException :
  public Exception
{
private:
  const PhyloTree* tree_;

public:
  /**
   * @brief Build a new PhyloTreeException.
   *
   * @param text A message to be passed to the exception hierarchy.
   * @param tree A const pointer toward the tree that threw the exception.
   */
  PhyloTreeException(const std::string& text, const PhyloTree* tree = 0);

  PhyloTreeException(const PhyloTreeException& tex) :
    Exception(tex),
    tree_(tex.tree_)
  {}

  PhyloTreeException& operator=(const PhyloTreeException& tex)
  {
    Exception::operator=(tex);
    tree_ = tex.tree_;
    return *this;
  }

  virtual ~PhyloTreeException() {}

public:
  /**
   * @brief Get the tree that threw the exception.
   *
   * @return The faulty tree
   */
  virtual const PhyloTree* getTree() const { return tree_; }
};

/**
 * @brief Exception thrown when a tree is expected to be rooted.
 */
class UnrootedPhyloTreeException :
  public PhyloTreeException
{
public:
  /**
   * @brief Build a new UnrootedPhyloTreeException.
   *
   * @param text A message to be passed to the exception hierarchy.
   * @param tree A const pointer toward the tree that threw the exception.
   */
  UnrootedPhyloTreeException(const std::string& text, const PhyloTree* tree = 0);

  virtual ~UnrootedPhyloTreeException() {}
};
} // end of namespace bpp.
#endif // BPP_PHYL_TREE_PHYLOTREEEXCEPTIONS_H
