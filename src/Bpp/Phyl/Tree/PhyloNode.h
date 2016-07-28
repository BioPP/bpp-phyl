//
// File: PhyloNode.h
// Created by: Thomas Bigot
// Created on: Thu Mar 13 12:03:18 2003
//

/*
 *  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)
 * 
 *  This software is a computer program whose purpose is to provide classes
 *  for phylogenetic data analysis.
 * 
 *  This software is governed by the CeCILL  license under French law and
 *  abiding by the rules of distribution of free software.  You can  use, 
 *  modify and/ or redistribute the software under the terms of the CeCILL
 *  license as circulated by CEA, CNRS and INRIA at the following URL
 *  "http://www.cecill.info". 
 * 
 *  As a counterpart to the access to the source code and  rights to copy,
 *  modify and redistribute granted by the license, users are provided only
 *  with a limited warranty  and the software's author,  the holder of the
 *  economic rights,  and the successive licensors  have only  limited
 *  liability. 
 * 
 *  In this respect, the user's attention is drawn to the risks associated
 *  with loading,  using,  modifying and/or developing or reproducing the
 *  software by the user in light of its specific status of free software,
 *  that may mean  that it is complicated to manipulate,  and  that  also
 *  therefore means  that it is reserved for developers  and  experienced
 *  professionals having in-depth computer knowledge. Users are therefore
 *  encouraged to load and test the software's suitability as regards their
 *  requirements in conditions enabling the security of their systems and/or 
 *  data to be ensured and,  more generally, to use and operate it in the 
 *  same conditions as regards security. 
 * 
 *  The fact that you are presently reading this means that you have had
 *  knowledge of the CeCILL license and that you accept its terms.
 */

#ifndef _PHYLONODE_H_
#define _PHYLONODE_H_

#include "PhyloTreeExceptions.h"

#include "PhyloBranch.h"
#include "PhyloTree.h"


namespace bpp
{
  
  class PhyloNode
  {
  private:
    // the tree this node is in
    PhyloTree *phyloTree_;
    
    // a name, if specified
    std::string name_;
    
    // Node properties
    mutable std::map<std::string, Clonable*> properties_;
    
    /**
     * @brief Build a new void Node object.
     * @param index the wanted index. If zero: the next available index
     */
    unsigned int setIndex_(unsigned int index = 0)
    {
      return phyloTree_->setNodeIndex(this,index);
    }
    
    
    
  public:
    /**
     * @brief Build a new void Node object.
     * @param tree 
     */
    PhyloNode(PhyloTree *tree) :
    phyloTree_(tree),
    name_(""),
    properties_()
    {
    }
    
    /**
     * @brief Build a new Node with specified id.
     */
    PhyloNode(PhyloTree *tree, unsigned int index) :
    phyloTree_(tree),
    name_(""),
    properties_()
    {
      tree->setNodeIndex(this,index);
    }
    
    /**
     * @brief Build a new Node with specified name.
     */
    PhyloNode(PhyloTree *tree, const std::string& name):
    phyloTree_(tree),
    name_(name),
    properties_()
    {}
    
    /**
     * @brief Build a new Node with specified id and name.
     */
    PhyloNode(PhyloTree *tree,unsigned int index, const std::string& name) :
    phyloTree_(tree),
    name_(name),
    properties_()
    {
      tree->setNodeIndex(this,index);
    }
    
    /**
     * @brief Copy constructor.
     *
     * @warning This constructor copies all fields, excepted father and son node pointers.
     *
     * @param node The node to copy.
     */
    PhyloNode(const PhyloNode& node);
    
    /**
     * @brief Assignation operator.
     *
     * @warning This operator copies all fields, excepted father and son node pointers.
     *
     * @param node the node to copy.
     * @return A reference toward this node.
     */
    PhyloNode& operator=(const PhyloNode& node);
    
    PhyloNode* clone() const { return new PhyloNode(*this); }
    
    virtual ~PhyloNode()
    {
      for (std::map<std::string, Clonable*>::iterator i = properties_.begin(); i != properties_.end(); i++)
      {
        delete i->second;
      }
    }
    
  public:
    /**
     * @name Identity
     *
     * @{
     */
    
    /**
     * @brief Get the node's id.
     *
     * @return The identity tag of this node.
     */

    int getId() const { return phyloTree_->getNodeIndex(this); }
    
    /**
     * @brief Set this node's id.
     *
     * @param id The new identity tag.
     */

    void setId(int id) { phyloTree_->setNodeIndex(this,id); }
    
    std::vector<int> getSonsId() const
    {
      std::vector<PhyloNode*> sons = getSons();
      std::vector<int> sonsId(sons.size());
      for (size_t i = 0; i < sons.size(); i++)
      {
        sonsId[i] = phyloTree_->getNodeIndex(sons[i]);
      }
      return sonsId;
    }
    
    /**
     * @brief Get the node leading to this branch.
     */

    PhyloBranch* getBranch() const
    {
        return phyloTree_->getBranchToFather(this);
    }
    
    
    /** @} */
    
    /**
     * @name Name:
     *
     * @{
     */
    
    /**
     * @brief Get the name associated to this node, if there is one,
     * otherwise throw a NodeException.
     *
     * @return The name associated to this node.
     */

    std::string getName() const 
    {
      if (!hasName()) throw PhyloNodePException("Node::getName: no name associated to this node.", this);
      return name_;
    }
    
    /**
     * @brief Give a name or update the name associated to the node.
     *
     * @param name The name to give to the node.
     */

    void setName(const std::string& name)
    {
      name_ = name;
    }
    
    /**
     * @brief Delete the name associated to this node (do nothing if there is no name).
     */

    void deleteName()
    {
      name_ = "";
    }
    
    /**
     * @brief Tell is this node has a name.
     *
     * @return True if name != 0.
     */

    bool hasName() const { return name_ != ""; }
    
    /** @} */
    
    /**
     * @name Distances:
     *
     * @{
     */
    
    /**
     * @brief Get the distance to the father node is there is one,
     * otherwise throw a NodeException.
     *
     * @return The distance to the father node.
     */

    double getDistanceToFather() const
    {
      if (!phyloTree_->isRooted())
        throw PhyloNodePException("PhyloNode::getDistanceToFather: unRooted tree has no father.", this);
      return getBranch()->getLength();
    }
    
    /**
     * @brief Set or update the distance toward the father node.
     *
     * Warning: a distance to the father node may be set even if no father node is specified.
     * This is used by several tree reconstruction methods.
     * It may also be useful for manipulating subtrees.
     *
     * @param distance The new distance to the father node.
     */

    void setDistanceToFather(double distance)
    {
      phyloTree_->getBranchToFather(this)->setLength(distance);
    }
    
    /**
     * @brief Delete the distance to the father node.
     */

    void deleteDistanceToFather()
    {
      phyloTree_->getBranchToFather(this)->deleteLength();
    }
    
    /**
     * @brief Tell is this node has a distance to the father.
     *
     * @return True if distanceToFather != 0.
     */

    bool hasDistanceToFather() const
    {
      return phyloTree_->getBranchToFather(this)->hasLength();
    }
    
    /** @} */
    
    /**
     * @name Father:
     *
     * @{
     */
    
    /**
     * @brief Get the father of this node is there is one.
     *
     * @return A pointer toward the father node, 0 if there is not.
     */

    const PhyloNode* getFather() const { return phyloTree_->getFather(this); }
    
    /**
     * @brief Get the father of this node is there is one.
     *
     * @return A pointer toward the father node, 0 if there is not.
     */

    PhyloNode* getFather() { return phyloTree_->getFather(this); }
    
    int getFatherId() const { return phyloTree_->getNodeIndex(phyloTree_->getFather(this)); }
    
    /**
     * @brief Set the father node of this node.
     *
     * @param node The father node.
     */

    void setFather(PhyloNode* node) 
    {
      if (!node)
        throw NullPointerException("Node::setFather(). Empty node given as input.");
      phyloTree_->setFather(this,node);
    }
    

    /**
     * @brief Tell if this node has a father node.
     */

    bool hasFather() const { return phyloTree_->hasFather(this); }
    
    /** @} */
    
    /**
     * @name Sons:
     *
     * @{
     */

    size_t getNumberOfSons() const { return phyloTree_->getNumberOfSons(this); }
    
    std::vector<PhyloNode*> getSons() const
    {
      return phyloTree_->getSons(this);
    }
    
    void addSon(PhyloNode* node) 
    {
      if (!node)
        throw NullPointerException("Node::addSon(). Empty node given as input.");
      phyloTree_->addSon(this,node);
    }
    
    void removeSon(PhyloNode* node) 
    {
      if (!node)
        throw NullPointerException("Node::removeSon(). Empty node given as input.");
      phyloTree_->addSon(this,node);
    }
    
    void removeSons()
    {
      phyloTree_->removeSons(this);
    }
            
    /** @} */
    
    // These functions must not be declared as virtual!!
    
    std::vector<const PhyloNode*> getNeighbors() const;
    
    std::vector<PhyloNode*> getNeighbors();
    
    size_t degree() const { return phyloTree_->getNeighbors(const_cast<PhyloNode*>(this)).size();} //FIXME: to replace by a GetDegree function.
    
    
    /**
     * @name Node properties:
     *
     * @{
     */
    
    /**
     * @brief Set/add a node property.
     *
     * If no property with the same name is found, the new property will be added to the list.
     * Conversely, the property will be deleted and replaced by the new one.
     * If you want to keep a copy of the old property, consider using the removeProperty function before.
     *
     * @param name The name of the property to set.
     * @param property The property object (will be cloned).
     */

    void setProperty(const std::string& name, const Clonable& property)
    {
      if (hasProperty(name))
        delete properties_[name];
      properties_[name] = property.clone();
    }
    
    Clonable* getProperty(const std::string& name) 
    {
      if (hasProperty(name))
        return properties_[name];
      else
        throw PhyloNodePropertyNotFoundException("", name, this);
    }
    
    const Clonable* getProperty(const std::string& name) const 
    {
      if (hasProperty(name))
        return const_cast<const Clonable*>(properties_[name]);
      else
        throw PhyloNodePropertyNotFoundException("", name, this);
    }
    
    Clonable* removeProperty(const std::string& name) 
    {
      if (hasProperty(name))
      {
        Clonable* removed = properties_[name];
        properties_.erase(name);
        return removed;
      }
      else
        throw PhyloNodePropertyNotFoundException("", name, this);
    }
    
    void deleteProperty(const std::string& name) 
    {
      if (hasProperty(name))
      {
        delete properties_[name];
        properties_.erase(name);
      }
      else
        throw PhyloNodePropertyNotFoundException("", name, this);
    }
    
    /**
     * @brief Remove all node properties.
     *
     * Attached objects will not be deleted.
     */
    void removeProperties()
    {
      properties_.clear();
    }
    
    /**
     * @brief Delete all node properties.
     */
    void deleteProperties()
    {
      for (std::map<std::string, Clonable*>::iterator i = properties_.begin(); i != properties_.end(); i++)
      {
        delete i->second;
      }
      properties_.clear();
    }
    
    bool hasProperty(const std::string& name) const { return properties_.find(name) != properties_.end(); }
    
    std::vector<std::string> getPropertyNames() const { return MapTools::getKeys(properties_); }
    /** @} */

    /**
     * The following is for backward compatibility
     *
     * @{
     */
    
    /**
     * @name Node properties:
     *
     * @{
     */
    
    void setNodeProperty(const std::string& name, const Clonable& property)
    {
      setProperty(name, property);
    }
    
    Clonable* getNodeProperty(const std::string& name) 
    {
      return getProperty(name);
    }
    
    const Clonable* getNodeProperty(const std::string& name) const 
    {
      return getProperty(name);
    }
    
    Clonable* removeNodeProperty(const std::string& name) 
    {
      return removeProperty(name);
    }
    
    void deleteNodeProperty(const std::string& name) 
    {
      deleteProperty(name);
    }
    
    /**
     * @brief Remove all node properties.
     *
     * Attached objects will not be deleted.
     */
    void removeNodeProperties()
    {
      removeProperties();
    }
    
    /**
     * @brief Delete all node properties.
     */
    void deleteNodeProperties()
    {
      deleteProperties();
    }
    
    bool hasNodeProperty(const std::string& name) const { return hasProperty(name); }
    
    std::vector<std::string> getNodePropertyNames() const { return getPropertyNames(); }
    
    /** @} */
    
    /**
     * @name Branch properties:
     *
     * @{
     */
    
    /**
     * @brief Set/add a branch property.
     *
     * If no property with the same name is found, the new property will be added to the list.
     * Conversely, the property will be deleted and replaced by the new one.
     * If you want to keep a copy of the old property, consider using the removeBranchProperty function before.
     *
     * @param name The name of the property to set.
     * @param property The property object (will be cloned).
     */
    void setBranchProperty(const std::string& name, const Clonable& property)
    {
      getBranch()->setProperty(name, property);
    }
    
    Clonable* getBranchProperty(const std::string& name) 
    {
      return getBranch()->getProperty(name);
    }
    
    const Clonable* getBranchProperty(const std::string& name) const 
    {
      return getBranch()->getProperty(name);
    }
    
    Clonable* removeBranchProperty(const std::string& name) 
    {
      return getBranch()->removeProperty(name);
    }
    
    void deleteBranchProperty(const std::string& name) 
    {
      getBranch()->deleteProperty(name);
    }
    
    /**
     * @brief Remove all branch properties.
     *
     * Attached objects will not be deleted.
     */
    void removeBranchProperties()
    {
      getBranch()->removeProperties();
    }
    
    /**
     * @brief Delete all branch properties.
     */
    void deleteBranchProperties()
    {
      getBranch()->deleteProperties();
    }
    
    bool hasBranchProperty(const std::string& name) const { return getBranch()->hasProperty(name); }
    
    std::vector<std::string> getBranchPropertyNames() const { return getBranch()->getPropertyNames(); }
    
    bool hasBootstrapValue() const;
    
    double getBootstrapValue() const;
    
    /** @} */

    /** @} */
    
    /** Equality operator: */
    
    bool operator==(const PhyloNode& node) const { return phyloTree_->getNodeGraphid(this) == phyloTree_->getNodeGraphid(&node); }
    
    // Tests:
    
    bool isLeaf() const { return (phyloTree_->getOutgoingNeighbors(this).size() == 0); }
    
    
  }; //end of class node 
  
  
} //end of namespace bpp.

#else
namespace bpp{ class PhyloNode; }

#endif  //_PHYLONODE_H_
