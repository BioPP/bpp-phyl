//
// File: Tree.h
// Created by: Julien Dutheil
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
    PhyloNode(PhyloTree *tree,int id, const std::string& name) :
    phyloTree_(tree),
    name_(name),
    properties_()
    {}
    
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
    virtual int getId() const { return phyloTree_->getNodeIndex(this); }
    
    /**
     * @brief Set this node's id.
     *
     * @param id The new identity tag.
     */
    virtual void setId(int id) { phyloTree_->setNodeIndex(this,id); }
    
    virtual std::vector<int> getSonsId() const
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
    virtual std::string getName() const throw (PhyloNodePException)
    {
      if (!hasName()) throw PhyloNodePException("Node::getName: no name associated to this node.", this);
      return name_;
    }
    
    /**
     * @brief Give a name or update the name associated to the node.
     *
     * @param name The name to give to the node.
     */
    virtual void setName(const std::string& name)
    {
      name_ = name;
    }
    
    /**
     * @brief Delete the name associated to this node (do nothing if there is no name).
     */
    virtual void deleteName()
    {
      name_ = "";
    }
    
    /**
     * @brief Tell is this node has a name.
     *
     * @return True if name != 0.
     */
    virtual bool hasName() const { return name_ != ""; }
    
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
    virtual double getDistanceToFather() const
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
    virtual void setDistanceToFather(double distance)
    {
      phyloTree_->getBranchToFather(this)->setLength(distance);
    }
    
    /**
     * @brief Delete the distance to the father node.
     */
    virtual void deleteDistanceToFather()
    {
      phyloTree_->getBranchToFather(this)->deleteLength();
    }
    
    /**
     * @brief Tell is this node has a distance to the father.
     *
     * @return True if distanceToFather != 0.
     */
    virtual bool hasDistanceToFather() const
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
    virtual const PhyloNode* getFather() const { return phyloTree_->getFather(this); }
    
    /**
     * @brief Get the father of this node is there is one.
     *
     * @return A pointer toward the father node, 0 if there is not.
     */
    virtual PhyloNode* getFather() { return phyloTree_->getFather(this); }
    
    virtual int getFatherId() const { return phyloTree_->getNodeIndex(phyloTree_->getFather(this)); }
    
    /**
     * @brief Set the father node of this node.
     *
     * @param node The father node.
     */
    virtual void setFather(PhyloNode* node) throw (NullPointerException)
    {
      if (!node)
        throw NullPointerException("Node::setFather(). Empty node given as input.");
      phyloTree_->setFather(this,node);
    }
    

    /**
     * @brief Tell if this node has a father node.
     */
    virtual bool hasFather() const { return phyloTree_->hasFather(this); }
    
    /** @} */
    
    /**
     * @name Sons:
     *
     * @{
     */
    virtual size_t getNumberOfSons() const { return phyloTree_->getNumberOfSons(this); }
    
    virtual std::vector<PhyloNode*> getSons() const
    {
      return phyloTree_->getSons(this);
    }
    
    virtual void addSon(PhyloNode* node) throw (NullPointerException, PhyloNodePException)
    {
      if (!node)
        throw NullPointerException("Node::addSon(). Empty node given as input.");
      phyloTree_->addSon(this,node);
    }
    
    virtual void removeSon(PhyloNode* node) throw (PhyloNodeNotFoundException, NullPointerException)
    {
      if (!node)
        throw NullPointerException("Node::removeSon(). Empty node given as input.");
      phyloTree_->addSon(this,node);
    }
    
    virtual void removeSons()
    {
      phyloTree_->removeSons(this);
    }
            
    /** @} */
    
    // These functions must not be declared as virtual!!
    
    std::vector<const PhyloNode*> getNeighbors() const;
    
    std::vector<PhyloNode*> getNeighbors();
    
    virtual size_t degree() const { return phyloTree_->getNeighbors(const_cast<PhyloNode*>(this)).size();} //FIXME: to replace by a GetDegree function.
    
    
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
     * If you want to keep a copy of the old property, consider using the removeNodeProperty function before.
     *
     * @param name The name of the property to set.
     * @param property The property object (will be cloned).
     */
    virtual void setNodeProperty(const std::string& name, const Clonable& property)
    {
      if (hasNodeProperty(name))
        delete properties_[name];
      properties_[name] = property.clone();
    }
    
    virtual Clonable* getNodeProperty(const std::string& name) throw (PropertyNotFoundException)
    {
      if (hasNodeProperty(name))
        return properties_[name];
      else
        throw PropertyNotFoundException("", name, this);
    }
    
    virtual const Clonable* getNodeProperty(const std::string& name) const throw (PropertyNotFoundException)
    {
      if (hasNodeProperty(name))
        return const_cast<const Clonable*>(properties_[name]);
      else
        throw PropertyNotFoundException("", name, this);
    }
    
    virtual Clonable* removeNodeProperty(const std::string& name) throw (PropertyNotFoundException)
    {
      if (hasNodeProperty(name))
      {
        Clonable* removed = properties_[name];
        properties_.erase(name);
        return removed;
      }
      else
        throw PropertyNotFoundException("", name, this);
    }
    
    virtual void deleteNodeProperty(const std::string& name) throw (PropertyNotFoundException)
    {
      if (hasNodeProperty(name))
      {
        delete properties_[name];
        properties_.erase(name);
      }
      else
        throw PropertyNotFoundException("", name, this);
    }
    
    /**
     * @brief Remove all node properties.
     *
     * Attached objects will not be deleted.
     */
    virtual void removeNodeProperties()
    {
      properties_.clear();
    }
    
    /**
     * @brief Delete all node properties.
     */
    virtual void deleteNodeProperties()
    {
      for (std::map<std::string, Clonable*>::iterator i = properties_.begin(); i != properties_.end(); i++)
      {
        delete i->second;
      }
      properties_.clear();
    }
    
    virtual bool hasNodeProperty(const std::string& name) const { return properties_.find(name) != properties_.end(); }
    
    virtual std::vector<std::string> getNodePropertyNames() const { return MapTools::getKeys(properties_); }
    
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
    virtual void setBranchProperty(const std::string& name, const Clonable& property)
    {
      getBranch()->setProperty(name, property);
    }
    
    virtual Clonable* getBranchProperty(const std::string& name) throw (PropertyNotFoundException)
    {
      return getBranch()->getProperty(name);
    }
    
    virtual const Clonable* getBranchProperty(const std::string& name) const throw (PropertyNotFoundException)
    {
      return getBranch()->getProperty(name);
    }
    
    virtual Clonable* removeBranchProperty(const std::string& name) throw (PropertyNotFoundException)
    {
      return getBranch()->removeProperty(name);
    }
    
    virtual void deleteBranchProperty(const std::string& name) throw (PropertyNotFoundException)
    {
      getBranch()->deleteProperty(name);
    }
    
    /**
     * @brief Remove all branch properties.
     *
     * Attached objects will not be deleted.
     */
    virtual void removeBranchProperties()
    {
      getBranch()->removeProperties();
    }
    
    /**
     * @brief Delete all branch properties.
     */
    virtual void deleteBranchProperties()
    {
      getBranch()->deleteProperties();
    }
    
    virtual bool hasBranchProperty(const std::string& name) const { return getBranch()->hasProperty(name); }
    
    virtual std::vector<std::string> getBranchPropertyNames() const { return getBranch()->getPropertyNames(); }
    
    virtual bool hasBootstrapValue() const;
    
    virtual double getBootstrapValue() const throw (PropertyNotFoundException);
    /** @} */
    // Equality operator:
    
    virtual bool operator==(const PhyloNode& node) const { return phyloTree_->getNodeGraphid(this) == phyloTree_->getNodeGraphid(&node); }
    
    // Tests:
    
    virtual bool isLeaf() const { return (phyloTree_->getOutgoingNeighbors(this).size() == 0); }
    
    
  }; //end of class node 
  
  
} //end of namespace bpp.

#else
namespace bpp{ class PhyloNode; }

#endif  //_PHYLONODE_H_
