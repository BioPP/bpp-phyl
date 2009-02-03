//
// File: Node.h
// Created by: Julien Dutheil
// Created on: Thu Mar 13 12:03:18 2003
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004)

This software is a computer program whose purpose is to provide classes
for phylogenetic data analysis.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef _NODE_H_
#define _NODE_H_

#include "TreeExceptions.h"

// From the STL:
#include <string>
#include <vector>
#include <map>
#include <iostream>

using namespace std;

// From Utils:
#include <Utils/Clonable.h>
#include <Utils/MapTools.h>
#include <Utils/BppString.h>
#include <Utils/Number.h>

namespace bpp
{

/**
 * @brief The phylogenetic node class.
 * 
 * This class is for use with the TreeTemplate class, an implementation of the Tree interface.
 * TreeTemplates are made made of nodes, instances of this class.
 * Since trees are oriented (rooted), each node has one <i>father node</i> and possibly
 * many <i>son nodes</i>. Leaves are nodes without descendant and root is defined has the without
 * father. Inner nodes will generally contain two descendants (the tree is then called
 * <i>bifurcating</i>), but mutlifurcating trees are also allowed with this kind of description.
 * In the rooted case, each inner node also defines a <i>subtree</i>.
 * This allows to work recursively on trees, which is very convenient in most cases.</p>
 * 
 * This class is made the more general as possible, while keeping it very simple. It contains:</p>
 * - An identity tag, to identity it in the tree;
 * - A name, necessary for leaf nodes, optionnal else;
 * - A pointer toward the father node;
 * - A vector of pointer toward son nodes;
 * - The distance from the father node:
 * - A property map, that may contain any information to link to each node, e.g. bootstrap
 * value or GC content.
 * 
 * Methods are provided to help the building of trees from scratch.
 * Trees are more easily built from root to leaves:
 * The addSon(Node) method adds a node to the list of direct descendants of a 
 * given node. The son node will also have its father set to the current node.
 * There is also a setFather method that enables you to change the pointer
 * toward the parent node. This will however not change the list of descendants
 * of the parent node, you will have to tune it manually.
 * 
 * @see Tree, TreeTemplate
 */
class Node:
  public Clonable
{
      
  protected:
    int                             _id;
    string                        * _name;
    vector<Node *>                  _sons;
    Node                          * _father;
    double                        * _distanceToFather;
    mutable map<string, Clonable *> _nodeProperties;
    mutable map<string, Clonable *> _branchProperties;

  public:
  
    /**
     * @brief Build a new void Node object.
     */
    Node() : _id(0), _name(NULL), _sons(), _father(NULL), _distanceToFather(NULL), _nodeProperties(), _branchProperties() {}
      
    /**
     * @brief Build a new Node with specified id.
     */
    Node(int id) : _id(id), _name(NULL), _sons(), _father(NULL), _distanceToFather(NULL), _nodeProperties(), _branchProperties() {}

    /**
     * @brief Build a new Node with specified name.
     */
    Node(const string & name) : _id(0), _name(new string(name)), _sons(), _father(NULL), _distanceToFather(NULL), _nodeProperties(), _branchProperties() {}

    /**
     * @brief Build a new Node with specified id and name.
     */
    Node(int id, const string & name) : _id(id), _name(new string(name)), _sons(), _father(NULL), _distanceToFather(NULL), _nodeProperties(), _branchProperties() {}

    /**
     * @brief Copy constructor.
     * 
     * @param node The node to copy.
     */
    Node(const Node & node);

    /**
     * @brief Assignation operator.
     *
     * @param node the node to copy.
     * @return A reference toward this node.
     */
    Node & operator=(const Node & node);

    virtual ~Node()
    {  
      delete _name;
      delete _distanceToFather;
      for(map<string, Clonable *>::iterator i = _nodeProperties.begin(); i != _nodeProperties.end(); i++)
        delete i->second;
      for(map<string, Clonable *>::iterator i = _branchProperties.begin(); i != _branchProperties.end(); i++)
        delete i->second;
    }

    Node * clone() const { return new Node(*this); }

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
    virtual int getId() const { return _id; }
    
    /**
     * @brief Set this node's id.
     *
     * @param id The new identity tag.
     */
    virtual void setId(int id) { _id = id; }

    virtual vector<int> getSonsId() const
    {
      vector<int> sonsId(_sons.size());
      for(unsigned int i = 0; i < _sons.size(); i++) {
        sonsId[i] = _sons[i] -> getId();
      }
      return sonsId;
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
    virtual string getName() const throw (NodeException)
    {
      if(!hasName()) throw NodeException("Node::getName: no name associated to this node.", this);
      return * _name;
    }
        
    /**
     * @brief Give a name or update the name associated to the node.
     * 
     * @param name The name to give to the node.
     */
    virtual void setName(const string & name)
    { 
      if(_name) delete _name;
      _name = new string(name);
    }
    
    /**
     * @brief Delete the name associated to this node (do nothing if there is no name).
     */
    virtual void deleteName()
    {
      if(_name) delete _name;
      _name = NULL;
    }
    
    /**
     * @brief Tell is this node has a name.
     * 
     * @return True if name != NULL.
     */
    virtual bool hasName() const { return _name != NULL; }
    
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
    virtual double getDistanceToFather() const throw (NodeException) 
    {
      if(!hasDistanceToFather()) throw NodeException("Node::getDistanceToFather: Node has no distance." , this);
      return * _distanceToFather;
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
    virtual void setDistanceToFather(double distance) { delete _distanceToFather; _distanceToFather = new double(distance); }
    
    /**
     * @brief Delete the distance to the father node.
     */
    virtual void deleteDistanceToFather() { delete _distanceToFather; _distanceToFather = NULL; }
        
    /**
     * @brief Tell is this node has a distance to the father.
     * 
     * @return True if distanceToFather != NULL.
     */
    virtual bool hasDistanceToFather() const { return _distanceToFather != NULL; }

    /** @} */
        
    /**
     * @name Father:
     *
     * @{
     */
      
    /**
     * @brief Get the father of this node is there is one.
     * 
     * @return A pointer toward the father node, NULL if there is not.
     */     
    virtual const Node * getFather() const { return _father; }
        
    /**
     * @brief Get the father of this node is there is one.
     * 
     * @return A pointer toward the father node, NULL if there is not.
     */     
    virtual Node * getFather() { return _father; }
    
    virtual int getFatherId() const { return _father->getId(); }
        
    /**
     * @brief Set the father node of this node.
     * 
     * @param node The father node.
     */
    virtual void setFather(Node & node) { _father = & node; }
        
    /**
     * @brief Remove the father of this node.
     */
    virtual Node * removeFather() { Node * f = _father; _father = NULL; return f; }
        
    /**
     * @brief Tell if this node has a father node.
     */
    virtual bool hasFather() const { return _father != NULL; }

    /** @} */

    /**
     * @name Sons:
     *
     * @{
     */
         
    virtual unsigned int getNumberOfSons() const { return _sons.size(); }

    virtual vector<Node *>& getSons()
    {
      return _sons;
    }

    virtual const Node * getSon(unsigned int pos) const throw (IndexOutOfBoundsException)
    {
      if(pos >= _sons.size()) throw IndexOutOfBoundsException("Node::getSon().", pos, 0, _sons.size()-1);
      return _sons[pos];
    }
      
    virtual Node * getSon(unsigned int pos) throw (IndexOutOfBoundsException)
    {
      if(pos >= _sons.size()) throw IndexOutOfBoundsException("Node::getSon().", pos, 0, _sons.size()-1);
      return _sons[pos];
    }
        
    virtual void addSon(unsigned int pos, Node & node)
    {
      _sons.insert(_sons.begin() + pos, & node);
      node._father = this;
    }

    virtual void addSon(Node & node)
    {
      _sons.push_back(& node);
      node._father = this;
    }

    virtual void setSon(unsigned int pos, Node & node) throw (IndexOutOfBoundsException)
    {
      if(pos >= _sons.size()) throw IndexOutOfBoundsException("Node::setSon(). Invalid node position.", pos, 0, _sons.size()-1);
      _sons[pos] = & node;
      node._father = this;
    }
        
    virtual Node * removeSon(unsigned int pos) throw (IndexOutOfBoundsException)
    {
      if(pos >= _sons.size()) throw IndexOutOfBoundsException("Node::removeSon(). Invalid node position.", pos, 0, _sons.size()-1);
      Node * node = _sons[pos];
      _sons.erase(_sons.begin() + pos);
      node->removeFather();
      return node;
    }
    
    virtual void removeSon(Node & node) throw (NodeNotFoundException)
    {
      for(unsigned int i = 0; i < _sons.size(); i++)
      {
        if(_sons[i] == &node)
        {
          _sons.erase(_sons.begin() + i);
          node.removeFather();
          return;
        }
      }
      throw NodeNotFoundException("Node::removeSon.", node.getId());
    }
    
    virtual void removeSons() {  while(_sons.size() != 0) removeSon(0); }
        
    virtual void swap(unsigned int branch1, unsigned int branch2) throw (IndexOutOfBoundsException);

    virtual unsigned int getSonPosition(const Node & son) const throw (NodeNotFoundException);

    /** @} */

    // These functions must not be declared as virtual!!
    
    vector<const Node *> getNeighbors() const;
    
    vector<Node *> getNeighbors();

    virtual unsigned int degree() const { return getNumberOfSons() + (hasFather() ? 1 : 0); }
    
    /**
     * @name Operators:
     * 
     * - a positive value send the corresponding son;
     * - a negative value send the father.
     *
     * @{
     */
         
    Node * operator[](int i) { return (i < 0) ? _father : _sons[i]; }
        
    const Node * operator[](int i) const { return (i < 0) ? _father : _sons[i]; }
    
    /** @} */
    
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
    virtual void setNodeProperty(const string & name, const Clonable & property)
    {
      if(hasNodeProperty(name))
        delete _nodeProperties[name];
      _nodeProperties[name] = property.clone();
    }
        
    virtual Clonable * getNodeProperty(const string & name) throw (PropertyNotFoundException)
    {
      if(hasNodeProperty(name))
        return _nodeProperties[name];
      else
        throw PropertyNotFoundException("", name, this);
    }
        
    virtual const Clonable * getNodeProperty(const string & name) const throw (PropertyNotFoundException)
    {
      if(hasNodeProperty(name))
        return const_cast<const Clonable *>(_nodeProperties[name]);
      else
        throw PropertyNotFoundException("", name, this);
    }
        
    virtual Clonable * removeNodeProperty(const string & name) throw (PropertyNotFoundException)
    {
      if(hasNodeProperty(name))
      {
        Clonable * removed = _nodeProperties[name];
        _nodeProperties.erase(name);
        return removed;
      }
      else
        throw PropertyNotFoundException("", name, this);
    }  
        
    virtual void deleteNodeProperty(const string & name) throw (PropertyNotFoundException)
    {
      if(hasNodeProperty(name))
      {
        delete _nodeProperties[name];
        _nodeProperties.erase(name);
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
      _nodeProperties.clear();
    }  

    /**
     * @brief Delete all node properties.
     */
    virtual void deleteNodeProperties()
    {
      for(map<string, Clonable *>::iterator i = _nodeProperties.begin(); i != _nodeProperties.end(); i++)
        delete i->second; 
      _nodeProperties.clear();
    }  
                
    virtual bool hasNodeProperty(const string & name) const { return _nodeProperties.find(name) != _nodeProperties.end(); }

    virtual vector<string> getNodePropertyNames() const { return MapTools::getKeys(_nodeProperties); }
    
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
    virtual void setBranchProperty(const string & name, const Clonable & property)
    {
      if(hasBranchProperty(name))
        delete _branchProperties[name];
      _branchProperties[name] = property.clone();
    }
        
    virtual Clonable * getBranchProperty(const string & name) throw (PropertyNotFoundException)
    {
      if(hasBranchProperty(name))
        return _branchProperties[name];
      else
        throw PropertyNotFoundException("", name, this);
    }
        
    virtual const Clonable * getBranchProperty(const string & name) const throw (PropertyNotFoundException)
    {
      if(hasBranchProperty(name))
        return const_cast<const Clonable *>(_branchProperties[name]);
      else
        throw PropertyNotFoundException("", name, this);
    }
        
    virtual Clonable * removeBranchProperty(const string & name) throw (PropertyNotFoundException)
    {
      if(hasBranchProperty(name))
      {
        Clonable * removed = _branchProperties[name];
        _branchProperties.erase(name);
        return removed;
      }
      else
        throw PropertyNotFoundException("", name, this);
    }  
    
    virtual void deleteBranchProperty(const string & name) throw (PropertyNotFoundException)
    {
      if(hasBranchProperty(name))
      {
        delete _branchProperties[name];
        _branchProperties.erase(name);
      }
      else
        throw PropertyNotFoundException("", name, this);
    }  
    
    /**
     * @brief Remove all branch properties.
     *
     * Attached objects will not be deleted.
     */
    virtual void removeBranchProperties()
    {
      _branchProperties.clear();
    }  

    /**
     * @brief Delete all branch properties.
     */
    virtual void deleteBranchProperties()
    {
      for(map<string, Clonable *>::iterator i = _branchProperties.begin(); i != _branchProperties.end(); i++)
        delete i->second; 
      _branchProperties.clear();
    }  
        
    virtual bool hasBranchProperty(const string & name) const { return _branchProperties.find(name) != _branchProperties.end(); }

    virtual vector<string> getBranchPropertyNames() const { return MapTools::getKeys(_branchProperties); }
   
    virtual double getBootstrapValue() const throw (PropertyNotFoundException);
    /** @} */
    // Equality operator:

    virtual bool operator==(const Node & node) const { return _id == node._id; }  
        
    // Tests:

    virtual bool isLeaf() const { return degree() == 1; }

};

} //end of namespace bpp.

#endif  //_NODE_H_

