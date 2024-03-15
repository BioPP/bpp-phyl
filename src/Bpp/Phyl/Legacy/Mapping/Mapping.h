// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _LEGACY_MAPPING_H_
#define _LEGACY_MAPPING_H_

#include <Bpp/Clonable.h>

#include "../../Tree/Tree.h"
#include "../../Tree/TreeTemplate.h"

//From the STL:
#include <vector>
#include <memory>

namespace bpp
{

  /**
   * @brief General interface for storing mapping data.
   */
  class LegacyMappingInterface:
    public virtual Clonable
  {

  public:
    LegacyMappingInterface() {}
    virtual ~LegacyMappingInterface() {}

    LegacyMappingInterface* clone() const override = 0;

  public:
    
    /**
     * @return Get the phylogenetic tree associated to this mapping.
     */
    virtual const Tree& tree() const = 0;
    
    /**
     * @return True is the map is empty, that is, if no tree is associated to the map yet.
     */
    virtual bool isEmpty() const = 0;

    /**
     * @return The number of sites mapped.
     */
    virtual size_t getNumberOfSites() const = 0;
    
    /**
     * @return The number of branches mapped.
     */
    virtual size_t getNumberOfBranches() const = 0;
    
    /**
     * @param index The site index.
     * @return The site position corresponding to the index.
     */
    virtual int getSitePosition(size_t index) const = 0;
    
    /**
     * @return A vector with all tree branch lengths.
     */
    virtual std::vector<double> getBranchLengths() const = 0;
    
    /**
     * @param nodeId An id of the node to look for in the map.
     * @return The mapping index for the specified node id.
     */
    virtual size_t getNodeIndex(int nodeId) const = 0;

    /**
     * @brief Set the position of a given site.
     *
     * @warning No index checking is performed, use with care!
     * @param index The site index.
     * @param position The position of the site.
     */
    virtual void setSitePosition(size_t index, int position) = 0;
  };






  /**
   * @brief Partial implementation of the mapping interface.
   *
   * This implementation copies the input tree in a TreeTemplate<Node> object.
   */
  class LegacyAbstractMapping:
    public virtual LegacyMappingInterface
  {
  private:
    std::unique_ptr<const TreeTemplate<Node>> tree_;
    std::vector<int> sitesPositions_;
    std::vector<const Node *> nodes_;
    size_t nbSites_;
    size_t nbBranches_;

  public:

    LegacyAbstractMapping(const Tree& tree) : tree_(new TreeTemplate<Node>(tree)), sitesPositions_(), nodes_(), nbSites_(0), nbBranches_(0)
    {
      nodes_ = tree_->getNodes();
      nodes_.pop_back(); // remove root node.
      nbBranches_ = nodes_.size();
    }

    LegacyAbstractMapping(const LegacyAbstractMapping& absm):
      tree_(dynamic_cast<const TreeTemplate<Node>*>(absm.tree_->clone())),
      sitesPositions_(absm.sitesPositions_),
      nodes_(),
      nbSites_(absm.nbSites_),
      nbBranches_(absm.nbBranches_)
    {
      nodes_ = tree_->getNodes();
      nodes_.pop_back(); // remove root node.
    }

    LegacyAbstractMapping& operator=(const LegacyAbstractMapping& absm)
    {
      tree_.reset(dynamic_cast<const TreeTemplate<Node>*>(absm.tree_->clone()));
      sitesPositions_ = absm.sitesPositions_;
      nbSites_        = absm.nbSites_;
      nbBranches_     = absm.nbBranches_;
      nodes_          = tree_->getNodes();
      nodes_.pop_back(); // remove root node.
      return *this;
    }

    LegacyAbstractMapping* clone() const override = 0;
  
    virtual ~LegacyAbstractMapping() {}
  
  public:
  
    bool isEmpty() const override { return !tree_; }

    const TreeTemplate<Node>& tree() const override
    {
      if (isEmpty()) throw Exception("LegacyAbstractSubstitutionMapping::getSitePosition. No tree is assigned to this map yet.");
      return *tree_;
    }
  
    void setTree(const Tree& tree)
    {
      tree_.reset(new TreeTemplate<Node>(tree));
      nodes_ = tree_->getNodes();
      nodes_.pop_back(); // remove root node.
      nbBranches_ = nodes_.size();
    }
 
    int getSitePosition(size_t index) const override
    {
      if (isEmpty()) throw Exception("LegacyAbstractMapping::getSitePosition. No tree is assigned to this map yet.");
      return sitesPositions_[index];
    }
    
    void setSitePosition(size_t index, int position) override
    {
      if (isEmpty()) throw Exception("LegacyAbstractMapping::setSitePosition. No tree is assigned to this map yet.");
      sitesPositions_[index] = position;
    }
		
    size_t getNumberOfSites() const override { return nbSites_; }

    size_t getNumberOfBranches() const override { return nbBranches_; }
    
    virtual const Node* getNode(size_t nodeIndex) const { return nodes_[nodeIndex]; }

    virtual void setNumberOfSites(size_t numberOfSites)
    {
      nbSites_ = numberOfSites;
      sitesPositions_.resize(numberOfSites);
      for (size_t i = 0; i < numberOfSites; i++)
        sitesPositions_[i] = static_cast<int>(i + 1); //Default: all sizes numbered for 1 to n.
    }

    virtual std::vector<double> getBranchLengths() const override
    {
      std::vector<double> brLen(nbBranches_);
      for (size_t i = 0; i < nbBranches_; i++)
        brLen[i] = nodes_[i]->getDistanceToFather();
      return brLen;
    }

    virtual size_t getNodeIndex(int nodeId) const override
    {
      for (size_t i = 0; i < nbBranches_; i++)
        if(nodes_[i]->getId() == nodeId) return i;
      throw NodeNotFoundException("LegacyAbstractMapping::getNodeIndex(nodeId).", TextTools::toString(nodeId));
    }


  };

} //end of namespace bpp.

#endif //_LEGACY_MAPPING_H_

