//
// File: SubstitutionMapping.h
// Created by: Julien Dutheil
// Created on: Wed Apr 5 09:51 2005
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004, 2005, 2006)

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

#ifndef _SUBSTITUTIONMAPPING_H_
#define _SUBSTITUTIONMAPPING_H_

#include "../Tree.h"
#include "../TreeTemplate.h"

#include <Bpp/Clonable.h>

//From the STL:
#include <vector>
#include <memory>

namespace bpp
{

/**
 * @brief General interface for storing mapping data.
 *
 * There are several kinds of mapping:
 * - Exact mapping, storing the positions of each substitution onto each branch,
 * - Probabilistic mapping, storing the number of substitutions onto each branch.
 *
 * Since only probabilistic substitution mapping is implemented for now, the basal 
 * interfac only contains one method.
 * More methods are expected to be added later.
 */
class SubstitutionMapping:
  public virtual Clonable
{

  public:
    SubstitutionMapping() {}
    virtual ~SubstitutionMapping() {}

#ifndef NO_VIRTUAL_COV
    SubstitutionMapping* clone() const = 0;
#endif

  public:
    
    /**
     * @return Get the phylogenetic tree associated to this mapping.
     */
    virtual const Tree& getTree() const = 0;
    
    /**
     * @return True is the map is empty, that is, if no tree is associated to the map yet.
     */
    virtual bool isEmpty() const = 0;

    /**
     * @return The number of sites mapped.
     */
    virtual unsigned int getNumberOfSites() const = 0;
    
    /**
     * @return The number of branches mapped.
     */
    virtual unsigned int getNumberOfBranches() const = 0;
    
    /**
     * @param index The site index.
     * @return The site position corresponding to the index.
     */
    virtual int getSitePosition(unsigned int index) const = 0;
    
    /**
     * @return A vector with all tree branch lengths.
     */
    virtual std::vector<double> getBranchLengths() const = 0;
    
    /**
     * @param nodeId An id of the node to look for in the map.
     * @return The mapping index for the specified node id.
     */
    virtual unsigned int getNodeIndex(int nodeId) const throw (NodeNotFoundException) = 0;

    /**
     * @brief Set the position of a given site.
     *
     * @warning No index checking is performed, use with care!
     * @param index The site index.
     * @param position The position of the site.
     */
    virtual void setSitePosition(unsigned int index, int position) = 0;

};






/**
 * @brief Partial implementation of the substitution mapping interface.
 *
 * This implementation copies the input tree in a TreeTemplate<Node> object.
 */
class AbstractSubstitutionMapping:
  public SubstitutionMapping
{
  private:
    std::auto_ptr<const TreeTemplate<Node> > tree_;
    std::vector<int> sitesPositions_;
    std::vector<const Node *> nodes_;
    unsigned int nbSites_;
    unsigned int nbBranches_;

  public:
    AbstractSubstitutionMapping() : tree_(0), sitesPositions_(), nodes_(), nbSites_(0), nbBranches_(0) {}

    AbstractSubstitutionMapping(const Tree& tree) : tree_(new TreeTemplate<Node>(tree)), sitesPositions_(), nodes_(), nbSites_(0), nbBranches_(0)
    {
      nodes_ = tree_->getNodes();
      nodes_.pop_back(); // remove root node.
      nbBranches_ = nodes_.size();
    }

    AbstractSubstitutionMapping(const AbstractSubstitutionMapping& absm):
      tree_(dynamic_cast<const TreeTemplate<Node>*>(absm.tree_->clone())),
      sitesPositions_(absm.sitesPositions_),
      nodes_(),
      nbSites_(absm.nbSites_),
      nbBranches_(absm.nbBranches_)
    {
      nodes_ = tree_->getNodes();
      nodes_.pop_back(); // remove root node.
    }

    AbstractSubstitutionMapping& operator=(const AbstractSubstitutionMapping& absm)
    {
      tree_.reset(dynamic_cast<const TreeTemplate<Node>*>(absm.tree_->clone()));
      sitesPositions_ = absm.sitesPositions_;
      nbSites_        = absm.nbSites_;
      nbBranches_     = absm.nbBranches_;
      nodes_          = tree_->getNodes();
      nodes_.pop_back(); // remove root node.
      return *this;
    }

    virtual ~AbstractSubstitutionMapping() {}

  public:

    bool isEmpty() const { return tree_.get() == 0; }

		const	TreeTemplate<Node>& getTree() const throw (Exception)
    {
      if (isEmpty()) throw Exception("AbstractSubstitutionMapping::getSitePosition. No tree is assigned to this map yet.");
      return *tree_.get();
    }

    void setTree(const Tree& tree)
    {
      tree_.reset(new TreeTemplate<Node>(tree));
      nodes_ = tree_->getNodes();
      nodes_.pop_back(); // remove root node.
      nbBranches_ = nodes_.size();
    }
 
    int getSitePosition(unsigned int index) const throw (Exception)
    {
      if (isEmpty()) throw Exception("AbstractSubstitutionMapping::getSitePosition. No tree is assigned to this map yet.");
      return sitesPositions_[index];
    }
    
    void setSitePosition(unsigned int index, int position) throw (Exception)
    {
      if (isEmpty()) throw Exception("AbstractSubstitutionMapping::setSitePosition. No tree is assigned to this map yet.");
      sitesPositions_[index] = position;
    }
		
    virtual unsigned int getNumberOfSites() const { return nbSites_; }

    virtual unsigned int getNumberOfBranches() const { return nbBranches_; }
     
    virtual const Node* getNode(unsigned int nodeIndex) const { return nodes_[nodeIndex]; }

    virtual void setNumberOfSites(unsigned int numberOfSites)
    {
      nbSites_ = numberOfSites;
      sitesPositions_.resize(numberOfSites);
      for (unsigned int i = 0; i < numberOfSites; i++)
        sitesPositions_[i] = i + 1; //Default: all sizes numbered for 1 to n.
    }

    virtual std::vector<double> getBranchLengths() const
    {
      std::vector<double> brLen(nbBranches_);
      for (unsigned int i = 0; i < nbBranches_; i++)
        brLen[i] = nodes_[i]->getDistanceToFather();
      return brLen;
    }

    virtual unsigned int getNodeIndex(int nodeId) const throw (NodeNotFoundException)
    {
      for (unsigned int i = 0; i < nbBranches_; i++)
        if(nodes_[i]->getId() == nodeId) return i;
      throw NodeNotFoundException("ProbabilisticSubstitutionMapping::getNodeIndex(nodeId).", TextTools::toString(nodeId));
    }


};

} //end of namespace bpp.

#endif //_SUBSTITUTIONMAPPING_H_

