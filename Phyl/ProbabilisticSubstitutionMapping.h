//
// File: ProbabilisticSubstitutionMapping.h
// Created by: Julien Dutheil
// Created on: Wed Apr 5 10:47 2006
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

#ifndef _PROBABILISTICSUBSTITUTIONMAPPING_H_
#define _PROBABILISTICSUBSTITUTIONMAPPING_H_

#include "SubstitutionMapping.h"
#include "TreeExceptions.h"

// From the STL:
#include <vector>

// From Utils:
#include <Utils/TextTools.h>

namespace bpp
{

/**
 * @brief Data storage class for probabilistic substitution mappings.
 *
 * A 'probabilistic' mapping contains the expected number of substitutions for all branches and all sites.
 */
class ProbabilisticSubstitutionMapping:
  public AbstractSubstitutionMapping
{
  protected:
    /**
     * @brief Substitution numbers storage.
     *
     * Numbers are stored by sites.
     */
    vector< vector<double> > _mapping;
    vector<const Node *> _nodes;
    unsigned int _nbSites;
    unsigned int _nbBranches;
  
  public:
    
    /**
     * @brief Build a new ProbabilisticSubstitutionMapping object.
     *
     * @param tree The tree object to use. It will be cloned for internal use.
     * @param numberOfSites The number of sites to map.
     */
    ProbabilisticSubstitutionMapping(const Tree & tree, unsigned int numberOfSites);
    /**
     * @brief Build a new ProbabilisticSubstitutionMapping object.
     *
     * @param tree The tree object to use. It will be cloned for internal use.
     */
    ProbabilisticSubstitutionMapping(const Tree & tree);
    
    /**
     * @brief Copy constructor: clone tree.
     */
    ProbabilisticSubstitutionMapping(const ProbabilisticSubstitutionMapping & psm):
      _mapping(psm._mapping),
      _nbSites(psm._nbSites),
      _nbBranches(psm._nbBranches)
    {
      if(psm._tree == NULL)
      {
        _tree = NULL;
      } else { 
        _tree = new TreeTemplate<Node>(*psm._tree);
        _nodes = _tree->getNodes();
        _nodes.pop_back(); // remove root node.
      }
    }

    ProbabilisticSubstitutionMapping & operator=(const ProbabilisticSubstitutionMapping & psm)
    {
      _mapping    = psm._mapping;
      _nbSites    = psm._nbSites;
      _nbBranches = psm._nbBranches;
      if(psm._tree == NULL)
      {
        _tree = NULL;
      } else {
        _tree = new TreeTemplate<Node>(*psm._tree);
        _nodes = _tree->getNodes();
        _nodes.pop_back(); // remove root node.
      }
      return *this;
    }

#ifdef NO_VIRTUAL_COV
    Clonable *
#else
    ProbabilisticSubstitutionMapping *
#endif
    clone() const { return new ProbabilisticSubstitutionMapping(*this); }

    virtual ~ProbabilisticSubstitutionMapping()
    {
      delete _tree;
    }

  public:

    unsigned int getNumberOfSites() const { return _nbSites; }

    unsigned int getNumberOfBranches() const { return _nbBranches; }
     
    virtual double getNumberOfSubstitutions(int nodeId, unsigned int siteIndex) const
    {
      return _mapping[siteIndex][getNodeIndex(nodeId)];
    }
    
    virtual const Node * getNode(unsigned int nodeIndex) const { return _nodes[nodeIndex]; }

    virtual vector<double> getBranchLengths() const
    {
      vector<double> brLen(_nbBranches);
      for(unsigned int i = 0; i < _nbBranches; i++)
        brLen[i] = _nodes[i]->getDistanceToFather();
      return brLen;
    }

    virtual unsigned int getNodeIndex(int nodeId) const throw (NodeNotFoundException)
    {
      for(unsigned int i = 0; i < _nodes.size(); i++)
        if(_nodes[i]->getId() == nodeId) return i;
      throw NodeNotFoundException("ProbabilisticSubstitutionMapping::getNodeIndex(nodeId).", TextTools::toString(nodeId));
    }

    /**
     * @brief (Re)-set the phylogenetic tree associated to this mapping.
     *
     * @param tree The new tree.
     */
    virtual void setTree(const Tree & tree);

    virtual void setNumberOfSites(unsigned int numberOfSites);
    
    /**
     * @brief Direct access to substitution numbers.
     *
     * @warning No index checking is performed, use with care!
     */
    virtual double & operator()(unsigned int nodeIndex, unsigned int siteIndex)
    {
      return _mapping[siteIndex][nodeIndex];
    }

    /**
     * @brief Direct access to substitution numbers.
     *
     * @warning No index checking is performed, use with care!
     */
    virtual const double & operator()(unsigned int nodeIndex, unsigned int siteIndex) const
    {
      return _mapping[siteIndex][nodeIndex];
    }
     
    /**
     * @brief Direct access to substitution numbers.
     *
     * @warning No index checking is performed, use with care!
     */
    vector<double> & operator[](unsigned int siteIndex)
    {
      return _mapping[siteIndex];
    }

    /**
     * @brief Direct access to substitution numbers.
     *
     * @warning No index checking is performed, use with care!
     */
    const vector<double> & operator[](unsigned int siteIndex) const
    {
      return _mapping[siteIndex];
    }

    /**
     * @brief Set the position of a given site.
     *
     * @warning No index checking is performed, use with care!
     * @param index The site index.
     * @param position The position of the site.
     */
    void setSitePosition(unsigned int index, int position) { _sitesPostions[index] = position; }
};

} //end of namespace bpp.

#endif //_PROBABILISTICSUBSTITUTIONMAPPING_H_

