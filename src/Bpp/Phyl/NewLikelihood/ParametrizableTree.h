//
// File: ParametrizableTree.h
// Created by: Julien Dutheil
// Created on: Wed Jul 11 20:36 2012
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

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

#ifndef _PARAMETRIZABLETREE_H_
#define _PARAMETRIZABLETREE_H_

#include <Bpp/Numeric/AbstractParametrizable.h>

//From the stl:
#include <string>

namespace bpp
{

  class ParametrizableTree:
    public AbstractParametrizable
  {
    private:
      TreeTemplate<Node> tree_;
      bool liveIndex_;
      std::map<int, unsigned int> index_;
      std::map<std::string, Node*> reverseIndex_; //Watch out when copying!
      bool isSynchronized_; //If no live index is performed, record any parameter change that has been made.

      double minimumBrLen_;
      double maximumBrLen_;
      std::auto_ptr<Constraint> brLenConstraint_;

    public:
      ParametrizableTree(const Tree& tree, bool liveIndex = false, const std::string& prefix = "");

      ParametrizableTree(const ParametrizableTree& pTree);

      ParametrizableTree& operator=(const ParametrizableTree& pTree);

    public:
      const Parameter& getBranchLengthParameter(int nodeId) const throw (NodeNotFoundException) {
        std::map<int, unsigned int>::const_iterator it = index_.find(nodeId);
        if (it != index_.end())
          return getParameter_(it->second);
        else
          throw NodeNotFoundException("ParametrizableTree::getBranchLengthParameter.", nodeId);
      }

      const TreeTemplate<Node>& getTree() const;

    virtual void setMinimumBranchLength(double minimum) throw (Exception)
    {
      if (minimum > maximumBrLen_)
        throw Exception("ParametrizableTree::setMinimumBranchLength. Minimum branch length sould be lower than the maximum one: " + TextTools::toString(maximumBrLen_));
      minimumBrLen_ = minimum;
      if (brLenConstraint_.get()) brLenConstraint_.release();
      brLenConstraint_.reset(new IntervalConstraint(minimumBrLen_, maximumBrLen_, true, true));
      resetParameters_();
      buildIndex_(*tree_.getRootNode());
    }

    virtual void setMaximumBranchLength(double maximum) throw (Exception)
    {
      if (maximum < minimumBrLen_)
        throw Exception("ParametrizableTree::setMaximumBranchLength. Maximum branch length sould be higher than the minimum one: " + TextTools::toString(minimumBrLen_));
      maximumBrLen_ = maximum;
      if (brLenConstraint_.get()) brLenConstraint_.release();
      brLenConstraint_.reset(new IntervalConstraint(minimumBrLen_, maximumBrLen_, true, true));
      resetParameters_();
      buildIndex_(*tree_.getRootNode());
    }

    virtual double getMinimumBranchLength() const { return minimumBrLen_; }
    virtual double getMaximumBranchLength() const { return maximumBrLen_; }

    private:
      void buildIndex_(Node& node);
      void buildReverseIndex_(Node* node);
      void updateTreeFromParameters_();
      void fireParameterChanged (const ParameterList& parameters);
  };

} //end of namespace bpp

#endif //_PARAMETRIZABLETREE_H_

