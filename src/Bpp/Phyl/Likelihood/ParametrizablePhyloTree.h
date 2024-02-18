//
// File: ParametrizablePhyloTree.h
// Authors:
//   Laurent Guéguen
// Created: jeudi 15 septembre 2016, ÃÂ  06h 40
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

#ifndef BPP_PHYL_LIKELIHOOD_PARAMETRIZABLEPHYLOTREE_H
#define BPP_PHYL_LIKELIHOOD_PARAMETRIZABLEPHYLOTREE_H

#include <Bpp/Graph/AssociationTreeGraphImplObserver.h>
#include <Bpp/Numeric/AbstractParametrizable.h>

#include "../Tree/PhyloBranchParam.h"
#include "../Tree/PhyloNode.h"
#include "../Tree/PhyloTree.h"

// From the stl:
#include <string>

namespace bpp
{
/**
 * @brief PhyloTree with Parametrizable Phylo Branches. They SHARE their
 * branch length parameters.
 */
class ParametrizablePhyloTree :
  public AssociationTreeGlobalGraphObserver<PhyloNode, PhyloBranchParam>,
  public AbstractParametrizable
{
private:
  double minimumBrLen_;
  double maximumBrLen_;
  std::shared_ptr<ConstraintInterface> brLenConstraint_;

public:
  ParametrizablePhyloTree(const PhyloTree& tree, const std::string& prefix = "");

  ParametrizablePhyloTree(const ParametrizablePhyloTree& pTree);

  ParametrizablePhyloTree& operator=(const ParametrizablePhyloTree& pTree);

  ParametrizablePhyloTree* clone() const { return new ParametrizablePhyloTree(*this); }

public:
  std::vector<std::string> getAllLeavesNames() const;

  Vdouble getBranchLengths() const;

  virtual void setMinimumBranchLength(double minimum)
  {
    if (minimum > maximumBrLen_)
      throw Exception("ParametrizablePhyloTree::setMinimumBranchLength. Minimum branch length sould be lower than the maximum one: " + TextTools::toString(maximumBrLen_));
    minimumBrLen_ = minimum;
    if (!brLenConstraint_)
    {
      brLenConstraint_ = std::make_shared<IntervalConstraint>(minimumBrLen_, maximumBrLen_, true, true);
      resetParameters_();
    }
    else
      dynamic_cast<IntervalConstraint&>(*brLenConstraint_).setLowerBound(minimumBrLen_, false);
  }

  virtual void setMaximumBranchLength(double maximum)
  {
    if (maximum < minimumBrLen_)
      throw Exception("ParametrizablePhyloTree::setMaximumBranchLength. Maximum branch length sould be higher than the minimum one: " + TextTools::toString(minimumBrLen_));
    maximumBrLen_ = maximum;
    if (!brLenConstraint_)
    {
      brLenConstraint_ = std::make_shared<IntervalConstraint>(minimumBrLen_, maximumBrLen_, true, true);
      resetParameters_();
    }
    else
      dynamic_cast<IntervalConstraint&>(*brLenConstraint_).setUpperBound(minimumBrLen_, false);
  }

  virtual double getMinimumBranchLength() const { return minimumBrLen_; }
  virtual double getMaximumBranchLength() const { return maximumBrLen_; }

private:
  void fireParameterChanged (const ParameterList& parameters);
};
} // end of namespace bpp
#endif // BPP_PHYL_LIKELIHOOD_PARAMETRIZABLEPHYLOTREE_H
