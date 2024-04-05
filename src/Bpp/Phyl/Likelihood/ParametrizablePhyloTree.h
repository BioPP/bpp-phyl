// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

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
      throw Exception("ParametrizablePhyloTree::setMinimumBranchLength. Minimum branch length should be lower than the maximum one: " + TextTools::toString(maximumBrLen_));
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
      throw Exception("ParametrizablePhyloTree::setMaximumBranchLength. Maximum branch length should be higher than the minimum one: " + TextTools::toString(minimumBrLen_));
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
