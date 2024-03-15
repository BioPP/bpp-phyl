// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_FORMULAOFPHYLOLIKELIHOOD_H
#define BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_FORMULAOFPHYLOLIKELIHOOD_H

#include <Bpp/Numeric/Function/Operators/ComputationTree.h>

#include "PhyloLikelihoodSet.h"

namespace bpp
{
/**
 * @brief The PhyloLikelihoodFormula class, for phylogenetic
 * likelihood on several independent data.
 *
 * WARNING: This formula applies on the log-likelihoods (ie getValues())
 *
 */
class PhyloLikelihoodFormula :
  public AbstractPhyloLikelihoodSet
{
private:
  std::unique_ptr<ComputationTree> compTree_;

  std::shared_ptr<LikelihoodCalculation> likCal_;

public:
  PhyloLikelihoodFormula(Context& context, std::shared_ptr<PhyloLikelihoodContainer> pC, const std::string& formula, bool inCollection = true);

  virtual ~PhyloLikelihoodFormula() {}

protected:

  PhyloLikelihoodFormula(const PhyloLikelihoodFormula& sd):
    AbstractPhyloLikelihood(sd),
    AbstractParametrizable(sd),
    AbstractPhyloLikelihoodSet(sd),
    compTree_(sd.compTree_->clone()),
    likCal_(sd.likCal_)
  {}

  PhyloLikelihoodFormula& operator=(const PhyloLikelihoodFormula& sd)
  {
    AbstractPhyloLikelihoodSet::operator=(sd);
    compTree_.reset(sd.compTree_->clone());
    likCal_ = sd.likCal_;
    return *this;
  }


  PhyloLikelihoodFormula* clone() const
  {
    return new PhyloLikelihoodFormula(*this);
  }

public:
  /**
   * @ input
   *
   */

  void readFormula(const std::string& formula, bool inCollection = true);

  /**
   * @ output
   *
   */

  std::string output() const;

  /**
   * @name The likelihood functions.
   */
  LikelihoodCalculation& likelihoodCalculation() const
  {
    return *likCal_;
  }

  std::shared_ptr<LikelihoodCalculation> getLikelihoodCalculation() const
  {
    return likCal_;
  }

private:
  /**
   * @brief Build the LikelihoodNode from the computation Tree
   *
   */
  ValueRef<DataLik> makeLikelihoods()
  {
    return makeLikelihoodsFromOperator(compTree_->getRoot());
  }

  /**
   * @brief Build the LikelihoodNode from a node of the computation Tree
   *
   */

  ValueRef<DataLik> makeLikelihoodsFromOperator(std::shared_ptr<Operator> op);
};
} // end of namespace bpp.
#endif // BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_FORMULAOFPHYLOLIKELIHOOD_H
