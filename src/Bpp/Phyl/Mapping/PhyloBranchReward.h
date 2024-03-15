// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MAPPING_PHYLOBRANCHREWARD_H
#define BPP_PHYL_MAPPING_PHYLOBRANCHREWARD_H

#include <Bpp/Clonable.h>
#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/Number.h>

#include "../Tree/PhyloBranch.h"

namespace bpp
{
/*
 * @brief A branch with countings.
 *
 * WARNING : this class does not know anything about site
 * compression, if any. If there are site patterns, they are
 * available in ProbabilisticSubstitutionReward class.
 *
 */

class PhyloBranchReward :
  public PhyloBranch
{
protected:
  /*
   * @brief rewards are stored by site
   *
   */

  Vdouble rewards_;

public:
  /**
   * @brief Constructors.
   *
   * @warning phyloTree_ does not know the edge exists.
   *
   */

  PhyloBranchReward() :
    PhyloBranch(),
    rewards_()
  {}

  PhyloBranchReward(double length) :
    PhyloBranch(length),
    rewards_()
  {}

  PhyloBranchReward(const PhyloBranch& branch) :
    PhyloBranch(branch),
    rewards_()
  {}

  /**
   * @brief Copy constructor.
   *
   * @param branch The branch to copy.
   */

  PhyloBranchReward(const PhyloBranchReward& branch) :
    PhyloBranch(branch),
    rewards_(branch.rewards_)
  {}

  /**
   * @brief Assignation operator.
   *
   * @param branch the branch to copy.
   * @return A reference toward this branch.
   */
  PhyloBranchReward& operator=(const PhyloBranchReward& branch)
  {
    PhyloBranch::operator=(branch);
    rewards_ = branch.rewards_;
    return *this;
  }

  PhyloBranchReward* clone() const { return new PhyloBranchReward(*this); }

  /**
   * @brief destructor. In Graph, nothing is changed.
   *
   */
  ~PhyloBranchReward()
  {}

  /**
   * @brief Sets a number of sites.
   */
  void setNumberOfSites(size_t nbSites)
  {
    rewards_.resize(nbSites);
  }

  /**
   * @brief Gets the number of sites.
   */
  size_t getNumberOfSites() const
  {
    return rewards_.size();
  }


  double getSiteReward(size_t site) const
  {
    if (site >= getNumberOfSites())
      throw BadSizeException("PhyloBranchReward::getSiteReward : bad site number", site, getNumberOfSites());
    return rewards_[site];
  }

  void setSiteReward(size_t site, double rew)
  {
    if (site >= getNumberOfSites())
      throw BadSizeException("PhyloBranchReward::setSiteReward : bad site number", site, getNumberOfSites());
    rewards_[site] = rew;
  }

  /**
   * @brief Gets the rewards at a given site on a given type
   *
   */

  /**
   * @brief Without check
   *
   */
  double operator()(size_t site) const
  {
    return rewards_[site];
  }

  double& operator()(size_t site)
  {
    return rewards_[site];
  }

  /**
   * @brief return rewards
   *
   */
  const Vdouble& getRewards() const
  {
    return rewards_;
  }

  Vdouble& getRewards()
  {
    return rewards_;
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MAPPING_PHYLOBRANCHREWARD_H
