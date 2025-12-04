// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_SIMULATION_GIVENDATASUBSTITUTIONPROCESSSITESIMULATOR_H
#define BPP_PHYL_SIMULATION_GIVENDATASUBSTITUTIONPROCESSSITESIMULATOR_H

#include <Bpp/Phyl/Likelihood/DataFlow/LikelihoodCalculationSingleProcess.h>

#include "SimpleSubstitutionProcessSiteSimulator.h"

namespace bpp
{
/**
 * @brief Site simulation under a unique substitution process, given data.
 *
 * So it will need a LikelihoodCalculationSingleProcess to have site
 * specific likelihoods and manage a posteriori transition
 * probabilities.
 *
 *
 * Transition probabilities are computed a posteriori: On an edge with
 * states x -> y
 *
 *  @f$P(y|x,D) = P(y|x) * P(D\_ | y) / (sum_y' P(y'|x) * P(D\_ | y'))@f$ where @f$D\_@f$ is the data below son node
 *
 * Mixture probabilities are computed a posteriori: On a node n, to
 * choose an edge e among all outgoing edges of n:
 *
 *  @f$P(e|D) = P(e) * P(D\_ | e) / (sum_e' P(e') * P(D\_ | e'))@f$
 *  where @f$D\_@f$ is the data below son node and @f$P(D\_ | e')@f$ is
 *  the likelihood of @f$D\_@f$ under the TOP of edge e':
 *
 *  @f$P(D\_ | e') = 1/N . \sum_x P(D\_ | x)@f$ for @f$x@f$ all states
 *  at the top of edge @f$e'@f$ and @f$N@f$ the number of states
 */
class GivenDataSubstitutionProcessSiteSimulator :
  public SimpleSubstitutionProcessSiteSimulator
{
private:
  std::shared_ptr<LikelihoodCalculationSingleProcess> calcul_;

  /**
   * @brief Position of the copied site, in SHRUNKED data
   */
  Eigen::Index pos_;

  /**
   * @brief Vector of branch indexes where the data is not used, and
   * substitution probabilities are computed from the prior process.
   */
  Vuint vPriorBranch_;

public:
  
  /**
   * @brief Build a Site Simulator of histories from the a posteriori likelihoods at a given site
   *
   * @param calcul the posterior likelihood calculation
   * @param pos the position of the site to imitate
   * @param shrunked if the given position is on the shrunked data (default: false)
   * @param vPrior the vector of branch numbers where the posterior is not used.
   */
  GivenDataSubstitutionProcessSiteSimulator(std::shared_ptr<LikelihoodCalculationSingleProcess> calcul, size_t pos, bool shrunked = false, Vuint vPrior = Vuint()) :
    SimpleSubstitutionProcessSiteSimulator(calcul->getSubstitutionProcess()),
    calcul_(calcul),
    pos_(shrunked ? Eigen::Index(pos) : Eigen::Index(calcul->getRootArrayPosition(pos))),
    vPriorBranch_(vPrior)
  {
    init();
    // Continuous rates not possible for this, since there is no a posteriori for all rates.
    continuousRates_ = false;
  }

  GivenDataSubstitutionProcessSiteSimulator(const GivenDataSubstitutionProcessSiteSimulator& nhss) :
    SimpleSubstitutionProcessSiteSimulator(nhss),
    calcul_(nhss.calcul_),
    pos_(nhss.pos_),
    vPriorBranch_(nhss.vPriorBranch_)
  {}

  GivenDataSubstitutionProcessSiteSimulator& operator=(const GivenDataSubstitutionProcessSiteSimulator& nhss)
  {
    SimpleSubstitutionProcessSiteSimulator::operator=(nhss);
    calcul_ = nhss.calcul_;
    pos_ = nhss.pos_;
    vPriorBranch_ = nhss.vPriorBranch_;

    return *this;
  }

  GivenDataSubstitutionProcessSiteSimulator* clone() const override
  {
    return new GivenDataSubstitutionProcessSiteSimulator(*this);
  }

private:
  /**
   * @brief Init all probabilities.
   *
   * Method called by constructors.
   */
  void init() override;
};
} // end of namespace bpp.
#endif // BPP_PHYL_SIMULATION_GIVENDATASUBSTITUTIONPROCESSSITESIMULATOR_H
