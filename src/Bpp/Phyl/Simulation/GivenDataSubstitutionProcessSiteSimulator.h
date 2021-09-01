//
// File: GivenDataSubstitutionProcessSiteSimulator.h
// Created by: Laurent Guéguen
// Created on: lundi 12 octobre 2020, à 06h 04
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

#ifndef _GIVEN_DATA_SIMPLE_SUBSTITUTION_PROCESS_SITE_SIMULATOR_H_
#define _GIVEN_DATA_SIMPLE_SUBSTITUTION_PROCESS_SITE_SIMULATOR_H_

#include <Bpp/Phyl/NewLikelihood/DataFlow/LikelihoodCalculationSingleProcess.h>
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
 *
 */

  class GivenDataSubstitutionProcessSiteSimulator:
    public SimpleSubstitutionProcessSiteSimulator
  {
  private:
    std::shared_ptr<LikelihoodCalculationSingleProcess> calcul_;

    /*
     * @brief Position of the copied site, in SHRUNKED data
     *
     */
    
    Eigen::Index pos_;
    
  public:
    /*
     * @brief Build a Site Simulator of histories from the a posteriori likelihoods at a given site
     *
     * @param calcul the a posteriori likelihood calculation
     * @param pos the position of the site to imitate
     * @paream shrunked if the given position is on the shrunked data (default: false)
     */
     
    GivenDataSubstitutionProcessSiteSimulator(std::shared_ptr<LikelihoodCalculationSingleProcess> calcul, size_t pos, bool shrunked = false) : SimpleSubstitutionProcessSiteSimulator(calcul->getSubstitutionProcess()),
                                                                                                                                               calcul_(calcul), pos_(shrunked?Eigen::Index(pos):Eigen::Index(calcul->getRootArrayPosition(pos)))
    {
      init();
      /*
       * Continuous rates not possible for this, since there is no a posteriori for all rates.
       *
       */
      
      continuousRates_ = false;
    }

    GivenDataSubstitutionProcessSiteSimulator(const GivenDataSubstitutionProcessSiteSimulator& nhss) :
      SimpleSubstitutionProcessSiteSimulator(nhss),
      calcul_(nhss.calcul_),
      pos_(nhss.pos_)
    {}

    GivenDataSubstitutionProcessSiteSimulator& operator=(const GivenDataSubstitutionProcessSiteSimulator& nhss)
    {
      SimpleSubstitutionProcessSiteSimulator::operator=(nhss);
      calcul_=nhss.calcul_;
      pos_=nhss.pos_;
      
      return *this;
    }

    GivenDataSubstitutionProcessSiteSimulator* clone() const { return new GivenDataSubstitutionProcessSiteSimulator(*this); }

  private:
    /**
     * @brief Init all probabilities.
     *
     * Method called by constructors.
     */
    void init();
  };

} //end of namespace bpp.

#endif //_GIVEN_DATA_SIMPLE_SUBSTITUTION_PROCESS_SITE_SIMULATOR_H_

