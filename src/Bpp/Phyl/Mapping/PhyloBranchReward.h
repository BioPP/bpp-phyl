//
// File: PhyloBranchReward.h
// Created by: Laurent Guéguen
// Created on: dimanche 8 octobre 2017, à 22h 14
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

#ifndef _PHYLOBRANCH_REWARD_H_
#define _PHYLOBRANCH_REWARD_H_

#include <Bpp/Numeric/Number.h>
#include <Bpp/Clonable.h>
#include <Bpp/Exceptions.h>

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
    
    PhyloBranchReward():
      PhyloBranch(),
      rewards_()
    {
    }

    PhyloBranchReward(double length):
      PhyloBranch(length),
      rewards_()
    {
    }

    PhyloBranchReward(const PhyloBranch& branch):
      PhyloBranch(branch),
      rewards_()
    {
    }
    
    /**
     * @brief Copy constructor.
     *
     * @param branch The branch to copy.
     */
    
    PhyloBranchReward(const PhyloBranchReward& branch):
      PhyloBranch(branch),
      rewards_(branch.rewards_)
    {
    }
    
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
    {
    }

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
      if (site>=getNumberOfSites())
        throw BadSizeException("PhyloBranchReward::getSiteReward : bad site number",site,getNumberOfSites());
      return rewards_[site];
    }

    double setSiteReward(size_t site, double rew)
    {
      if (site>=getNumberOfSites())
        throw BadSizeException("PhyloBranchReward::setSiteReward : bad site number",site,getNumberOfSites());
      rewards_[site]=rew;
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
    
    Vdouble getRewards() const
    {
      return rewards_;
    }

    Vdouble& getRewards()
    {
      return rewards_;
    }


  }; 


} //end of namespace bpp.

#endif  //_PHYLOBRANCH_REWARD_H_
