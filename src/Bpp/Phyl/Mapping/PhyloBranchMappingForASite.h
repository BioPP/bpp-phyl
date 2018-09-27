//
// File: PhyloBranchMappingForASite.h
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

#ifndef _PHYLOBRANCH_MAPPING_FOR_A_SITEH_
#define _PHYLOBRANCH_MAPPING_FOR_A_SITEH_

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
   * available in ProbabilisticSubstitutionMapping class.
   *
   */
  
  class PhyloBranchMappingForASite :
    public PhyloBranch
  {
  protected:
    /*
     * @brief counts are stored by type
     *
     */
    
    Vdouble counts_;
    
  public:
    /**
     * @brief Constructors.
     *
     * @warning phyloTree_ does not know the edge exists.
     *
     */
    
    PhyloBranchMappingForASite():
      PhyloBranch(),
      counts_()
    {
    }

    PhyloBranchMappingForASite(double length):
      PhyloBranch(length),
      counts_()
    {
    }

    PhyloBranchMappingForASite(const PhyloBranch& branch):
      PhyloBranch(branch),
      counts_()
    {
    }
    
    /**
     * @brief Copy constructor.
     *
     * @param branch The branch to copy.
     */
    
    PhyloBranchMappingForASite(const PhyloBranchMappingForASite& branch):
      PhyloBranch(branch),
      counts_(branch.counts_)
    {
    }
    
    /**
     * @brief Assignation operator.
     *
     * @param branch the branch to copy.
     * @return A reference toward this branch.
     */

    PhyloBranchMappingForASite& operator=(const PhyloBranchMappingForASite& branch)
    {
      PhyloBranch::operator=(branch);
      counts_ = branch.counts_;
      return *this;
      
    }
    
    PhyloBranchMappingForASite* clone() const { return new PhyloBranchMappingForASite(*this); }
    
    /**
     * @brief destructor. In Graph, nothing is changed.
     *
     */
    
    ~PhyloBranchMappingForASite()
    {
    }

    /**
     * @brief Define a number of types.
     */
    
    void setNumberOfTypes(size_t nbTypes)
    {
      counts_.resize(nbTypes);
    }

    /**
     * @brief Gets the number of types.
     */
    
    size_t getNumberOfTypes() const
    {
      return counts_.size();
    }
    
    /**
     * @brief Gets the counts
     *
     */
    
    Vdouble& getCounts()
    {
      return counts_;
    }

    const Vdouble& getCounts() const
    {
      return counts_;
    }

    /**
     * @brief Gets the counts on a given type
     *
     */

    /**
     * @brief With check
     *
     */
    
    double getTypeCount(size_t type) const
    {
      if (type>=getNumberOfTypes())
        throw BadSizeException("PhyloBranchMappingForASite::getSiteTypeCount : bad site number",type,getNumberOfTypes());
      return counts_[type];
    }

    /**
     * @brief Sets the counts at a given site on a given type
     *
     */
    
    void setTypeCount(size_t type, double value)
    {
      if (type>=getNumberOfTypes())
        throw BadSizeException("PhyloBranchMappingForASite::setSiteTypeCount : bad type number",type,getNumberOfTypes());
      counts_[type]=value;
    }

    
    /**
     * @brief Without check
     *
     */
    
    double operator()(size_t type) const
    {
      return counts_[type];
    }

    double& operator()(size_t type)
    {
      return counts_[type];
    }

  }; 


} //end of namespace bpp.

#endif  //_PHYLOBRANCH_MAPPING_FOR_A_SITE_H_
