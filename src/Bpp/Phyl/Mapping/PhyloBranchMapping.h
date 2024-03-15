// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MAPPING_PHYLOBRANCHMAPPING_H
#define BPP_PHYL_MAPPING_PHYLOBRANCHMAPPING_H

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
 * available in ProbabilisticSubstitutionMapping class.
 *
 */

class PhyloBranchMapping :
  public PhyloBranch
{
protected:
  /*
   * @brief counts are stored by site / type
   *
   */

  VVdouble counts_;

public:
  /**
   * @brief Constructors.
   *
   * @warning phyloTree_ does not know the edge exists.
   *
   */

  PhyloBranchMapping() :
    PhyloBranch(),
    counts_()
  {}

  PhyloBranchMapping(double length) :
    PhyloBranch(length),
    counts_()
  {}

  PhyloBranchMapping(const PhyloBranch& branch) :
    PhyloBranch(branch),
    counts_()
  {}

  /**
   * @brief Copy constructor.
   *
   * @param branch The branch to copy.
   */

  PhyloBranchMapping(const PhyloBranchMapping& branch) :
    PhyloBranch(branch),
    counts_(branch.counts_)
  {}

  /**
   * @brief Assignation operator.
   *
   * @param branch the branch to copy.
   * @return A reference toward this branch.
   */
  PhyloBranchMapping& operator=(const PhyloBranchMapping& branch)
  {
    PhyloBranch::operator=(branch);
    counts_ = branch.counts_;
    return *this;
  }

  PhyloBranchMapping* clone() const { return new PhyloBranchMapping(*this); }

  /**
   * @brief destructor. In Graph, nothing is changed.
   *
   */
  ~PhyloBranchMapping()
  {}

  /**
   * @brief Sets a number of sites. If the number of types is
   * already defined, it is kept.
   */
  void setNumberOfSites(size_t nbSites)
  {
    VectorTools::resize2(counts_, nbSites, (counts_.size() != 0 ? counts_[0].size() : 0));
  }

  /**
   * @brief Define a number of types.
   */
  void setNumberOfTypes(size_t nbTypes)
  {
    VectorTools::resize2(counts_, counts_.size(), nbTypes);
  }

  /**
   * @brief Define a number of types.
   */
  void setNumberOfSitesAndTypes(size_t nbSites, size_t nbTypes)
  {
    VectorTools::resize2(counts_, nbSites, nbTypes);
  }


  /**
   * @brief Gets the number of sites.
   */
  size_t getNumberOfSites() const
  {
    return counts_.size();
  }

  /**
   * @brief Gets the number of types.
   */
  size_t getNumberOfTypes() const
  {
    return counts_.size() ? counts_[0].size() : 0;
  }

  /**
   * @brief Gets the counts at a given site
   *
   */
  Vdouble& getSiteCount(size_t site)
  {
    return counts_[site];
  }

  const Vdouble& getSiteCount(size_t site) const
  {
    return counts_[site];
  }

  /**
   * @brief Gets the counts at a given site on a given type
   *
   */

  /**
   * @brief With check
   *
   */
  double getSiteTypeCount(size_t site, size_t type) const
  {
    if (site >= getNumberOfSites())
      throw BadSizeException("PhyloBranchMapping::getSiteTypeCount : bad site number", site, getNumberOfSites());
    if (type >= getNumberOfTypes())
      throw BadSizeException("PhyloBranchMapping::getSiteTypeCount : bad site number", type, getNumberOfTypes());
    return counts_[site][type];
  }

  /**
   * @brief Sets the counts at a given site on a given type
   *
   */
  void setSiteTypeCount(size_t site, size_t type, double value)
  {
    if (site >= getNumberOfSites())
      throw BadSizeException("PhyloBranchMapping::setSiteTypeCount : bad site number", site, getNumberOfSites());
    if (type >= getNumberOfTypes())
      throw BadSizeException("PhyloBranchMapping::setSiteTypeCount : bad type number", type, getNumberOfTypes());
    counts_[site][type] = value;
  }


  /**
   * @brief Without check
   *
   */
  double operator()(size_t site, size_t type) const
  {
    return counts_[site][type];
  }

  double& operator()(size_t site, size_t type)
  {
    return counts_[site][type];
  }

  /**
   * @brief return counts
   *
   */
  const VVdouble& getCounts() const
  {
    return counts_;
  }

  VVdouble& getCounts()
  {
    return counts_;
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MAPPING_PHYLOBRANCHMAPPING_H
