// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_SITEPATTERNS_H
#define BPP_PHYL_SITEPATTERNS_H

#include <Bpp/Clonable.h>
#include <Bpp/Numeric/VectorTools.h>


// From bpp-seq:
#include <Bpp/Seq/CoreSite.h>
#include <Bpp/Seq/Container/SiteContainer.h>

// From the STL:
#include <map>
#include <vector>
#include <string>

#include <Eigen/Core>

namespace bpp
{
/**
 * @brief Data structure for site patterns.
 *
 * 'names' are the sequence names
 * 'sites' points toward a unique site
 * 'weights' is the number of sites identical to this sites
 * 'indices' are the positions in the original container
 */
class SitePatterns :
  public virtual Clonable
{
public:
  typedef Eigen::Matrix<size_t, 1, -1> IndicesType;

private:
  /**
   * @brief Class used for site pattern sorting.
   */
  class SortableSite
  {
public:
    std::string siteS;
    const CoreSiteInterface* siteP;
    size_t originalPosition;

public:
    SortableSite() : siteS(), siteP(0), originalPosition(0) {}
    SortableSite(const SortableSite& ss) : siteS(ss.siteS), siteP(ss.siteP), originalPosition(ss.originalPosition) {}
    SortableSite& operator=(const SortableSite& ss)
    {
      siteS = ss.siteS;
      siteP = ss.siteP;
      originalPosition = ss.originalPosition;
      return *this;
    }

    bool operator<(const SortableSite& ss) const { return siteS < ss.siteS; }

    virtual ~SortableSite() {}
  };

private:
  std::vector<std::string> names_;
  std::vector<std::shared_ptr<const CoreSiteInterface>> sites_;
  std::vector<unsigned int> weights_;
  IndicesType indices_;
  std::shared_ptr<const Alphabet> alpha_;

public:
  /**
   * @brief Build a new SitePattern object.
   *
   * Look for patterns (unique sites) within a site container.
   *
   * @param sequences The container to look in.
   * @param names The vector of the names of the sequences that are
   * effectively used in the computation. If not empty, the sites are
   * filtered with the sequences names.
   */
  SitePatterns(
      const AlignmentDataInterface& sequences,
      std::vector<std::string> names = {});

  SitePatterns* clone() const { return new SitePatterns(*this); }

private:
  void init_(
      const AlignmentDataInterface& sequences,
      std::vector<std::string> names = {});

public:
  /**
   * @return The number of times each unique site was found.
   */
  const std::vector<unsigned int>& getWeights() const { return weights_; }
  /**
   * @return The position of each unique site.
   */
  const IndicesType& getIndices() const { return indices_; }

  /**
   * @return A new container with each unique site.
   */
  std::unique_ptr<AlignmentDataInterface> getSites() const;
};
} // end of namespace bpp.
#endif // BPP_PHYL_SITEPATTERNS_H
