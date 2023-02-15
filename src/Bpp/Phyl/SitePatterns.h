//
// File: SitePatterns.h
// Authors:
//   Julien Dutheil
// Created: 2005-11-29 15:37:00
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
  std::vector< std::shared_ptr<const CoreSiteInterface> > sites_;
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
      std::vector<std::string> names= {});

  
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
