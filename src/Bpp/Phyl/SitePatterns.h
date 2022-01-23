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
#include <Bpp/Seq/Site.h>
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
    const CruxSymbolListSite* siteP;
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
  std::vector<const CruxSymbolListSite* > sites_;
  std::vector<unsigned int> weights_;
  IndicesType indices_;
  const Alphabet* alpha_;
  bool own_;

public:
  /**
   * @brief Build a new SitePattern object.
   *
   * Look for patterns (unique sites) within a site container.
   *
   * @param sequences The container to look in.
   * @param names The vector of the names of the sequences that are
   * effectively used in the computation. If not empty, the sites are
   * filtered with the sequences names and belong to the sequences will
   * be deleted together with this instance.
   */
  
  SitePatterns(const AlignedValuesContainer* sequences, std::vector<std::string> names= {});

  /**
   * @brief Build a new SitePattern object.
   *
   * Look for patterns (unique sites) within a site container.
   *
   * @param sequences The container to look in.   
   * @param own if the SitePatterns will own the sites of the
   * sequences. In which case the sites will be deleted together with
   * this instance.
   */
  
  SitePatterns(const AlignedValuesContainer* sequences, bool own);

private:
  void init_(const AlignedValuesContainer* sequences, std::vector<std::string> names= {});

  
public:
  
  ~SitePatterns()
  {
    if (own_)
      for (auto si : sites_)
      {
        delete si;
      }
  }

  SitePatterns(const SitePatterns& patterns) :
    names_(patterns.names_),
    sites_(),
    weights_(patterns.weights_),
    indices_(patterns.indices_),
    alpha_(patterns.alpha_),
    own_(patterns.own_)
  {
    if (!patterns.own_)
      sites_ = patterns.sites_;
    else
    {
      for (auto sit : patterns.sites_)
      {
        sites_.push_back(sit->clone());
      }
    }
  }

  SitePatterns& operator=(const SitePatterns& patterns)
  {
    names_     = patterns.names_;
    weights_   = patterns.weights_;
    indices_   = patterns.indices_;
    if (own_)
      for (auto si : sites_)
      {
        delete si;
      }

    sites_.clear();

    if (!patterns.own_)
      sites_ = patterns.sites_;
    else
    {
      for (auto si : patterns.sites_)
      {
        sites_.push_back(si->clone());
      }
    }

    alpha_     = patterns.alpha_;
    own_       = patterns.own_;
    return *this;
  }

  SitePatterns* clone() const { return new SitePatterns(*this); }

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
  AlignedValuesContainer* getSites() const;
};
} // end of namespace bpp.
#endif // BPP_PHYL_SITEPATTERNS_H
