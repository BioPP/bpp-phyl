// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_SITEPARTITION_H
#define BPP_PHYL_LIKELIHOOD_SITEPARTITION_H


namespace bpp
{
/**
 * @brief This is the interface for classes describing a site partition,
 * each partition being intended to have its own substitution model.
 *
 * A site partition defines the distinct patterns to be found in the data,
 * A pattern being a unique combination of site * model.
 *
 * @warning This interface is still under construction. Breaks are expected to occur!
 */
class SitePartition :
  public virtual Clonable
{
public:
  SitePartition* clone() const = 0;

public:
  virtual size_t getNumberOfPartitions() const = 0;
  virtual size_t getNumberOfPatternsForPartition(size_t partitionIndex) const = 0;
};

/**
 * @brief Trivial site partition: all sites belong to the same, unique partition.
 */
class TrivialSitePartition :
  public virtual SitePartition
{
private:
  size_t nbDistinctSites_;

public:
  TrivialSitePartition(size_t nbDistinctSites) : nbDistinctSites_(nbDistinctSites) {}

  TrivialSitePartition* clone() const { return new TrivialSitePartition(*this); }

public:
  size_t getNumberOfPartitions() const { return 1; }
  size_t getNumberOfPatternsForPartition(size_t partitionIndex) const { return nbDistinctSites_; }
};
} // end of namespace bpp
#endif // BPP_PHYL_LIKELIHOOD_SITEPARTITION_H
