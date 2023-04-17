//
// File: SubstitutionDistance.h
// Authors:
//   Laurent GuÃÂ©guen
// Created: lundi 26 mars 2018, ÃÂ  15h 08
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004, 2005, 2006)
  
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

#ifndef BPP_PHYL_MAPPING_SUBSTITUTIONDISTANCE_H
#define BPP_PHYL_MAPPING_SUBSTITUTIONDISTANCE_H


#include "SubstitutionCount.h"

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/AlphabetIndex/AlphabetIndex2.h>

namespace bpp
{
/**
 * @brief Interface allowing for using distances between states in
 * substitution counts.
 *
 * These distances are used for integration of substitution distances
 * before the mapping process (see WeightedSubstitutionCount for
 * weight on the final counts, after the mapping process).
 *
 */

class SubstitutionDistance :
  public virtual SubstitutionCountInterface
{
public:
  virtual void setDistances(std::shared_ptr<const AlphabetIndex2> index) = 0;
  virtual bool hasDistances() const = 0;
  virtual std::shared_ptr<const AlphabetIndex2> getDistances() const = 0;
};

/**
 * @brief Partial implementation of the SubstitutionDistance interface.
 */
class AbstractSubstitutionDistance :
  public virtual SubstitutionDistance
{
protected:
  std::shared_ptr<const AlphabetIndex2> distances_;

public:
  AbstractSubstitutionDistance(std::shared_ptr<const AlphabetIndex2> distances) :
    distances_(distances)
  {}

  AbstractSubstitutionDistance(const AbstractSubstitutionDistance& index) :
    distances_(index.distances_)
  {}

  AbstractSubstitutionDistance& operator=(const AbstractSubstitutionDistance& index)
  {
    distances_ = index.distances_;

    return *this;
  }

  virtual ~AbstractSubstitutionDistance()
  {}

public:
  void setDistances(std::shared_ptr<const AlphabetIndex2> distances);
  bool hasDistances() const { return distances_.get() != 0; }
  std::shared_ptr<const AlphabetIndex2> getDistances() const { return distances_; }

protected:
  virtual void distancesHaveChanged() = 0;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MAPPING_SUBSTITUTIONDISTANCE_H
