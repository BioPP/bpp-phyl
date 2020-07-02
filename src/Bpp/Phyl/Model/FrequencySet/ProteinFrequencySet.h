//
// File: ProteinFrequencySet.h
// Created by: Bastien Boussau
//             Julien Dutheil
// Created on: Tue Aug 21 2007
//

/*
  Copyright or (c) or Copr. Bio++ Development Team, (November 16, 2004)

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

#ifndef _PROTEINFREQUENCYSET_H_
#define _PROTEINFREQUENCYSET_H_

#include <Bpp/Seq/Alphabet/ProteicAlphabet.h>
#include "FrequencySet.h"

namespace bpp
{
/**
 * @brief Parametrize a set of state frequencies for proteins.
 */
  class ProteinFrequencySet :
    public virtual FrequencySet
  {
  public:
  
    ProteinFrequencySet* clone() const = 0;

    const ProteicAlphabet* getAlphabet() const = 0;
  };

/**
 * @brief Protein FrequencySet using 19 independent parameters to
 * model the 20 frequencies.
 *
 * The parameters are called @f$ \theta_{i \in 1..19} @f$, and are
 * initialized so that all frequencies are equal to 0.005. The
 * parametrization depends on the method used. Default
 * method is 1 (ie global ratio).
 *
 * @see Simplex
 */
  class FullProteinFrequencySet :
    public virtual ProteinFrequencySet,
    public FullFrequencySet
  {
  public:
    FullProteinFrequencySet(const ProteicAlphabet* alphabet, bool allowNullFreqs = false, unsigned short method = 1, const std::string& name = "Full") :
      FullFrequencySet(std::shared_ptr<const StateMap>(new CanonicalStateMap(alphabet, false)), allowNullFreqs, method, name) {}
    FullProteinFrequencySet(const ProteicAlphabet* alphabet, const std::vector<double>& initFreqs, bool allowNullFreqs = false, unsigned short method = 1, const std::string& name = "Full") :
      FullFrequencySet(std::shared_ptr<const StateMap>(new CanonicalStateMap(alphabet, false)), initFreqs, allowNullFreqs, method, name) {}

    FullProteinFrequencySet* clone() const { return new FullProteinFrequencySet(*this); }

  public:
    const ProteicAlphabet* getAlphabet() const
    {
      return dynamic_cast<const ProteicAlphabet*>(AbstractFrequencySet::getAlphabet());
    }
  };

/**
 * @brief FrequencySet useful for homogeneous and stationary models, protein implementation
 *
 * This set contains no parameter.
 */
  class FixedProteinFrequencySet :
    public virtual ProteinFrequencySet,
    public FixedFrequencySet
  {
  public:
    FixedProteinFrequencySet(const ProteicAlphabet* alphabet, const std::vector<double>& initFreqs, const std::string& name = "Fixed") :
      FixedFrequencySet(std::shared_ptr<const StateMap>(new CanonicalStateMap(alphabet, false)), initFreqs, name) {}

    /**
     * @brief Construction with uniform frequencies on the letters of
     * the alphabet.
     */
    FixedProteinFrequencySet(const ProteicAlphabet* alphabet, const std::string& name = "Fixed") :
      FixedFrequencySet(std::shared_ptr<const StateMap>(new CanonicalStateMap(alphabet, false)), name) {}

    FixedProteinFrequencySet* clone() const { return new FixedProteinFrequencySet(*this); }

    const ProteicAlphabet* getAlphabet() const
    {
      return dynamic_cast<const ProteicAlphabet*>(AbstractFrequencySet::getAlphabet());
    }
  };

  /**
   * @brief FrequencySet from file
   *
   * This set contains no parameter.
   */
  
  class UserProteinFrequencySet :
    public virtual ProteinFrequencySet,
    public UserFrequencySet
  {
  public:
    UserProteinFrequencySet(const ProteicAlphabet* alphabet, const std::string& path, size_t nCol=1) :
      UserFrequencySet(std::shared_ptr<const StateMap>(new CanonicalStateMap(alphabet, false)), path, nCol) {}
    
    UserProteinFrequencySet* clone() const { return new UserProteinFrequencySet(*this); }

    const ProteicAlphabet* getAlphabet() const
    {
      return dynamic_cast<const ProteicAlphabet*>(AbstractFrequencySet::getAlphabet());
    }
  };


} // end of namespace bpp.

#endif // _PROTEINFREQUENCYSET_H_


