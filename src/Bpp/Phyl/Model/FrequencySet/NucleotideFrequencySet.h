//
// File: NucleotideFrequencySet.h
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

#ifndef _NUCLEOTIDEFREQUENCYSET_H_
#define _NUCLEOTIDEFREQUENCYSET_H_

#include "FrequencySet.h"
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>

namespace bpp
{
/**
 * @brief Parametrize a set of state frequencies for nucleotides.
 */
class NucleotideFrequencySet :
  public virtual FrequencySet
{
public:
  NucleotideFrequencySet* clone() const = 0;

  const NucleicAlphabet* getAlphabet() const = 0;
};

/**
 * @brief Nucleotide FrequencySet using only one parameter, the GC content.
 */
class GCFrequencySet :
  public virtual NucleotideFrequencySet,
  public AbstractFrequencySet
{
public:
  GCFrequencySet(const NucleicAlphabet* alphabet) :
    AbstractFrequencySet(std::shared_ptr<const StateMap>(new CanonicalStateMap(alphabet, false)), "GC.", "GC")
  {
    addParameter_(new Parameter("GC.theta", 0.5, Parameter::PROP_CONSTRAINT_IN));
    getFreq_(0) = getFreq_(1) = getFreq_(2) = getFreq_(3) = 0.25;
  }

  GCFrequencySet(const NucleicAlphabet* alphabet, double theta) :
    AbstractFrequencySet(std::shared_ptr<const StateMap>(new CanonicalStateMap(alphabet, false)), "GC.", "GC")
  {
    addParameter_(new Parameter("GC.theta", theta, Parameter::PROP_CONSTRAINT_IN));
    getFreq_(0) = getFreq_(3) = (1. - theta) / 2.;
    getFreq_(1) = getFreq_(2) = theta / 2.;
  }

  GCFrequencySet* clone() const
  {
    return new GCFrequencySet(*this);
  }

  GCFrequencySet(const GCFrequencySet& gcf) :
    AbstractFrequencySet(gcf)
  {}

public:
  const NucleicAlphabet* getAlphabet() const
  {
    return dynamic_cast<const NucleicAlphabet*>(AbstractFrequencySet::getAlphabet());
  }

  void setFrequencies(const std::vector<double>& frequencies);

protected:
  void fireParameterChanged(const ParameterList& parameters);
};

/**
 * @brief Nucleotide FrequencySet using three independent parameters
 * (theta, theta1, theta2) to modelize the four frequencies:
 *
 * \f[
 * \begin{cases}
 * \theta = \pi_C + \pi_G\\
 * \theta_1 = \frac{\pi_A}{1 - \theta} = \frac{\pi_A}{\pi_A + \pi_T}\\
 * \theta_2 = \frac{\pi_G}{\theta} = \frac{\pi_G}{\pi_C + \pi_G}\\
 * \end{cases}
 * \Longleftrightarrow
 * \begin{cases}
 * \pi_A = \theta_1 (1 - \theta)\\
 * \pi_C = (1 - \theta_2) \theta\\
 * \pi_G = \theta_2 \theta\\
 * \pi_T = (1 - \theta_1)(1 - \theta).
 * \end{cases}
 * \f]
 *
 * with \f$\pi_x\f$ the frequency of nucleotide \f$x\f$.
 */
class FullNucleotideFrequencySet :
  public virtual NucleotideFrequencySet,
  public AbstractFrequencySet
{
public:
  FullNucleotideFrequencySet(const NucleicAlphabet* alphabet, bool allowNullFreqs = false, const std::string& name = "Full");

  FullNucleotideFrequencySet(const NucleicAlphabet* alphabet, double theta, double theta1, double theta2, bool allowNullFreqs = false, const std::string& name = "Full");

  FullNucleotideFrequencySet* clone() const { return new FullNucleotideFrequencySet(*this); }

public:
  const NucleicAlphabet* getAlphabet() const
  {
    return dynamic_cast<const NucleicAlphabet*>(AbstractFrequencySet::getAlphabet());
  }

  void setFrequencies(const std::vector<double>& frequencies);

protected:
  void fireParameterChanged(const ParameterList& parameters);
};


/**
 * @brief FrequencySet useful for homogeneous and stationary models, nucleotide implementation
 *
 * This set contains no parameter.
 */
class FixedNucleotideFrequencySet :
  public virtual NucleotideFrequencySet,
  public FixedFrequencySet
{
public:
  FixedNucleotideFrequencySet(const NucleicAlphabet* alphabet, const std::vector<double>& initFreqs, const std::string& name = "Fixed") :
    FixedFrequencySet(std::shared_ptr<const StateMap>(new CanonicalStateMap(alphabet, false)), initFreqs, name) {}

  /**
   * @brief Construction with uniform frequencies on the letters of
   * the alphabet.
   */
  FixedNucleotideFrequencySet(const NucleicAlphabet* alphabet, const std::string& name = "Fixed") :
    FixedFrequencySet(std::shared_ptr<const StateMap>(new CanonicalStateMap(alphabet, false)), name) {}

  FixedNucleotideFrequencySet* clone() const { return new FixedNucleotideFrequencySet(*this); }

  const NucleicAlphabet* getAlphabet() const
  {
    return dynamic_cast<const NucleicAlphabet*>(AbstractFrequencySet::getAlphabet());
  }
};

/**
 * @brief FrequencySet useful for homogeneous and stationary models, nucleotide implementation
 *
 * This set contains no parameter.
 */
class UserNucleotideFrequencySet :
  public virtual NucleotideFrequencySet,
  public UserFrequencySet
{
public:
  UserNucleotideFrequencySet(const NucleicAlphabet* alphabet, const std::string& path, size_t nCol = 1) :
    UserFrequencySet(std::shared_ptr<const StateMap>(new CanonicalStateMap(alphabet, false)), path, nCol) {}

  UserNucleotideFrequencySet* clone() const { return new UserNucleotideFrequencySet(*this); }

  const NucleicAlphabet* getAlphabet() const
  {
    return dynamic_cast<const NucleicAlphabet*>(AbstractFrequencySet::getAlphabet());
  }
};
} // end of namespace bpp.

#endif// _NUCLEOTIDEFREQUENCYSET_H_
