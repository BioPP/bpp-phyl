//
// File: StateMap.h
// Authors:
//   Julien Dutheil
// Created: 2012-06-13 15:03:00
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

#ifndef BPP_PHYL_MODEL_STATEMAP_H
#define BPP_PHYL_MODEL_STATEMAP_H

#include <Bpp/Clonable.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>


// From the STL:
#include <vector>
#include <string>

namespace bpp
{
/**
 * @brief Map the states of a given alphabet which have a model state.
 */
class StateMapInterface :
  public virtual Clonable
{
public:
  virtual ~StateMapInterface() {}
  virtual StateMapInterface* clone() const override = 0;

public:
  /**
   * @return The associated alphabet.
   */
  virtual const Alphabet& alphabet() const = 0;

  /**
   * @return A shared_ptr towards the associated alphabet.
   */
  virtual std::shared_ptr<const Alphabet> getAlphabet() const = 0;

  /**
   * @return The number of states supported by the model.
   */
  virtual size_t getNumberOfModelStates() const = 0;

  /**
   * @return A string describing the model state.
   * @param index The state index.
   */
  virtual std::string getStateDescription(size_t index) const = 0;

  /**
   * @return A vector with the corresponding alphabet states for each model state.
   * the size of the vector is the number of model states, not the number of supported alphabet states,
   * as distinct model states can correspond to a single alphabet state.
   */
  virtual const std::vector<int>& getAlphabetStates() const = 0;

  /**
   * @param index The model state.
   * @return The corresponding alphabet state as character code.
   */
  virtual std::string getAlphabetStateAsChar(size_t index) const = 0;

  /**
   * @param index The model state.
   * @return The corresponding alphabet state as int code.
   */
  virtual int getAlphabetStateAsInt(size_t index) const = 0;

  /**
   * @param code The character code of the alphabet state to check.
   * @return The corresponding model states, is any.
   */
  virtual std::vector<size_t> getModelStates(const std::string& code) const = 0;

  /**
   * @param code The int code of the alphabet state to check.
   * @return The corresponding model states, is any.
   */
  virtual std::vector<size_t> getModelStates(int code) const = 0;

  /**
   * @brief Comparison operator==
   *
   * Two StateMaps are identical if they share the same alphabet and states.
   */
  virtual bool operator==(const StateMapInterface& sm) const = 0;

  /**
   * @brief Comparison operator!=
   *
   * Two StateMaps are not identical if they do not share the same alphabet or states.
   */
  virtual bool operator!=(const StateMapInterface& sm) const = 0;
};

/**
 * @brief A convenience partial implementation of the StateMap interface.
 *
 * Model states are stored as their corresponding int codes, stored in a vector 'states_'.
 * This vector has to be initialized and filled by the derived class.
 */
class AbstractStateMap :
  public virtual StateMapInterface
{
protected:
  std::shared_ptr<const Alphabet> alphabet_;
  std::vector<int> states_;

public:
  AbstractStateMap(std::shared_ptr<const Alphabet> alphabet) :
    alphabet_(alphabet),
    states_()
  {}

  AbstractStateMap(const AbstractStateMap& absm) :
    alphabet_(absm.alphabet_),
    states_(absm.states_)
  {}

  AbstractStateMap& operator=(const AbstractStateMap& absm)
  {
    alphabet_ = absm.alphabet_;
    states_ = absm.states_;
    return *this;
  }

public:
  const Alphabet& alphabet() const override { return *alphabet_; }
  std::shared_ptr<const Alphabet> getAlphabet() const override { return alphabet_; }
  size_t getNumberOfModelStates() const override { return states_.size(); }
  const std::vector<int>& getAlphabetStates() const override { return states_; }
  int getAlphabetStateAsInt(size_t index) const override { return states_[index]; }
  std::string getAlphabetStateAsChar(size_t index) const override { return alphabet_->intToChar(states_[index]); }
  std::vector<size_t> getModelStates(int code) const override
  {
    return VectorTools::whichAll(states_, code);
  }
  std::vector<size_t> getModelStates(const std::string& code) const override
  {
    return VectorTools::whichAll(states_, alphabet_->charToInt(code));
  }
  
  bool operator==(const StateMapInterface& sm) const override
  {
    if (alphabet_->getAlphabetType() != sm.alphabet().getAlphabetType()) return false;
    if (states_ != sm.getAlphabetStates()) return false;
    return true;
  }

  bool operator!=(const StateMapInterface& sm) const override
  {
    if (alphabet_->getAlphabetType() != sm.alphabet().getAlphabetType()) return true;
    if (states_ != sm.getAlphabetStates()) return true;
    return false;
  }
};

/**
 * @brief This class implements a state map where all resolved states are modeled.
 *
 * For nucleotides, the underlying states are for instance: A (0), C (1), G (2), T/U (3).
 * Optionally, gaps can be modeled.
 */
class CanonicalStateMap :
  public AbstractStateMap
{
public:

  CanonicalStateMap(std::shared_ptr<const Alphabet> alphabet, bool includeGaps);

  /**
   * @brief this contructors takes an existing StateMap and adds one model states for gaps.
   * If the original StateMap alread had a state for gaps, a new one will be appended.
   */
  CanonicalStateMap(const StateMapInterface& sm, bool includeGaps);
    
  virtual CanonicalStateMap* clone() const { return new CanonicalStateMap(*this); }

  virtual std::string getStateDescription(size_t index) const { return getAlphabetStateAsChar(index); }
};


/**
 * @brief This class implements a state map for Markov modulated models.
 *
 * For nucleotides with two classes, the underlying states are for instance:
 * A (0), C (1), G (2), T/U (3), A (4), C (5), G (6), T/U (7).
 */
class MarkovModulatedStateMap :
  public AbstractStateMap
{
private:
  unsigned int nbClasses_;

public:
  MarkovModulatedStateMap(const StateMapInterface& unitMap, unsigned int nbClasses);
  virtual MarkovModulatedStateMap* clone() const { return new MarkovModulatedStateMap(*this); }

  virtual std::string getStateDescription(size_t index) const { return getAlphabetStateAsChar(index) + TextTools::toString(index % nbClasses_); }
};
}// end of namespace bpp
#endif // BPP_PHYL_MODEL_STATEMAP_H
