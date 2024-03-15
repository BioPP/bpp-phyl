// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_FREQUENCYSET_FREQUENCYSET_H
#define BPP_PHYL_MODEL_FREQUENCYSET_FREQUENCYSET_H


#include "../StateMap.h"

// From bpp-core:
#include <Bpp/Numeric/ParameterAliasable.h>
#include <Bpp/Numeric/AbstractParameterAliasable.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Prob/Simplex.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>

namespace bpp
{
/**
 * @brief Parametrize a set of state frequencies.
 *
 * Frequencies are ordered according to alphabet states.
 */
class FrequencySetInterface :
  public virtual ParameterAliasable
{
public:
  FrequencySetInterface* clone() const = 0;

public:
  /**
   * @return The alphabet associated to this set.
   */
  virtual std::shared_ptr<const Alphabet> getAlphabet() const = 0;

  /**
   * @return The alphabet associated to this set.
   */
  virtual const Alphabet& alphabet() const = 0;

  /**
   * @return The mapping of model states with alphabet states.
   */
  virtual const StateMapInterface& stateMap() const = 0;

  /**
   * @return A shared_ptr toward the mapping of model states with alphabet states.
   */
  virtual std::shared_ptr<const StateMapInterface> getStateMap() const = 0;

  /**
   * @return The frequencies values of the set.
   */
  virtual const Vdouble& getFrequencies() const = 0;

  /**
   * @return The frequencies of each alphabet states according to this model.
   */
  virtual const std::map<int, double> getAlphabetStatesFrequencies() const = 0;

  /**
   * @brief Set the parameters in order to match a given set of frequencies.
   *
   * @param frequencies The set of frequencies to match.
   * @throw DimensionException If the number of frequencies does not match the size of the alphabet.
   * @throw Exception If the frequencies do not sum to 1.
   */
  virtual void setFrequencies(const std::vector<double>& frequencies) = 0;

  /**
   * @brief Set the Frequencies from the one of the map which keys
   *  match with a letter of the Alphabet.
   *  The frequencies are normalized so that the matching values sum 1.
   *
   * @param frequencies The set of frequencies to match.
   */
  virtual void setFrequenciesFromAlphabetStatesFrequencies(const std::map<int, double>& frequencies) = 0;

  virtual std::string getName() const = 0;

  /**
   * @return The number of frequencies in the set. This is equivalent to getStateMap().getNumberOfModelStates().
   */
  virtual size_t getNumberOfFrequencies() const = 0;

public:
  static std::shared_ptr<IntervalConstraint> FREQUENCE_CONSTRAINT_SMALL;
  static std::shared_ptr<IntervalConstraint> FREQUENCE_CONSTRAINT_MILLI;
  static std::shared_ptr<IntervalConstraint> FREQUENCE_CONSTRAINT_CENTI;
};

/**
 * @brief Basic implementation of the FrequencySet interface.
 */

class AbstractFrequencySet :
  public virtual FrequencySetInterface,
  public AbstractParameterAliasable
{
private:
  std::shared_ptr<const Alphabet> alphabet_;
  std::shared_ptr<const StateMapInterface> stateMap_;
  std::vector<double> freq_;
  std::string name_;

public:
  AbstractFrequencySet(
      std::shared_ptr<const StateMapInterface> stateMap,
      const std::string& prefix,
      const std::string& name) :
    AbstractParameterAliasable(prefix),
    alphabet_(stateMap->getAlphabet()),
    stateMap_(stateMap),
    freq_(stateMap->getNumberOfModelStates()),
    name_(name)
  {}

  AbstractFrequencySet(const AbstractFrequencySet& af) :
    AbstractParameterAliasable(af),
    alphabet_(af.alphabet_),
    stateMap_(af.stateMap_),
    freq_(af.freq_),
    name_(af.name_)
  {}

  AbstractFrequencySet& operator=(const AbstractFrequencySet& af)
  {
    AbstractParameterAliasable::operator=(af);
    alphabet_ = af.alphabet_;
    stateMap_ = af.stateMap_;
    freq_ = af.freq_;
    name_ = af.name_;
    return *this;
  }

public:
  std::shared_ptr<const Alphabet> getAlphabet() const override { return alphabet_; }
  
  const Alphabet& alphabet() const override { return *alphabet_; }

  const StateMapInterface& stateMap() const override { return *stateMap_; }
  
  std::shared_ptr<const StateMapInterface> getStateMap() const override { return stateMap_; }

  const Vdouble& getFrequencies() const override { return freq_; }

  const std::map<int, double> getAlphabetStatesFrequencies() const override;

  /**
   * @brief Set the Frequencies from the one of the map which keys
   *  match with a letter of the Alphabet.
   *  The frequencies are normalized so that the matching values sum 1.
   *
   * In this implementation, all model states with the same alphabet state are given the same frequency.
   *
   * @param frequencies The set of frequencies to match.
   */
  void setFrequenciesFromAlphabetStatesFrequencies(const std::map<int, double>& frequencies) override;

  size_t getNumberOfFrequencies() const override { return freq_.size(); }

  std::string getName() const override { return name_; }

  void normalize()
  {
    double x = 0;
    for (auto f : freq_)
    {
      x += f;
    }
    freq_ /= x;
  }

protected:
  std::vector<double>& getFrequencies_() { return freq_; }
  double& getFreq_(size_t i) { return freq_[i]; }
  const double& getFreq_(size_t i) const { return freq_[i]; }
  void setFrequencies_(const std::vector<double>& frequencies) { freq_ = frequencies; }
};



/**
 * @brief A generic FrequencySet allowing for the estimation of all frequencies.
 *
 * The FrequencySet has hence n-1 parameters, where n is the size of
 * the input alphabet.
 *
 * The parametrization depends on the method used.
 * Default method is 1 (ie global ratio).
 *
 * @see Simplex
 */
class FullFrequencySet :
  public AbstractFrequencySet
{
private:
  /**
   * @brief Simplex to handle the probabilities and the parameters.
   */
  Simplex sFreq_;

public:
  /**
   * @brief Construction with uniform frequencies on the states of
   * the alphabet.
   */
  FullFrequencySet(
      std::shared_ptr<const StateMapInterface> stateMap,
      bool allowNullFreqs = false,
      unsigned short method = 1,
      const std::string& name = "Full");

  FullFrequencySet(
      std::shared_ptr<const StateMapInterface> stateMap,
      const std::vector<double>& initFreqs,
      bool allowNullFreqs = false,
      unsigned short method = 1,
      const std::string& name = "Full");

  FullFrequencySet* clone() const override
  {
    return new FullFrequencySet(*this);
  }

public:
  void setFrequencies(const std::vector<double>& frequencies) override;

  unsigned short getMethod() const { return sFreq_.getMethod();}

  void setNamespace(const std::string& nameSpace) override;

protected:
  void fireParameterChanged(const ParameterList& parameters) override;

private:
  void updateFreq_();
};

class TransitionModelInterface;

/**
 * @brief FrequencySet defined from the equilibrium distribution
 * of a given model.
 *
 * Its parameters are the parameters of the model.
 */
class FromModelFrequencySet :
  public AbstractFrequencySet
{
private:
  std::shared_ptr<TransitionModelInterface> model_;

public:
  FromModelFrequencySet(std::shared_ptr<TransitionModelInterface> model);

  FromModelFrequencySet(const FromModelFrequencySet& fmfs);

  FromModelFrequencySet& operator=(const FromModelFrequencySet& fmfs);

  FromModelFrequencySet* clone() const override { return new FromModelFrequencySet(*this); }

  virtual ~FromModelFrequencySet();

public:
  const TransitionModelInterface& model() const
  {
    return *model_;
  }

  std::shared_ptr<const TransitionModelInterface> getModel() const
  {
    return model_;
  }

  void setFrequencies(const std::vector<double>& frequencies) override;

  void fireParameterChanged(const ParameterList& pl) override;

  void setNamespace(const std::string& name) override;
};


/**
 * @brief FrequencySet to be used with a Markov-modulated substitution model.
 *
 * This implementation uses one parameter per character state frequency.
 * The rate states are assumed to be fixed and are passed as an argument to the constructor, together with a 'regular'
 * FrequencySet. The number of parameters hence do not depends on the number of rates used.
 */
class MarkovModulatedFrequencySet :
  public AbstractFrequencySet
{
private:
  std::unique_ptr<FrequencySetInterface> freqSet_;
  std::vector<double> rateFreqs_;

public:
  MarkovModulatedFrequencySet(
      std::unique_ptr<FrequencySetInterface> freqSet,
      const std::vector<double>& rateFreqs);

  MarkovModulatedFrequencySet(const MarkovModulatedFrequencySet& mmfs) :
    AbstractFrequencySet(mmfs),
    freqSet_(mmfs.freqSet_->clone()),
    rateFreqs_(mmfs.rateFreqs_)
  {}

  MarkovModulatedFrequencySet& operator=(const MarkovModulatedFrequencySet& mmfs)
  {
    AbstractFrequencySet::operator=(mmfs);
    freqSet_.reset(mmfs.freqSet_->clone());
    rateFreqs_ = mmfs.rateFreqs_;
    return *this;
  }

  MarkovModulatedFrequencySet* clone() const override { return new MarkovModulatedFrequencySet(*this); }

  virtual ~MarkovModulatedFrequencySet() {}

public:
  void setFrequencies(const std::vector<double>& frequencies) override
  {
    // Just forward this method to the sequence state frequencies set. This may change in the future...
    freqSet_->setFrequencies(frequencies);
  }

  void fireParameterChanged(const ParameterList& pl) override
  {
    freqSet_->matchParametersValues(pl);
    setFrequencies_(VectorTools::kroneckerMult(rateFreqs_, freqSet_->getFrequencies()));
  }

  const FrequencySetInterface& getStatesFrequencySet() const { return *freqSet_; }
};


/**
 * @brief FrequencySet useful for homogeneous and stationary models.
 *
 * This set contains no parameter.
 */
class FixedFrequencySet :
  public AbstractFrequencySet
{
public:
  /**
   * @brief Construction with user-defined frequencies on the states of the model.
   *
   * @param stateMap The model states for which frequencies should be built.
   * @param initFreqs The frequencies to use. The size of the vector should match the number of model states.
   * @param name The name of the set.
   * @throw Exception In case the number of frequencies does not match the number of model states.
   */
  FixedFrequencySet(
      std::shared_ptr<const StateMapInterface> stateMap,
      const std::vector<double>& initFreqs,
      const std::string& name = "Fixed");

  /**
   * @brief Construction with uniform frequencies on the states of the model.
   *
   * @param stateMap The model states for which frequencies should be built.
   * @param name The name of the set.
   */
  FixedFrequencySet(
      std::shared_ptr<const StateMapInterface> stateMap,
      const std::string& name = "Fixed");

  FixedFrequencySet* clone() const override { return new FixedFrequencySet(*this); }

public:
  void setFrequencies(const std::vector<double>& frequencies) override;
};


/**
 * @brief FrequencySet to be read in a file. More specifically, a
 * frequency set is read in a column of a given file, which column
 * number is given in argument (default: 1).
 */
class UserFrequencySet :
  public AbstractFrequencySet
{
private:
  std::string path_;
  size_t nCol_;

public:
  UserFrequencySet(
      std::shared_ptr<const StateMapInterface> stateMap,
      const std::string& path,
      size_t nCol = 1);

  UserFrequencySet(const UserFrequencySet& fmfs);

  UserFrequencySet& operator=(const UserFrequencySet& fmfs);

  UserFrequencySet* clone() const override { return new UserFrequencySet(*this); }

  virtual ~UserFrequencySet(){}

public:
  const std::string& getPath() const { return path_; }

  size_t getColumnNumber() const
  {
    return nCol_;
  }

  void setFrequencies(const std::vector<double>& frequencies) override;

protected:
  void readFromFile_();
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_FREQUENCYSET_FREQUENCYSET_H
