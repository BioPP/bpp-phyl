// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_FREQUENCYSET_WORDFREQUENCYSET_H
#define BPP_PHYL_MODEL_FREQUENCYSET_WORDFREQUENCYSET_H

#include <Bpp/Seq/Alphabet/CodonAlphabet.h>
#include <Bpp/Seq/Alphabet/WordAlphabet.h>

#include "FrequencySet.h"

namespace bpp
{
/*********************************************************************/
/****   Frequencies Set in Words *****/
/*********************************************************************/


/**
 * @brief Frequencies in words computed from the frequencies on
 * letters. The parameters are the parameters of the Frequencies on
 * letters.
 * The WordFrequencySet owns the std::shared_ptr<FrequencySet> it is built on.
 * Interface class.
 * @author Laurent Guéguen
 */
class WordFrequencySetInterface :
  public virtual FrequencySetInterface
{
protected:
  virtual size_t getSizeFromVector(const std::vector<std::unique_ptr<FrequencySetInterface> >& freqVector) = 0;

public:
  WordFrequencySetInterface* clone() const override = 0;

  virtual std::shared_ptr<const CoreWordAlphabet> getWordAlphabet() const = 0;

  /**
   * @brief Returns the n-th FrequencySet.
   */
  virtual const FrequencySetInterface& frequencySetForLetter(size_t i) const = 0;

  /**
   * @brief Returns the length of the words
   */
  virtual size_t getLength() const = 0;
};


class AbstractWordFrequencySet :
  public virtual WordFrequencySetInterface,
  public AbstractFrequencySet
{
protected:
  size_t getSizeFromVector(const std::vector<std::unique_ptr<FrequencySetInterface> >& freqVector) override;

public:
  AbstractWordFrequencySet(std::shared_ptr<const StateMapInterface> stateMap, const std::string& prefix = "", const std::string& name = "");

  AbstractWordFrequencySet* clone() const override = 0;

  AbstractWordFrequencySet(const AbstractWordFrequencySet& af) :
    AbstractFrequencySet(af) {}

  AbstractWordFrequencySet& operator=(const AbstractWordFrequencySet& af)
  {
    AbstractFrequencySet::operator=(af);
    return *this;
  }

  std::shared_ptr<const CoreWordAlphabet> getWordAlphabet() const override
  {
    return std::dynamic_pointer_cast<const CoreWordAlphabet>(getAlphabet());
  }

  virtual ~AbstractWordFrequencySet();

  /**
   *@ brief Return the length of the words
   */
  size_t getLength() const override;
};


/**
 * @brief the Frequencies in words are the product of Independent Frequencies in letters
 *
 * @author Laurent Guéguen
 */
class WordFromIndependentFrequencySet :
  public AbstractWordFrequencySet
{
protected:
  std::vector<std::unique_ptr<FrequencySetInterface>> vFreq_;
  std::vector<std::string> vNestedPrefix_;

public:
  /**
   * @brief Constructor from a WordAlphabet* and a vector of different std::shared_ptr<FrequencySet>.
   * Throws an Exception if their lengths do not match.
   */
  WordFromIndependentFrequencySet(
      std::shared_ptr<const WordAlphabet> pWA,
      std::vector<std::unique_ptr<FrequencySetInterface>>& freqVector,
      const std::string& prefix = "",
      const std::string& name = "WordFromIndependent");

  WordFromIndependentFrequencySet(
      std::shared_ptr<const CodonAlphabet> pWA,
      std::vector<std::unique_ptr<FrequencySetInterface>>& freqVector,
      const std::string& prefix = "",
      const std::string& name = "WordFromIndependent");

  WordFromIndependentFrequencySet(const WordFromIndependentFrequencySet& iwfs);

  virtual ~WordFromIndependentFrequencySet();

  WordFromIndependentFrequencySet& operator=(const WordFromIndependentFrequencySet& iwfs);

  WordFromIndependentFrequencySet* clone() const override
  {
    return new WordFromIndependentFrequencySet(*this);
  }

public:
  void fireParameterChanged(const ParameterList& pl) override;

  virtual void updateFrequencies();

  /**
   * @brief Independent letter frequencies from given word frequencies.
   * The frequencies of a letter at a position is the sum of the
   *    frequencies of the words that have this letter at this
   *    position.
   */
  virtual void setFrequencies(const std::vector<double>& frequencies) override;

  /**
   * @brief Return the n-th FrequencySet&
   */
  const FrequencySetInterface& frequencySetForLetter(size_t i) const override
  { 
    return *vFreq_[i]; 
  }

  /**
   * @brief Return the length of the words
   */
  virtual size_t getLength() const override;

  void setNamespace(const std::string& prefix) override;

  std::string getDescription() const;
};

class WordFromUniqueFrequencySet :
  public AbstractWordFrequencySet
{

protected:
  std::unique_ptr<FrequencySetInterface> pFreq_;
  std::string NestedPrefix_;
  size_t length_;

public:
  /**
   * @brief Constructor from a WordAlphabet* and a std::shared_ptr<FrequencySet>
   *  repeated as many times as the length of the words.
   */
  WordFromUniqueFrequencySet(
      std::shared_ptr<const WordAlphabet> pWA,
      std::unique_ptr<FrequencySetInterface> pabsfreq,
      const std::string& prefix = "",
      const std::string& name = "WordFromUnique");

  WordFromUniqueFrequencySet(
      std::shared_ptr<const CodonAlphabet> pWA,
      std::unique_ptr<FrequencySetInterface> pabsfreq,
      const std::string& prefix = "",
      const std::string& name = "WordFromUnique");

  WordFromUniqueFrequencySet(const WordFromUniqueFrequencySet& iwfs);

  WordFromUniqueFrequencySet& operator=(const WordFromUniqueFrequencySet& iwfs);

  virtual ~WordFromUniqueFrequencySet();

  WordFromUniqueFrequencySet* clone() const override
  {
    return new WordFromUniqueFrequencySet(*this);
  }

public:
  virtual void fireParameterChanged(const ParameterList& pl) override;

  /**
   * @brief letter frequencies from given word frequencies. The
   * frequencies of a letter at a position is the sum of the
   * frequencies of the words that have this letter at this position.
   * The frequencies of each letter is the average of the frequencies
   * of that letter at all positions.
   */
  virtual void setFrequencies(const std::vector<double>& frequencies) override;

  virtual void updateFrequencies();

  /**
   * @brief Return the n-th FrequencySet&
   */
  const FrequencySetInterface& frequencySetForLetter(size_t i) const override
  { 
    return *pFreq_;
  }

  size_t getLength() const override { return length_; }

  void setNamespace(const std::string& prefix) override;

  std::string getDescription() const;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_FREQUENCYSET_WORDFREQUENCYSET_H
