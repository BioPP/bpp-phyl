//
// File: WordFrequencySet.h
// Created by: Laurent Gueguen
// Created on: lundi 2 avril 2012, à 13h 59
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

#ifndef _WORDFREQUENCYSET_H_
#define _WORDFREQUENCYSET_H_

#include <Bpp/Seq/Alphabet/WordAlphabet.h>
#include <Bpp/Seq/Alphabet/CodonAlphabet.h>
#include "FrequencySet.h"

namespace bpp
{

  /*********************************************************************/
/****   Frequencies Set in Words *****/
/*********************************************************************/


/**
 * @brief Frequencies in words computed from the  frequencies on
 * letters. The parameters are the parameters of the Frequencies on
 * letters.
 * The WordFrequencySet owns the std::shared_ptr<FrequencySet> it is built on.
 * Interface class.
 * @author Laurent Guéguen
 */

class WordFrequencySet :
  public virtual FrequencySet
{
protected:
  
  virtual size_t getSizeFromVector(const std::vector<std::shared_ptr<FrequencySet>>& freqVector) = 0;
  
public:
  WordFrequencySet* clone() const = 0;

  virtual const CoreWordAlphabet* getWordAlphabet() const = 0;

  /**
   *@ brief Returns the n-th FrequencySet&
   */
  virtual const std::shared_ptr<FrequencySet> getFrequencySetForLetter(size_t i) const = 0;

  /**
   *@ brief Returns the length of the words
   */
  virtual size_t getLength() const = 0;
};


class AbstractWordFrequencySet :
  public virtual WordFrequencySet,
  public AbstractFrequencySet
{
protected:
  size_t getSizeFromVector(const std::vector<std::shared_ptr<FrequencySet>>& freqVector);
  
public:
  AbstractWordFrequencySet(std::shared_ptr<const StateMap> stateMap, const std::string& prefix = "", const std::string& name="");

  AbstractWordFrequencySet* clone() const = 0;

  AbstractWordFrequencySet(const AbstractWordFrequencySet& af) :
    AbstractFrequencySet(af) {}

  AbstractWordFrequencySet & operator=(const AbstractWordFrequencySet& af)
  {
    AbstractFrequencySet::operator=(af);
    return *this;
  }

  const CoreWordAlphabet* getWordAlphabet() const
  {
    return dynamic_cast<const CoreWordAlphabet*>(AbstractFrequencySet::getAlphabet());
  }

  virtual ~AbstractWordFrequencySet();
  
  /**
   *@ brief Return the length of the words
   */
  size_t getLength() const;
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
  std::vector<std::shared_ptr<FrequencySet>> vFreq_;
  std::vector<std::string> vNestedPrefix_;

public:
  /**
   * @brief Constructor from a WordAlphabet* and a vector of different std::shared_ptr<FrequencySet>.
   * Throws an Exception if their lengths do not match.
   */
  WordFromIndependentFrequencySet(const WordAlphabet* pWA, const std::vector<std::shared_ptr<FrequencySet>>& freqVector, const std::string& prefix = "", const std::string& name="WordFromIndependent");

  WordFromIndependentFrequencySet(const CodonAlphabet* pWA, const std::vector<std::shared_ptr<FrequencySet>>& freqVector, const std::string& prefix = "", const std::string& name="WordFromIndependent");

  WordFromIndependentFrequencySet(const WordFromIndependentFrequencySet& iwfs);

  ~WordFromIndependentFrequencySet();

  WordFromIndependentFrequencySet& operator=(const WordFromIndependentFrequencySet& iwfs);

  WordFromIndependentFrequencySet* clone() const { return new WordFromIndependentFrequencySet(*this); }

public:
  void fireParameterChanged(const ParameterList& pl);

  virtual void updateFrequencies();

  /**
   *@ brief Independent letter frequencies from given word frequencies.
   * The frequencies of a letter at a position is the sum of the
   *    frequencies of the words that have this letter at this
   *    position.
   */
  virtual void setFrequencies(const std::vector<double>& frequencies);

  /**
   *@ brief Return the n-th FrequencySet&
   **/
  const std::shared_ptr<FrequencySet> getFrequencySetForLetter(size_t i) const { return vFreq_[i]; }

  /**
   *@ brief Return the length of the words
   **/

  virtual size_t getLength() const;

  void setNamespace(const std::string& prefix);

  std::string getDescription() const;
};

class WordFromUniqueFrequencySet :
  public AbstractWordFrequencySet
{
protected:
  std::shared_ptr<FrequencySet> pFreq_;
  std::string NestedPrefix_;
  size_t length_;

public:
  /**
   * @brief Constructor from a WordAlphabet* and a std::shared_ptr<FrequencySet>
   *  repeated as many times as the length of the words.
   */
  WordFromUniqueFrequencySet(const WordAlphabet* pWA, std::shared_ptr<FrequencySet> pabsfreq, const std::string& prefix = "", const std::string& name = "WordFromUnique");

  WordFromUniqueFrequencySet(const CodonAlphabet* pWA, std::shared_ptr<FrequencySet> pabsfreq, const std::string& prefix = "", const std::string& name = "WordFromUnique");

  WordFromUniqueFrequencySet(const WordFromUniqueFrequencySet& iwfs);

  WordFromUniqueFrequencySet& operator=(const WordFromUniqueFrequencySet& iwfs);

  ~WordFromUniqueFrequencySet();

  WordFromUniqueFrequencySet* clone() const { return new WordFromUniqueFrequencySet(*this); }

public:
  virtual void fireParameterChanged(const ParameterList& pl);

  /**
   *@ brief letter frequencies from given word frequencies. The
   * frequencies of a letter at a position is the sum of the
   * frequencies of the words that have this letter at this position.
   * The frequencies of each letter is the average of the frequencies
   * of that letter at all positions.
   */
  virtual void setFrequencies(const std::vector<double>& frequencies);

  virtual void updateFrequencies();

  /**
   *@ brief Return the n-th FrequencySet&
   **/
  const std::shared_ptr<FrequencySet> getFrequencySetForLetter(size_t i) const { return pFreq_; }

  size_t getLength() const { return length_; }

  void setNamespace(const std::string& prefix);

  std::string getDescription() const;
};

} // end of namespace bpp.

#endif // _WORDFREQUENCYSET_H_


