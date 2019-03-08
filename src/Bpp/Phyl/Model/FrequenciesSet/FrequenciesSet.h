//
// File: FrequenciesSet.h
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

#ifndef _FREQUENCIESSET_H_
#define _FREQUENCIESSET_H_

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

  class TransitionModel;
  
  class FrequenciesSet :
    public virtual ParameterAliasable
  {
  public:
    
    FrequenciesSet* clone() const = 0;

  public:
    /**
     * @return The alphabet associated to this set.
     */
    virtual const Alphabet* getAlphabet() const = 0;

    /**
     * @return The mapping of model states with alphabet states.
     */
    virtual const StateMap& getStateMap() const = 0;

    /**
     * @return Share the mapping of model states with alphabet states.
     */

    virtual std::shared_ptr<const StateMap> shareStateMap() const = 0;

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
  };

/**
 * @brief Basic implementation of the FrequenciesSet interface.
 */

  class AbstractFrequenciesSet :
    public virtual FrequenciesSet,
    public AbstractParameterAliasable
  {
  private:
    const Alphabet* alphabet_;
    std::shared_ptr<const StateMap> stateMap_;
    std::vector<double> freq_;
    std::string name_;

  public:
    AbstractFrequenciesSet(std::shared_ptr<const StateMap> stateMap, const std::string& prefix, const std::string& name) :
      AbstractParameterAliasable(prefix),
      alphabet_(stateMap->getAlphabet()),
      stateMap_(stateMap),
      freq_(stateMap->getNumberOfModelStates()),
      name_(name)
    {
    }

    AbstractFrequenciesSet* clone() const = 0;

    AbstractFrequenciesSet(const AbstractFrequenciesSet& af) :
      AbstractParameterAliasable(af),
      alphabet_(af.alphabet_),
      stateMap_(af.stateMap_),
      freq_(af.freq_),
      name_(af.name_)
    {}

    AbstractFrequenciesSet& operator=(const AbstractFrequenciesSet& af)
    {
      AbstractParameterAliasable::operator=(af);
      alphabet_ = af.alphabet_;
      stateMap_ = af.stateMap_;
      freq_ = af.freq_;
      name_ = af.name_;
      return *this;
    }

  public:
    const Alphabet* getAlphabet() const { return alphabet_; }

    const StateMap& getStateMap() const { return *stateMap_; }

    std::shared_ptr<const StateMap> shareStateMap() const { return stateMap_; }

    const Vdouble& getFrequencies() const { return freq_; }
  
    const std::map<int, double> getAlphabetStatesFrequencies() const;

    /**
     * @brief Set the Frequencies from the one of the map which keys
     *  match with a letter of the Alphabet.
     *  The frequencies are normalized so that the matching values sum 1.
     *
     * In this implementation, all model states with the same alphabet state are given the same frequency.
     * 
     * @param frequencies The set of frequencies to match.
     */
    void setFrequenciesFromAlphabetStatesFrequencies(const std::map<int, double>& frequencies);

    size_t getNumberOfFrequencies() const { return freq_.size(); }

    std::string getName() const { return(name_); }

    void normalize()
    {
      double x = 0;
      for (size_t i = 0; i < freq_.size(); i++)
        x += freq_[i];
      freq_ /= x;
    }

  protected:
    std::vector<double>& getFrequencies_() { return freq_; }
    double& getFreq_(size_t i) { return freq_[i]; }
    const double& getFreq_(size_t i) const { return freq_[i]; }
    void setFrequencies_(const std::vector<double>& frequencies) { freq_ = frequencies; }
  };

/**
 * @brief A generic FrequenciesSet allowing for the estimation of all frequencies.
 *
 * The FrequenciesSet has hence n-1 parameters, where n is the size of
 * the input alphabet.
 *
 * The parametrization depends on the method used.
 * Default method is 1 (ie global ratio).
 *
 * @see Simplex
 */

  class FullFrequenciesSet :
    public AbstractFrequenciesSet
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
    FullFrequenciesSet(std::shared_ptr<const StateMap> stateMap, bool allowNullFreqs = false, unsigned short method = 1, const std::string& name = "Full.");
    FullFrequenciesSet(std::shared_ptr<const StateMap> stateMap, const std::vector<double>& initFreqs, bool allowNullFreqs = false, unsigned short method = 1, const std::string& name = "Full.");

    FullFrequenciesSet* clone() const { return new FullFrequenciesSet(*this); }

  public:
    void setFrequencies(const std::vector<double>& frequencies);

    unsigned short getMethod() const { return sFreq_.getMethod();}

    void setNamespace(const std::string& nameSpace);
  
  protected:
    void fireParameterChanged(const ParameterList& parameters);

  private:
    void updateFreq_();
  };

  /**
   * @brief FrequenciesSet defined from the equilibrium distribution
   * of a given model.
   *
   * Its parameters are the parameters of the model.
   */
  
  class FromModelFrequenciesSet :
    public AbstractFrequenciesSet
  {
  private:
    TransitionModel* model_;

  public:
    FromModelFrequenciesSet(TransitionModel* model);

    FromModelFrequenciesSet(const FromModelFrequenciesSet& fmfs);
    
    FromModelFrequenciesSet& operator=(const FromModelFrequenciesSet& fmfs);

    FromModelFrequenciesSet* clone() const { return new FromModelFrequenciesSet(*this); }

    ~FromModelFrequenciesSet();

  public:

    const TransitionModel* getModel() const
    {
      return model_;
    }
    
    void setFrequencies(const std::vector<double>& frequencies);

    void fireParameterChanged(const ParameterList& pl);

    void setNamespace(const std::string& name);
    
  };


/**
 * @brief FrequenciesSet to be used with a Markov-modulated substitution model.
 *
 * This implementation uses one parameter per character state frequency.
 * The rate states are assumed to be fixed and are passed as an argument to the constructor, together with a 'regular'
 * FrequenciesSet. The number of parameters hence do not depends on the number of rates used.
 */
  class MarkovModulatedFrequenciesSet :
    public AbstractFrequenciesSet
  {
  private:
    FrequenciesSet* freqSet_;
    std::vector<double> rateFreqs_;

  public:
    MarkovModulatedFrequenciesSet(FrequenciesSet* freqSet, const std::vector<double>& rateFreqs);

    MarkovModulatedFrequenciesSet(const MarkovModulatedFrequenciesSet& mmfs) :
      AbstractFrequenciesSet(mmfs),
      freqSet_(mmfs.freqSet_->clone()),
      rateFreqs_(mmfs.rateFreqs_)
    {}

    MarkovModulatedFrequenciesSet& operator=(const MarkovModulatedFrequenciesSet& mmfs)
    {
      AbstractFrequenciesSet::operator=(mmfs);
      freqSet_ = mmfs.freqSet_->clone();
      rateFreqs_ = mmfs.rateFreqs_;
      return *this;
    }

    MarkovModulatedFrequenciesSet* clone() const { return new MarkovModulatedFrequenciesSet(*this); }

    virtual ~MarkovModulatedFrequenciesSet() { delete freqSet_; }

  public:
    void setFrequencies(const std::vector<double>& frequencies)
    {
      // Just forward this method to the sequence state frequencies set. This may change in the future...
      freqSet_->setFrequencies(frequencies);
    }

    void fireParameterChanged(const ParameterList& pl)
    {
      freqSet_->matchParametersValues(pl);
      setFrequencies_(VectorTools::kroneckerMult(rateFreqs_, freqSet_->getFrequencies()));
    }

    const FrequenciesSet& getStatesFrequenciesSet() const { return *freqSet_; }

  };


/**
 * @brief FrequenciesSet useful for homogeneous and stationary models.
 *
 * This set contains no parameter.
 */
  class FixedFrequenciesSet :
    public AbstractFrequenciesSet
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
    FixedFrequenciesSet(std::shared_ptr<const StateMap> stateMap, const std::vector<double>& initFreqs, const std::string& name = "Fixed");

    /**
     * @brief Construction with uniform frequencies on the states of the model.
     *
     * @param stateMap The model states for which frequencies should be built.
     * @param name The name of the set.
     */
    FixedFrequenciesSet(std::shared_ptr<const StateMap> stateMap, const std::string& name = "Fixed");

    FixedFrequenciesSet* clone() const { return new FixedFrequenciesSet(*this); }

  public:
    void setFrequencies(const std::vector<double>& frequencies);

  };


  /**
   * @brief FrequenciesSet to be read in a file. More specifically, a
   * frequency set is read in a column of a given file, which column
   * number is given in argument (default: 1).
   *
   */
  
  class UserFrequenciesSet :
    public AbstractFrequenciesSet
  {
  private:
    std::string path_;
    size_t nCol_;

  public:
    UserFrequenciesSet(std::shared_ptr<const StateMap> stateMap, const std::string& path, size_t nCol=1);

    UserFrequenciesSet(const UserFrequenciesSet& fmfs);
    
    UserFrequenciesSet& operator=(const UserFrequenciesSet& fmfs);

    UserFrequenciesSet* clone() const { return new UserFrequenciesSet(*this); }

    ~UserFrequenciesSet(){};

  public:

    const std::string& getPath() const { return path_; }

    size_t getColumnNumber() const 
    {
      return nCol_;
    }
    
    void setFrequencies(const std::vector<double>& frequencies);

    // void fireParameterChanged(const ParameterList& pl){};

  protected:
    void readFromFile_();
  };


} // end of namespace bpp.

#endif // _FREQUENCIESSET_H_


