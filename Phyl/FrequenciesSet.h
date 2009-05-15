//
// File: FrequenciesSet.h
// Created by: Bastien Boussau
//             Julien Dutheil
// Created on: Tue Aug 21 2007
//

/*
Copyright or (c) or Copr. CNRS, (November 16, 2004)

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

// From NumCalc:
#include <NumCalc/Parametrizable.h>
#include <NumCalc/VectorTools.h>
#include <NumCalc/AbstractParametrizable.h>

// From SeqLib:
#include <Seq/Alphabet.h>
#include <Seq/NucleicAlphabet.h>
#include <Seq/ProteicAlphabet.h>

using namespace std;

namespace bpp
{

/**
 * @brief Parametrize a set of state frequencies.
 */
class FrequenciesSet:
  public virtual Parametrizable
{
  public:
#ifndef NO_VIRTUAL_COV
    FrequenciesSet * clone() const = 0;
#endif
  public:
    /**
     * @return The alphabet associated to this set.
     */  
    virtual const Alphabet * getAlphabet() const = 0;

    /**
     * @return The frequencies values of the set.
     */ 
    virtual const vector<double>& getFrequencies() const = 0;

    /**
     * @brief Set the parameters in order to match a given set of frequencies.
     *
     * @param frequencies The set of frequencies to match.
     * @throw DimensionException If the number of frequencies does not match the size of the alphabet.
     * @throw Exception If the frequencies do not sum to 1.
     */
    virtual void setFrequencies(const vector<double> & frequencies) throw (DimensionException, Exception) = 0;

};

/**
 * @brief Parametrize a set of state frequencies for nucleotides.
 */
class NucleotideFrequenciesSet:
  public virtual FrequenciesSet
{
  public:
#ifndef NO_VIRTUAL_COV
    NucleotideFrequenciesSet * clone() const = 0;
#endif

};

/**
 * @brief Parametrize a set of state frequencies for proteins.
 */
class ProteinFrequenciesSet:
  public virtual FrequenciesSet
{
  public:
#ifndef NO_VIRTUAL_COV
    ProteinFrequenciesSet * clone() const = 0;
#endif

};



/**
 * @brief Basic implementation of the FrequenciesSet interface.
 */
class AbstractFrequenciesSet:
  public virtual FrequenciesSet,
  public AbstractParametrizable
{
  private:
    const Alphabet * alphabet_;
    vector<double> freq_;

  public:
    AbstractFrequenciesSet(unsigned int n, const Alphabet * alphabet, const string& prefix) :
      AbstractParametrizable(prefix), alphabet_(alphabet), freq_(n) {}
  
  public:
    const Alphabet * getAlphabet() const { return alphabet_; }
    const vector<double>& getFrequencies() const { return freq_; }

  protected:
    vector<double>& getFrequencies_() { return freq_; }
    double& getFreq_(unsigned int i) { return freq_[i]; }
    const double& getFreq_(unsigned int i) const { return freq_[i]; }
    void setFrequencies_(const vector<double>& frequencies) { freq_ = frequencies; }
};

/**
 * @brief FrequenciesSet using one parameter per frequency.
 */
class FullFrequenciesSet:
  public AbstractFrequenciesSet
{
  public:
    FullFrequenciesSet(const Alphabet * alphabet, const string& name, const string& prefix = "");
    FullFrequenciesSet(const Alphabet * alphabet, const vector<double> & initFreqs, const string& name, const string& prefix = "") throw (Exception);

#ifndef NO_VIRTUAL_COV
    FullFrequenciesSet *
#else
    Clonable *
#endif
    clone() const { return new FullFrequenciesSet(*this); }

  public:
    void fireParameterChanged(const ParameterList & pl)
    {
      for(unsigned int i = 0; i < getAlphabet()->getSize(); i++)
      {
        getFreq_(i) = getParameter_(i).getValue();
      }
    }

    void setFrequencies(const vector<double> & frequencies) throw (DimensionException, Exception)
    {
      if(frequencies.size() != getAlphabet()->getSize()) throw DimensionException("FullFrequenciesSet::setFrequencies", frequencies.size(), getAlphabet()->getSize());
      double sum = 0.0;
      for(unsigned int i = 0; i < frequencies.size(); i++)
        sum += frequencies[i];
      if(fabs(1.-sum) > 0.000001)
        throw Exception("FullFrequenciesSet::setFrequencies. Frequencies must equal 1 (sum = " + TextTools::toString(sum) + ").");
      for(unsigned int i = 0; i < getNumberOfParameters(); i++)
      {
        getParameter_(i).setValue(frequencies[i]);
        getFreq_(i) = frequencies[i];
      }
    }
};

/**
 * @brief Nucleotide FrequenciesSet using only one parameter, the GC content.
 */
class GCFrequenciesSet:
  public NucleotideFrequenciesSet, public AbstractFrequenciesSet
{
  public:
    GCFrequenciesSet(const NucleicAlphabet * alphabet, const string & prefix = ""):
      AbstractFrequenciesSet(4, alphabet, prefix)
    {
      Parameter p(getNamespace() + "theta", 0.5, &Parameter::PROP_CONSTRAINT_IN);
      addParameter_(p);
      getFreq_(0) = getFreq_(1) = getFreq_(2) = getFreq_(3) = 0.25;
    }
    GCFrequenciesSet(const NucleicAlphabet * alphabet, double theta, const string & prefix = ""):
      AbstractFrequenciesSet(4, alphabet, prefix)
    {
      Parameter p(getNamespace() + "theta", theta, &Parameter::PROP_CONSTRAINT_IN);
      addParameter_(p);
      getFreq_(0) = getFreq_(3) = (1. - theta) / 2.;
      getFreq_(1) = getFreq_(2) = theta / 2.;
    }
#ifndef NO_VIRTUAL_COV
    GCFrequenciesSet *
#else
    Clonable *
#endif
    clone() const { return new GCFrequenciesSet(*this); }

  public:
    void fireParameterChanged(const ParameterList & pl)
    {
      double theta = getParameter_(0).getValue();
      getFreq_(0) = getFreq_(3) = (1. - theta) / 2.;
      getFreq_(1) = getFreq_(2) = theta / 2.;
    }

    void setFrequencies(const vector<double> & frequencies) throw (DimensionException)
    {
      if(frequencies.size() != 4) throw DimensionException("GCFrequenciesSet::setFrequencies", frequencies.size(), 4);
      double sum = 0.0;
      for(unsigned int i = 0; i < 4; i++)
        sum += frequencies[i];
      if(fabs(1.-sum) > 0.000001)
        throw Exception("GCFrequenciesSet::setFrequencies. Frequencies must equal 1 (sum = " + TextTools::toString(sum) + ").");
      getParameter_(0).setValue(frequencies[1] + frequencies[2]);
      setFrequencies_(frequencies);
    }
};

/**
 * @brief Nucleotide FrequenciesSet using three indpeendent parameters to modelize the four frequencies.
 */
class FullNAFrequenciesSet:
  public NucleotideFrequenciesSet, public AbstractFrequenciesSet
{
  public:
    FullNAFrequenciesSet(const NucleicAlphabet * alphabet, const string & prefix = "");

    FullNAFrequenciesSet(const NucleicAlphabet * alphabet, double theta, double theta1, double theta2, const string & prefix = "");

#ifndef NO_VIRTUAL_COV
    FullNAFrequenciesSet *
#else
    Clonable *
#endif
    clone() const { return new FullNAFrequenciesSet(*this); }

  public:
    void setFrequencies(const vector<double> & frequencies) throw (DimensionException, Exception)
    {
      if(frequencies.size() != 4) throw DimensionException("FullNAFrequenciesSet::setFrequencies", frequencies.size(), 4);
      double sum = 0.0;
      for(unsigned int i = 0; i < 4; i++)
        sum += frequencies[i];
      if(fabs(1.-sum) > 0.000001)
        throw Exception("FullNAFrequenciesSet::setFrequencies. Frequencies must equal 1 (sum = " + TextTools::toString(sum) + ").");
      double theta = frequencies[1] + frequencies[2];
      getParameter_(0).setValue(theta);
      getParameter_(1).setValue(frequencies[0] / (1 - theta));
      getParameter_(2).setValue(frequencies[2] / theta);
      setFrequencies_(frequencies);
    }

    void fireParameterChanged(const ParameterList & pl);
};

/**
 * @brief Protein FrequenciesSet using 19 independent parameters to modelize the 20 frequencies.
 *
 * The parameters are called @f$ \theta_{i \in 1..19} @f$, and are initialized so that all frequencies are equal to  0.005, that is
 * @f[ \theta_i = \frac{0.05}{0.956{i-1}},\quad i = 1..19 @f] or according to a user-specified vector of initial values.
 */
class FullProteinFrequenciesSet:
  public ProteinFrequenciesSet, public AbstractFrequenciesSet
{
  public:
    FullProteinFrequenciesSet(const ProteicAlphabet * alphabet, const string & prefix = "");
    FullProteinFrequenciesSet(const ProteicAlphabet * alphabet, const vector<double> & initFreqs, const string & prefix = "") throw (Exception);

#ifndef NO_VIRTUAL_COV
    FullProteinFrequenciesSet *
#else
    Clonable *
#endif
    clone() const { return new FullProteinFrequenciesSet(*this); }

  public:
    void setFrequencies(const vector<double> & frequencies) throw (DimensionException, Exception);

    void fireParameterChanged(const ParameterList & pl);
};

/**
 * @brief FrequenciesSet to be used with a Markov-modulated substitution model.
 * 
 * This implementation uses one parameter per character state frequency.
 * The rate states are assumed to be fixed and are passed as an argument to the constructor, together with a 'regular'
 * FrequenciesSet. The number of parameters hence do not depends on the number of rates used.
 */
class MarkovModulatedFrequenciesSet:
  public AbstractFrequenciesSet
{
  protected:
    FrequenciesSet * _freqSet;
    vector<double> _rateFreqs;
  public:
    MarkovModulatedFrequenciesSet(const MarkovModulatedFrequenciesSet & mmfs):
      AbstractFrequenciesSet(mmfs)
    {
      _freqSet = dynamic_cast<FrequenciesSet *>(mmfs._freqSet->clone());
      _rateFreqs = mmfs._rateFreqs;
    }
    MarkovModulatedFrequenciesSet & operator=(const MarkovModulatedFrequenciesSet & mmfs)
    {
      AbstractFrequenciesSet::operator=(mmfs);
      _freqSet = dynamic_cast<FrequenciesSet *>(mmfs._freqSet->clone());
      _rateFreqs = mmfs._rateFreqs;
      return *this;
    }
    MarkovModulatedFrequenciesSet(FrequenciesSet * freqSet, const vector<double> & rateFreqs):
      AbstractFrequenciesSet(getAlphabet()->getSize() * rateFreqs.size(), freqSet->getAlphabet(), ""), _freqSet(freqSet), _rateFreqs(rateFreqs)
    {
      addParameters_(_freqSet->getParameters());
      setFrequencies_(VectorTools::kroneckerMult(rateFreqs, _freqSet->getFrequencies()));
    }
#ifndef NO_VIRTUAL_COV
    MarkovModulatedFrequenciesSet *
#else
    Clonable *
#endif
    clone() const { return new MarkovModulatedFrequenciesSet(*this); }

    virtual ~MarkovModulatedFrequenciesSet() { delete _freqSet; }

  public:
    void setFrequencies(const vector<double> & frequencies) throw (DimensionException)
    {
      //Just forward this method to the sequence state frequencies set. This may change in the future...
      _freqSet->setFrequencies(frequencies);
    }

    void fireParameterChanged(const ParameterList & pl)
    {
      _freqSet->matchParametersValues(pl);
      setFrequencies_(VectorTools::kroneckerMult(_rateFreqs, _freqSet->getFrequencies()));
    }

    const FrequenciesSet* getStatesFrequenciesSet() const { return _freqSet; }
};

/**
 * @brief FrequenciesSet useful for homogeneous and stationary models.
 *
 * This set contains no parameter.
 */
class FixedFrequenciesSet:
  public AbstractFrequenciesSet
{
  public:
    FixedFrequenciesSet(const Alphabet * alphabet, const vector<double>& initFreqs, const string & prefix = "");

#ifndef NO_VIRTUAL_COV
    FixedFrequenciesSet *
#else
    Clonable *
#endif
    clone() const { return new FixedFrequenciesSet(*this); }

  public:
    void setFrequencies(const vector<double> & frequencies) throw (DimensionException)
    {
      if(frequencies.size() != getAlphabet()->getSize()) throw DimensionException("FixedFrequenciesSet::setFrequencies", frequencies.size(), getAlphabet()->getSize());
      double sum = 0.0;
      for(unsigned int i = 0; i < frequencies.size(); i++)
        sum += frequencies[i];
      if(fabs(1.-sum) > 0.000001)
        throw Exception("FixedFrequenciesSet::setFrequencies. Frequencies must equal 1 (sum = " + TextTools::toString(sum) + ").");
      setFrequencies_(frequencies);
    }

   void fireParameterChanged(const ParameterList & pl) {}
};

} //end of namespace bpp.

#endif //_FREQUENCIESSET_H_


