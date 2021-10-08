//
// File: CodonFrequencySet.h
// Created by: laurent Gueguen
// Created on: lundi 2 avril 2012, à 14h 03
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

#ifndef _CODONFREQUENCYSET_H_
#define _CODONFREQUENCYSET_H_

#include "WordFrequencySet.h"
#include "FrequencySet.h"
#include "ProteinFrequencySet.h"

#include <Bpp/Seq/Alphabet/CodonAlphabet.h>
#include <Bpp/Seq/GeneticCode/GeneticCode.h>
#include <Bpp/Numeric/Prob/Simplex.h>

namespace bpp
{
/**
 * @brief Parametrize a set of state frequencies for codons.
 */
class CodonFrequencySet :
  public virtual FrequencySet
{
public:
  CodonFrequencySet* clone() const = 0;

  virtual const CodonAlphabet* getCodonAlphabet() const = 0;

public:
  /**
   * @return The associated genetic code.
   */
  virtual const GeneticCode* getGeneticCode() const = 0;

  /**
   * @brief A helper function that provide frequencies set for codon models
   * according to PAML option.
   *
   * @param option A code describing the option, one of F61, F1X4 or F3X4.
   * @param gCode The genetic code to use. The underlying codon alphabet object will be passed to the FrequencySet instance.
   * @param mgmtStopCodon the optional way the frequencies assigned
   * to the stop codons are redistributed to the other codons, with
   * F1X4 and F3X4 options. The available values are:
   *  - uniform : each stop frequency is distributed evenly
   *  - linear : each stop frequency is distributed to the neighbour
   *     codons (ie 1 substitution away), in proportion to each
   *     target codon frequency.
   *  - quadratic (default): each stop frequency is distributed to the
   *     neighbour codons (ie 1 substitution away), in proportion to
   *     the square of each target codon frequency.
   * @param method The parametrization used for F61. Default method
   * is 1 (ie global ratio).
   *
   * @see Simplex
   */
  static std::shared_ptr<FrequencySet> getFrequencySetForCodons(short option, const GeneticCode* gCode, const std::string& mgmtStopCodon = "quadratic", unsigned short method = 1);

  static const short F0;
  static const short F1X4;
  static const short F3X4;
  static const short F61;
};


/**
 * @brief A generic FrequencySet for Full Codon alphabets.
 *
 * It is very similar to FullFrequencySet, but only the non-stop codon
 *   frequencies are parameterized.
 */
class FullCodonFrequencySet :
  public virtual CodonFrequencySet,
  public AbstractFrequencySet
{
protected:
  const GeneticCode* pgc_;

private:
  /**
   * @brief Simplex to handle the probabilities and the parameters.
   *
   */

  Simplex sFreq_;

public:
  /**
   * @brief Construction with uniform frequencies on the letters of
   * the alphabet. The stop codon frequencies are null.
   */
  FullCodonFrequencySet(const GeneticCode* gCode, bool allowNullFreqs = false, unsigned short method = 1, const std::string& name = "Full");
  FullCodonFrequencySet(const GeneticCode* gCode, const std::vector<double>& initFreqs, bool allowNullFreqs = false, unsigned short method = 1, const std::string& name = "Full");

  FullCodonFrequencySet(const FullCodonFrequencySet& fcfs);
  FullCodonFrequencySet& operator=(const FullCodonFrequencySet& fcfs);

  FullCodonFrequencySet* clone() const { return new FullCodonFrequencySet(*this); }

public:
  const GeneticCode* getGeneticCode() const { return pgc_; }

  /**
   * @brief the given frequencies are normalized such that the sum of
   * the frequencies on the non-stop codons equals 1.
   *
   */
  void setFrequencies(const std::vector<double>& frequencies);

  const CodonAlphabet* getCodonAlphabet() const
  {
    return dynamic_cast<const CodonAlphabet*>(AbstractFrequencySet::getAlphabet());
  }

  void setNamespace(const std::string& nameSpace);

  unsigned short getMethod() const
  {
    return sFreq_.getMethod();
  }

protected:
  void fireParameterChanged(const ParameterList& parameters);

  void updateFreq_();
};


/**
 * @brief FrequencySet useful for homogeneous and stationary models, codon implementation
 *
 * This set contains no parameter.
 */
class FixedCodonFrequencySet :
  public virtual CodonFrequencySet,
  public AbstractFrequencySet
{
protected:
  const GeneticCode* pgc_;

public:
  FixedCodonFrequencySet(const GeneticCode* gCode, const std::vector<double>& initFreqs, const std::string& name = "Fixed");

  /**
   * @brief Construction with uniform frequencies on the letters of
   * the alphabet. The stop codon frequencies are null.
   */
  FixedCodonFrequencySet(const GeneticCode* gCode, const std::string& name = "Fixed");

  FixedCodonFrequencySet(const FixedCodonFrequencySet& fcfs) :
    AbstractFrequencySet(fcfs),
    pgc_(fcfs.pgc_)
  {}

  FixedCodonFrequencySet& operator=(const FixedCodonFrequencySet& fcfs)
  {
    AbstractFrequencySet::operator=(fcfs);
    pgc_ = fcfs.pgc_;
    return *this;
  }

  FixedCodonFrequencySet* clone() const { return new FixedCodonFrequencySet(*this); }

public:
  const GeneticCode* getGeneticCode() const { return pgc_; }

  const CodonAlphabet* getCodonAlphabet() const
  {
    return dynamic_cast<const CodonAlphabet*>(AbstractFrequencySet::getAlphabet());
  }

  /**
   * @brief the given frequencies are normalized such thaat the sum of
   * the frequencies on the non-stop codons equals 1.
   *
   */
  void setFrequencies(const std::vector<double>& frequencies);

protected:
  void fireParameterChanged(const ParameterList& parameters) {}
};

class UserCodonFrequencySet :
  public virtual CodonFrequencySet,
  public UserFrequencySet
{
protected:
  const GeneticCode* pgc_;

public:
  UserCodonFrequencySet(const GeneticCode* gCode, const std::string& path, size_t nCol = 1);

  UserCodonFrequencySet(const UserCodonFrequencySet& fcfs) :
    UserFrequencySet(fcfs),
    pgc_(fcfs.pgc_)
  {}

  UserCodonFrequencySet& operator=(const UserCodonFrequencySet& fcfs)
  {
    UserFrequencySet::operator=(fcfs);
    pgc_ = fcfs.pgc_;
    return *this;
  }

  UserCodonFrequencySet* clone() const { return new UserCodonFrequencySet(*this); }

public:
  const GeneticCode* getGeneticCode() const { return pgc_; }

  const CodonAlphabet* getCodonAlphabet() const
  {
    return dynamic_cast<const CodonAlphabet*>(AbstractFrequencySet::getAlphabet());
  }

  /**
   * @brief the given frequencies are normalized such thaat the sum of
   * the frequencies on the non-stop codons equals 1.
   *
   */
  void setFrequencies(const std::vector<double>& frequencies);

protected:
  void fireParameterChanged(const ParameterList& parameters) {}
};

/**
 * @brief FrequencySet integrating ProteinFrequencySet inside
 * CodonFrequencySet. In this case, FrequencieSet defined inside
 * each amino acid is parametrized as a FullFrequencySet. Hence
 * there are 61-20=41 parameters in addition of the parameters of the
 * ProteinFrequencySet.
 *
 * The parametrization depends on the method used.
 * Default method is 1 (ie global ratio).
 *
 * @see Simplex
 *
 */
class FullPerAACodonFrequencySet :
  public virtual CodonFrequencySet,
  public AbstractFrequencySet
{
private:
  const GeneticCode* pgc_;
  std::shared_ptr<ProteinFrequencySet> ppfs_;

  /**
   * @ brief vector of the simplexes, one for each AA
   */
  std::vector<Simplex> vS_;

  void updateFrequencies();

public:
  /**
   * @brief Create a new FullPerAACodonFrequencySet object.
   *
   * @param gencode The genetic code to use.
   * @param ppfs The protein frequencies to use. The codon
   * frequencies set will own the instance of the protein
   * frequencies set.
   * @param method the method used for parametrization.
   */
  FullPerAACodonFrequencySet(const GeneticCode* gencode, std::shared_ptr<ProteinFrequencySet> ppfs, unsigned short method = 1);

  /**
   * @brief Construction with fixed uniform frequencies on the amino acids.
   * The stop codon frequencies are null.
   * @param gencode The genetic code to use.
   * @param method the method used for parametrization.
   */

  FullPerAACodonFrequencySet(const GeneticCode* gencode, unsigned short method = 1);

  FullPerAACodonFrequencySet(const FullPerAACodonFrequencySet& ffs);

  FullPerAACodonFrequencySet& operator=(const FullPerAACodonFrequencySet& ffs);

  virtual ~FullPerAACodonFrequencySet() {}

  FullPerAACodonFrequencySet* clone() const { return new FullPerAACodonFrequencySet(*this); }

public:
  const CodonAlphabet* getCodonAlphabet() const
  {
    return dynamic_cast<const CodonAlphabet*>(AbstractFrequencySet::getAlphabet());
  }

  const GeneticCode* getGeneticCode() const { return pgc_; }

  /**
   * @brief the given frequencies are normalized such thaat the sum of
   * the frequencies on the non-stop codons equals 1.
   *
   */
  void setFrequencies(const std::vector<double>& frequencies);

  void setNamespace(const std::string& prefix);

  const std::shared_ptr<ProteinFrequencySet> getProteinFrequencySet() const
  {
    return ppfs_;
  }

  unsigned short getMethod() const
  {
    return vS_.size() > 0 ? vS_[0].getMethod() : static_cast<unsigned short>(1);
  }

protected:
  void fireParameterChanged(const ParameterList& parameters);
};


/**
 * @brief the Frequencies in codons are the product of Independent
 * Frequencies in letters with the frequencies of stop codons set to
 * zero.
 *
 * @author Laurent Guéguen
 */
class CodonFromIndependentFrequencySet :
  public virtual CodonFrequencySet,
  public WordFromIndependentFrequencySet
{
private:
  // a map associating stop codons numbers with numbers of neighbour non-stop codons
  std::map<int, Vint> mStopNeigh_;

  unsigned short mgmtStopCodon_;

  const GeneticCode* pgc_;

public:
  /**
   * @brief Constructor from a CodonAlphabet* and a vector of different std::shared_ptr<FrequencySet>.
   * Throws an Exception if their lengths do not match.
   *
   * @param gCode a pointer to the genetic code to use.
   * @param freqvector a vector of pointers to the phase specific FrequencySets
   * @param name the optional name of the FrequencySet (default codon)
   * @param mgmtStopCodon the optional way the frequencies assigned to the
   * stop codons are redistributed to the other codons. The
   * available values are:
   *  - uniform : each stop frequency is distributed evenly
   *  - linear : each stop frequency is distributed to the neighbour
   *     codons (ie 1 substitution away), in proportion to each
   *     target codon frequency.
   *  - quadratic (default): each stop frequency is distributed to the
   *     neighbour codons (ie 1 substitution away), in proportion to
   *     the square of each target codon frequency.
   *
   */
  CodonFromIndependentFrequencySet(const GeneticCode* gCode, const std::vector<std::shared_ptr<FrequencySet> >& freqvector, const std::string& name = "Codon", const std::string& mgmtStopCodon = "quadratic");

  CodonFromIndependentFrequencySet(const CodonFromIndependentFrequencySet& iwfs);

  virtual ~CodonFromIndependentFrequencySet(){}

  CodonFromIndependentFrequencySet& operator=(const CodonFromIndependentFrequencySet& iwfs);

  CodonFromIndependentFrequencySet* clone() const { return new CodonFromIndependentFrequencySet(*this); }

  const CodonAlphabet* getCodonAlphabet() const;

  const GeneticCode* getGeneticCode() const { return pgc_; }

  /**
   * @brief Update the frequencies given the parameters.
   */
  void updateFrequencies();

  /**
   * @brief Retrieve the mgmt method for the stop codons.
   */
  std::string getMgmtStopCodon() const
  {
    switch (mgmtStopCodon_)
    {
    case 0:
      return "uniform";
    case 1:
      return "linear";
    case 2:
      return "quadratic";
    }
    return "";
  }
};


/**
 * @brief the Frequencies in codons are the product of the frequencies
 * for a unique FrequencySet in letters, with the frequencies of
 * stop codons set to zero.
 *
 * @author Laurent Guéguen
 */

class CodonFromUniqueFrequencySet :
  public virtual CodonFrequencySet,
  public WordFromUniqueFrequencySet
{
private:
  // a map associating stop codons numbers with numbers of neighbour non-stop codons
  std::map<int, Vint> mStopNeigh_;

  unsigned short mgmtStopCodon_;

  const GeneticCode* pgc_;

public:
  /**
   * @brief Constructor from a CodonAlphabet* and a std::shared_ptr<FrequencySet>
   *  repeated three times.
   *
   * @param gCode a pointer to a genetic code.
   * @param pfreq a pointer to the nucleotidic FrequencySet
   * @param name the optional name of the FrequencySet (default codon)
   * @param mgmtStopCodon the optional way the frequencies assigned to the
   * stop codons are redistributed to the other codons. The
   * available values are:
   *  - uniform : each stop frequency is distributed evenly
   *  - linear : each stop frequency is distributed to the neighbour
   *      codons (ie 1 substitution away), in proportion to each
   *      target codon frequency.
   *  - quadratic (default): each stop frequency is distributed to the
   *      neighbour codons (ie 1 substitution away), in proportion to
   *      the square of each target codon frequency.
   */
  CodonFromUniqueFrequencySet(
    const GeneticCode* gCode,
    std::shared_ptr<FrequencySet> pfreq,
    const std::string& name = "Codon",
    const std::string& mgmtStopCodon = "quadratic");

  CodonFromUniqueFrequencySet(const CodonFromUniqueFrequencySet& iwfs);

  virtual ~CodonFromUniqueFrequencySet() {}

  CodonFromUniqueFrequencySet& operator=(const CodonFromUniqueFrequencySet& iwfs);

  CodonFromUniqueFrequencySet* clone() const { return new CodonFromUniqueFrequencySet(*this); }

  const CodonAlphabet* getCodonAlphabet() const;

  const GeneticCode* getGeneticCode() const { return pgc_; }

  /**
   * @brief Update the frequencies given the parameters.
   *
   */
  void updateFrequencies();

  /**
   * @brief Retrieve the mgmt method for the stop codons.
   */
  std::string getMgmtStopCodon() const
  {
    switch (mgmtStopCodon_)
    {
    case 0:
      return "uniform";
    case 1:
      return "linear";
    case 2:
      return "quadratic";
    }
    return "";
  }
};
} // end of namespace bpp.

#endif// _CODONFREQUENCYSET_H_
