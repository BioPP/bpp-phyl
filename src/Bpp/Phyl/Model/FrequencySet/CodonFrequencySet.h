// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_FREQUENCYSET_CODONFREQUENCYSET_H
#define BPP_PHYL_MODEL_FREQUENCYSET_CODONFREQUENCYSET_H

#include <Bpp/Numeric/Prob/Simplex.h>
#include <Bpp/Seq/Alphabet/CodonAlphabet.h>
#include <Bpp/Seq/GeneticCode/GeneticCode.h>

#include "FrequencySet.h"
#include "ProteinFrequencySet.h"
#include "WordFrequencySet.h"

namespace bpp
{
/**
 * @brief Parametrize a set of state frequencies for codons.
 */
class CodonFrequencySetInterface :
  public virtual FrequencySetInterface
{
public:
  CodonFrequencySetInterface* clone() const override = 0;

  virtual std::shared_ptr<const CodonAlphabet> getCodonAlphabet() const = 0;

public:
  /**
   * @return The associated genetic code.
   */
  virtual std::shared_ptr<const GeneticCode> getGeneticCode() const = 0;

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
  static std::unique_ptr<CodonFrequencySetInterface> getFrequencySetForCodons(
    short option,
    std::shared_ptr<const GeneticCode> gCode,
    const std::string& mgmtStopCodon = "quadratic",
    unsigned short method = 1);

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
  public virtual CodonFrequencySetInterface,
  public AbstractFrequencySet
{
protected:
  std::shared_ptr<const GeneticCode> pgc_;

private:
  /**
   * @brief Simplex to handle the probabilities and the parameters.
   */
  Simplex sFreq_;

public:
  /**
   * @brief Construction with uniform frequencies on the letters of
   * the alphabet. The stop codon frequencies are null.
   */
  FullCodonFrequencySet(
      std::shared_ptr<const GeneticCode> gCode,
      bool allowNullFreqs = false,
      unsigned short method = 1,
      const std::string& name = "Full");

  FullCodonFrequencySet(
      std::shared_ptr<const GeneticCode> gCode,
      const std::vector<double>& initFreqs,
      bool allowNullFreqs = false,
      unsigned short method = 1,
      const std::string& name = "Full");

  FullCodonFrequencySet(const FullCodonFrequencySet& fcfs);
  FullCodonFrequencySet& operator=(const FullCodonFrequencySet& fcfs);

  FullCodonFrequencySet* clone() const override { return new FullCodonFrequencySet(*this); }

public:
  std::shared_ptr<const GeneticCode> getGeneticCode() const override { return pgc_; }

  /**
   * @brief the given frequencies are normalized such that the sum of
   * the frequencies on the non-stop codons equals 1.
   */
  void setFrequencies(const std::vector<double>& frequencies) override;

  std::shared_ptr<const CodonAlphabet> getCodonAlphabet() const override
  {
    return std::dynamic_pointer_cast<const CodonAlphabet>(getAlphabet());
  }

  void setNamespace(const std::string& nameSpace) override;

  unsigned short getMethod() const
  {
    return sFreq_.getMethod();
  }

protected:
  void fireParameterChanged(const ParameterList& parameters) override;

  void updateFreq_();
};


/**
 * @brief FrequencySet useful for homogeneous and stationary models, codon implementation
 *
 * This set contains no parameter.
 */
class FixedCodonFrequencySet :
  public virtual CodonFrequencySetInterface,
  public AbstractFrequencySet
{
protected:
  std::shared_ptr<const GeneticCode> pgc_;

public:
  FixedCodonFrequencySet(
      std::shared_ptr<const GeneticCode> gCode,
      const std::vector<double>& initFreqs,
      const std::string& name = "Fixed");

  /**
   * @brief Construction with uniform frequencies on the letters of
   * the alphabet. The stop codon frequencies are null.
   */
  FixedCodonFrequencySet(
      std::shared_ptr<const GeneticCode> gCode,
      const std::string& name = "Fixed");

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

  FixedCodonFrequencySet* clone() const override { return new FixedCodonFrequencySet(*this); }

public:
  std::shared_ptr<const GeneticCode> getGeneticCode() const override { return pgc_; }

  std::shared_ptr<const CodonAlphabet> getCodonAlphabet() const override
  {
    return std::dynamic_pointer_cast<const CodonAlphabet>(getAlphabet());
  }

  /**
   * @brief the given frequencies are normalized such thaat the sum of
   * the frequencies on the non-stop codons equals 1.
   */
  void setFrequencies(const std::vector<double>& frequencies) override;

protected:
  void fireParameterChanged(const ParameterList& parameters) override {}
};

class UserCodonFrequencySet :
  public virtual CodonFrequencySetInterface,
  public UserFrequencySet
{
protected:
  std::shared_ptr<const GeneticCode> pgc_;

public:
  UserCodonFrequencySet(
      std::shared_ptr<const GeneticCode> gCode,
      const std::string& path,
      size_t nCol = 1);

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

  UserCodonFrequencySet* clone() const override { return new UserCodonFrequencySet(*this); }

public:
  std::shared_ptr<const GeneticCode> getGeneticCode() const override { return pgc_; }

  std::shared_ptr<const CodonAlphabet> getCodonAlphabet() const override
  {
    return std::dynamic_pointer_cast<const CodonAlphabet>(getAlphabet());
  }

  /**
   * @brief the given frequencies are normalized such thaat the sum of
   * the frequencies on the non-stop codons equals 1.
   */
  void setFrequencies(const std::vector<double>& frequencies) override;

protected:
  void fireParameterChanged(const ParameterList& parameters) override {}
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
  public virtual CodonFrequencySetInterface,
  public AbstractFrequencySet
{
private:
  std::shared_ptr<const GeneticCode> pgc_;
  std::unique_ptr<ProteinFrequencySetInterface> ppfs_;

  /**
   * @brief vector of the simplexes, one for each AA
   */
  std::vector<Simplex> vS_;

  void updateFrequencies_();

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
  FullPerAACodonFrequencySet(
      std::shared_ptr<const GeneticCode> gencode,
      std::unique_ptr<ProteinFrequencySetInterface> ppfs,
      unsigned short method = 1);

  /**
   * @brief Construction with fixed uniform frequencies on the amino acids.
   * The stop codon frequencies are null.
   * @param gencode The genetic code to use.
   * @param method the method used for parametrization.
   */
  FullPerAACodonFrequencySet(
      std::shared_ptr<const GeneticCode> gencode,
      unsigned short method = 1);

  FullPerAACodonFrequencySet(const FullPerAACodonFrequencySet& ffs);

  FullPerAACodonFrequencySet& operator=(const FullPerAACodonFrequencySet& ffs);

  virtual ~FullPerAACodonFrequencySet() {}

  FullPerAACodonFrequencySet* clone() const override { return new FullPerAACodonFrequencySet(*this); }

public:
  std::shared_ptr<const CodonAlphabet> getCodonAlphabet() const override
  {
    return std::dynamic_pointer_cast<const CodonAlphabet>(getAlphabet());
  }

  std::shared_ptr<const GeneticCode> getGeneticCode() const override { return pgc_; }

  /**
   * @brief the given frequencies are normalized such thaat the sum of
   * the frequencies on the non-stop codons equals 1.
   */
  void setFrequencies(const std::vector<double>& frequencies) override;

  void setNamespace(const std::string& prefix) override;

  bool hasProteinFrequencySet() const
  {
    return ppfs_ != nullptr;
  }

  const ProteinFrequencySetInterface& proteinFrequencySet() const
  {
    return *ppfs_;
  }

  unsigned short getMethod() const
  {
    return vS_.size() > 0 ? vS_[0].getMethod() : static_cast<unsigned short>(1);
  }

protected:
  void fireParameterChanged(const ParameterList& parameters) override;
};


/**
 * @brief the Frequencies in codons are the product of Independent
 * Frequencies in letters with the frequencies of stop codons set to
 * zero.
 *
 * @author Laurent Guéguen
 */
class CodonFromIndependentFrequencySet :
  public virtual CodonFrequencySetInterface,
  public WordFromIndependentFrequencySet
{
private:
  // a map associating stop codons numbers with numbers of neighbour non-stop codons
  std::map<int, Vint> mStopNeigh_;

  unsigned short mgmtStopCodon_;

  std::shared_ptr<const GeneticCode> pgc_;

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
   */
  CodonFromIndependentFrequencySet(
      std::shared_ptr<const GeneticCode> gCode,
      std::vector<std::unique_ptr<FrequencySetInterface>>& freqvector,
      const std::string& name = "Codon",
      const std::string& mgmtStopCodon = "quadratic");

  CodonFromIndependentFrequencySet(const CodonFromIndependentFrequencySet& iwfs);

  virtual ~CodonFromIndependentFrequencySet(){}

  CodonFromIndependentFrequencySet& operator=(const CodonFromIndependentFrequencySet& iwfs);

  CodonFromIndependentFrequencySet* clone() const override { return new CodonFromIndependentFrequencySet(*this); }

  std::shared_ptr<const CodonAlphabet> getCodonAlphabet() const override;

  std::shared_ptr<const GeneticCode> getGeneticCode() const override { return pgc_; }

  /**
   * @brief Update the frequencies given the parameters.
   */
  void updateFrequencies() override;

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
  public virtual CodonFrequencySetInterface,
  public WordFromUniqueFrequencySet
{
private:
  // a map associating stop codons numbers with numbers of neighbour non-stop codons
  std::map<int, Vint> mStopNeigh_;

  unsigned short mgmtStopCodon_;

  std::shared_ptr<const GeneticCode> pgc_;

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
      std::shared_ptr<const GeneticCode> gCode,
      std::unique_ptr<FrequencySetInterface> pfreq,
      const std::string& name = "Codon",
      const std::string& mgmtStopCodon = "quadratic");

  CodonFromUniqueFrequencySet(const CodonFromUniqueFrequencySet& iwfs);

  virtual ~CodonFromUniqueFrequencySet() {}

  CodonFromUniqueFrequencySet& operator=(const CodonFromUniqueFrequencySet& iwfs);

  CodonFromUniqueFrequencySet* clone() const override { return new CodonFromUniqueFrequencySet(*this); }

  std::shared_ptr<const CodonAlphabet> getCodonAlphabet() const override;

  std::shared_ptr<const GeneticCode> getGeneticCode() const override { return pgc_; }

  /**
   * @brief Update the frequencies given the parameters.
   */
  void updateFrequencies() override;

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
#endif // BPP_PHYL_MODEL_FREQUENCYSET_CODONFREQUENCYSET_H
