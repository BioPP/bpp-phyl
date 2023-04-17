//
// File: DSO78.h
// Authors:
//   Julien Dutheil
// Created: 2004-10-05 18:49:44
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

#ifndef BPP_PHYL_MODEL_PROTEIN_DSO78_H
#define BPP_PHYL_MODEL_PROTEIN_DSO78_H


#include "../AbstractSubstitutionModel.h"
#include "../FrequencySet/ProteinFrequencySet.h"
#include "ProteinSubstitutionModel.h"

// From SeqLib:
#include <Bpp/Seq/Alphabet/ProteicAlphabet.h>

namespace bpp
{
/**
 * @brief The Dayhoff, Schwartz and Orcutt substitution model for proteins.
 *
 * Exchangeabilities have been computed using the DCMut method of Kosiol and Goldman.
 * The exchangability matrix is normalized so that \f$Q = S . \pi\f$ and \f$\sum_i Q_{i,i}\pi_i = -1\f$.
 * The original frequencies can be used, or alternatively a parametrized version, corresponding to the
 * so-called JTT92+F model.
 * Eigen values and vectors are obtained numerically.
 *
 * References:
 * - Dayhoff MO, Schwartz RM and Orcutt BC (1978), _A model of evolutionary change in proteins_, 5(3) 345-352, in _Atlas of Protein Sequence and Structure_.
 * - Kosiol C and Goldman N (2005), _Molecular Biology And Evolution_ 22(2) 193-9.
 */
class DSO78 :
  public AbstractReversibleProteinSubstitutionModel
{
private:
  std::unique_ptr<ProteinFrequencySetInterface> freqSet_;

public:
  /**
   * @brief Build a simple DSO78 model, with original equilibrium frequencies.
   *
   * @param alpha A proteic alphabet.
   */
  DSO78(std::shared_ptr<const ProteicAlphabet> alpha);

  /**
   * @brief Build a DSO78 model with special equilibrium frequencies.
   *
   * @param alpha A proteic alphabet.
   * @param freqSet A pointer toward a protein frequencies set, which will be owned by this instance.
   * @param initFreqs Tell if the frequency set should be initialized with the original JTT92 values.
   * Otherwise, the values of the set will be used.
   */
  DSO78(std::shared_ptr<const ProteicAlphabet> alpha,
        std::unique_ptr<ProteinFrequencySetInterface> freqSet,
       	bool initFreqs = false);

  DSO78(const DSO78& model) :
    AbstractParameterAliasable(model),
    AbstractReversibleProteinSubstitutionModel(model),
    freqSet_(model.freqSet_->clone())
  {}

  DSO78& operator=(const DSO78& model)
  {
    AbstractParameterAliasable::operator=(model);
    AbstractReversibleProteinSubstitutionModel::operator=(model);
    freqSet_.reset(model.freqSet_->clone());
    return *this;
  }

  virtual ~DSO78() {}

  DSO78* clone() const override { return new DSO78(*this); }

public:
  std::string getName() const override
  {
    if (freqSet_->getNamespace().find("DSO78+F.") != std::string::npos)
      return "DSO78+F";
    else
      return "DSO78";
  }

  void fireParameterChanged(const ParameterList& parameters) override
  {
    freqSet_->matchParametersValues(parameters);
    freq_ = freqSet_->getFrequencies();
    AbstractReversibleSubstitutionModel::fireParameterChanged(parameters);
  }

  void setNamespace(const std::string& prefix) override
  {
    AbstractParameterAliasable::setNamespace(prefix);
    freqSet_->setNamespace(prefix + freqSet_->getName() + ".");
  }

  void setFrequencySet(const ProteinFrequencySetInterface& freqSet)
  {
    freqSet_.reset(freqSet.clone());
    resetParameters_();
    addParameters_(freqSet_->getParameters());
  }

  const FrequencySetInterface& frequencySet() const override
  {
    if (freqSet_)
      return *freqSet_;
    throw NullPointerException("DSO78::frequencySet(). No associated FrequencySet.");
  }
  
  void setFreqFromData(const SequenceDataInterface& data, double pseudoCount = 0) override;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_PROTEIN_DSO78_H
