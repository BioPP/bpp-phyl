//
// File: WAG01.h
// Authors:
//   Laurent Gueguen
// Created: mardi 28 septembre 2010, Ã  14h 43
//

/*
  Copyright or Â© or Copr. Bio++ Development Team, (November 16, 2004)
  
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

#ifndef BPP_PHYL_MODEL_PROTEIN_WAG01_H
#define BPP_PHYL_MODEL_PROTEIN_WAG01_H


#include "../AbstractSubstitutionModel.h"
#include "../FrequencySet/ProteinFrequencySet.h"
#include "ProteinSubstitutionModel.h"

// From SeqLib:
#include <Bpp/Seq/Alphabet/ProteicAlphabet.h>

namespace bpp
{
/**
 * @brief The Whelan and Goldman substitution model for proteins.
 *
 * The exchangability matrix is normalized so that \f$Q = S . \pi\f$ and
 * \f$\sum_i Q_{i,i}\pi_i = -1\f$.
 *
 * The original frequencies can be used, or alternatively a
 * parametrized version, corresponding to the so-called WAG01+F
 * model. Eigen values and vectors are obtained numerically.
 *
 * Reference:
 *
 * Whelan, S. and N. Goldman. 2001. A general empirical model of
 * protein evolution derived from multiple protein families using a
 * maximum likelihood approach. Molecular Biology and Evolution 18:691-699.
 *
 */
class WAG01 :
  public AbstractReversibleProteinSubstitutionModel
{
private:
  std::shared_ptr<ProteinFrequencySet> freqSet_;

public:
  /**
   * @brief Build a simple WAG01 model, with original equilibrium frequencies.
   *
   * @param alpha A proteic alphabet.
   */
  WAG01(const ProteicAlphabet* alpha);

  /**
   * @brief Build a WAG01 model with special equilibrium frequencies.
   *
   * @param alpha A proteic alphabet.
   * @param freqSet A pointer toward a protein frequencies set, which will be owned by this instance.
   * @param initFreqs Tell if the frequency set should be initialized with the original WAG01 values.
   * Otherwise, the values of the set will be used.
   */
  WAG01(const ProteicAlphabet* alpha, std::shared_ptr<ProteinFrequencySet> freqSet, bool initFreqs = false);

  WAG01(const WAG01& model) :
    AbstractParameterAliasable(model),
    AbstractReversibleProteinSubstitutionModel(model),
    freqSet_(dynamic_cast<ProteinFrequencySet*>(model.freqSet_->clone()))
  {}

  WAG01& operator=(const WAG01& model)
  {
    AbstractParameterAliasable::operator=(model);
    AbstractReversibleProteinSubstitutionModel::operator=(model);
    freqSet_.reset(dynamic_cast<ProteinFrequencySet*>(model.freqSet_->clone()));
    return *this;
  }

  virtual ~WAG01() {}

  WAG01* clone() const { return new WAG01(*this); }

public:
  std::string getName() const
  {
    if (freqSet_->getNamespace().find("WAG01+F.") != std::string::npos)
      return "WAG01+F";
    else
      return "WAG01";
  }

  void fireParameterChanged(const ParameterList& parameters)
  {
    freqSet_->matchParametersValues(parameters);
    freq_ = freqSet_->getFrequencies();
    AbstractReversibleSubstitutionModel::fireParameterChanged(parameters);
  }

  void setNamespace(const std::string& prefix)
  {
    AbstractParameterAliasable::setNamespace(prefix);
    freqSet_->setNamespace(prefix + freqSet_->getName() + ".");
  }

  void setFrequencySet(const ProteinFrequencySet& freqSet)
  {
    freqSet_ = std::shared_ptr<ProteinFrequencySet>(dynamic_cast<ProteinFrequencySet*>(freqSet.clone()));
    resetParameters_();
    addParameters_(freqSet_->getParameters());
  }

  const std::shared_ptr<FrequencySet> getFrequencySet() const { return freqSet_; }

  void setFreqFromData(const SequencedValuesContainer& data, double pseudoCount = 0);
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_PROTEIN_WAG01_H
