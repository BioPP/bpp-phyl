// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_NUCLEOTIDE_GBGC_H
#define BPP_PHYL_MODEL_NUCLEOTIDE_GBGC_H

#include <Bpp/Numeric/Constraints.h>

#include "../AbstractSubstitutionModel.h"
#include "NucleotideSubstitutionModel.h"

// From SeqLib:
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>

// From Utils:

using namespace std;

namespace bpp
{
/**
 * @brief gBGC model.
 *
 * Modelling of GC biased gene-conversion.
 *
 * This model adds strand symmetric GC biased gene conversion to a
 * given nucleotidic substitution model.
 *
 * In addition to the parameters of the basic nucleic model, the
 * biased gene conversion effect is parametrized by @f$ B @f$, that
 * corresponds to fixation, and stands for a dubious selection.
 *
 * With this term, the mutation rates from A and T to C and G are
 * multiplied by @f$ B . \frac{1}{1-exp(-B)}@f$, and
 * the mutation rates from C and G to A and T are multiplied by
 * @f$ B . \frac{1}{exp(B)-1}@f$.
 *
 * @see AbstractSubstitutionModel
 *
 * Reference:
 * - Galtier & al (2009), Trends in Genetics, 25(1), doi:10.1016/j.tig.2008.10.011
 */

class gBGC :
  public AbstractNucleotideSubstitutionModel
{
private:
  std::unique_ptr<NucleotideSubstitutionModelInterface> model_;
  std::string nestedPrefix_;

  /**
   * @brief the value of the bias.
   */
  double B_;

public:
  /**
   * @brief Build a new gBGC substitution model.
   */
  gBGC(
      std::shared_ptr<const NucleicAlphabet>,
      std::unique_ptr<NucleotideSubstitutionModelInterface>,
      double B = 0);

  gBGC(const gBGC&);

  gBGC& operator=(const gBGC& gbgc);

  gBGC* clone() const override { return new gBGC(*this); }

  virtual ~gBGC() {}

public:
  std::string getName() const override
  {
    return model_->getName() + "+gBGC";
  }

  size_t getNumberOfStates() const override { return model_->getNumberOfStates(); }

  void fireParameterChanged(const ParameterList&) override;

  const SubstitutionModelInterface& nestedModel() const { return *model_; }

  void setNamespace(const std::string&) override;

protected:
  void updateMatrices_() override;
};
}
#endif // BPP_PHYL_MODEL_NUCLEOTIDE_GBGC_H
