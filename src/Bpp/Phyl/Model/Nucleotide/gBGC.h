//
// File: gBGC.h
// Authors:
//   Laurent Gueguen
// Created: lundi 13 fÃ©vrier 2012, Ã  09h 43
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
 * @author Laurent Guéguen
 *
 * modelling of GC biased gene-conversion.
 *
 * This model adds strand symetric GC biased gene conversion to a
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
  std::unique_ptr<NucleotideSubstitutionModelInterface>  model_;
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
  std::string getName() const override;

  size_t getNumberOfStates() const override { return model_->getNumberOfStates(); }

  void fireParameterChanged(const ParameterList&) override;

  const SubstitutionModelInterface& nestedModel() const {return *model_; }

  void updateMatrices() override;

  void setNamespace(const std::string&) override;
};
}
#endif // BPP_PHYL_MODEL_NUCLEOTIDE_GBGC_H
