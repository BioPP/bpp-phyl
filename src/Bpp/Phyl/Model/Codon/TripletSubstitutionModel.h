// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_TRIPLETSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_CODON_TRIPLETSUBSTITUTIONMODEL_H


#include "../Nucleotide/NucleotideSubstitutionModel.h"
#include "../WordSubstitutionModel.h"

// From SeqLib:
#include <Bpp/Seq/Alphabet/CodonAlphabet.h>

namespace bpp
{
/**
 * @brief Class for neutral substitution models on triplets,
 * which correspond to codons that do not have any significance
 * (whether they are STOP or functional).
 * @author Laurent GuÃÂ©guen
 *
 * Objects of this class are built from three substitution
 * models of NucleicAlphabets. No model is directly accessible. </p>
 */
class TripletSubstitutionModel :
  public WordSubstitutionModel
{
public:
  /**
   * @brief Build a new TripletSubstitutionModel object from
   * a pointer to a NucleotideSubstitutionModel.
   *
   * @param palph pointer to a CodonAlphabet
   * @param pmod  pointer to the NucleotideSubstitutionModel to be used
   *       in the three positions. It is owned by the instance.
   */
  TripletSubstitutionModel(
      std::shared_ptr<const CodonAlphabet> palph,
      std::unique_ptr<NucleotideSubstitutionModelInterface> pmod);

  /**
   * @brief Build a new TripletSubstitutionModel object
   * from three pointers to NucleotideSubstitutionModels.
   *
   * @param palph pointer to a CodonAlphabet
   * @param pmod1, pmod2, pmod3 pointers to the
   *   NucleotideSubstitutionModels to use in the three positions.
   */
  TripletSubstitutionModel(
      std::shared_ptr<const CodonAlphabet> palph,
      std::unique_ptr<NucleotideSubstitutionModelInterface> pmod1,
      std::unique_ptr<NucleotideSubstitutionModelInterface> pmod2,
      std::unique_ptr<NucleotideSubstitutionModelInterface> pmod3);

  virtual ~TripletSubstitutionModel() {}

  TripletSubstitutionModel* clone() const override { return new TripletSubstitutionModel(*this);}

public:
  std::string getName() const override;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_TRIPLETSUBSTITUTIONMODEL_H
