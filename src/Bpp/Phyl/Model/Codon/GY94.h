//
// File: GY94.h
// Authors:
//   Laurent Gueguen
// Created: 2009-07-08 00:00:00
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

#ifndef BPP_PHYL_MODEL_CODON_GY94_H
#define BPP_PHYL_MODEL_CODON_GY94_H

#include <Bpp/Seq/AlphabetIndex/GranthamAAChemicalDistance.h>

#include "../AbstractBiblioSubstitutionModel.h"
#include "CodonDistanceFrequenciesSubstitutionModel.h"

namespace bpp
{
/**
 * @brief The Goldman and Yang (1994) substitution model for codons.
 * @author Laurent GuÃÂ©guen
 *
 * This model has one rate of transitions and one rate of
 * transversion. It also allows distinct equilibrium frequencies
 * between codons. A multiplicative factor accounts for the selective
 * restraints at the amino acid level. This factor applies on the
 * distance @f$d@f$ between amino acids given by Grantham (1974).
 *
 * For codons @f$i=i_1i_2i_3@f$ and @f$j=j_1j_2j_3@f$, the generator
 * term @f$Q_{ij} (i \neq j)@f$ is:
 *
 * 0 if 2 or 3 of the pair @f$(i_1,j_1) (i_2,j_2) (i_3,j_3) @f$ are
 * different.
 *
 * @f$\mu \pi_j \exp(-d_{aa_i,aa_j}/V)@f$ if exactly 1 of the pairs
 * @f$(i_1,j_1) (i_2,j_2) (i_3,j_3) @f$ is different, and that
 * difference is a transversion.
 *
 * @f$\mu \kappa \pi_j \exp(-d_{aa_i,aa_j}/V)@f$ if exactly 1 of the
 * pairs @f$(i_1,j_1) (i_2,j_2) (i_3,j_3) @f$ is different, and that
 * difference is a transition.
 *
 * @f$\mu@f$ is a normalization factor.
 *
 * This model includes 2 parameters (@f$\kappa@f$ and @f$V@f$). The
 * codon frequencies @f$\pi_j@f$ are either observed or infered.
 *
 * Reference:
 * - Goldman N. and Yang Z. (1994), _Molecular Biology And Evolution_ 11(5) 725--736.
 */
class GY94 :
  public AbstractBiblioSubstitutionModel,
  public virtual ReversibleSubstitutionModelInterface
{

private:

  std::shared_ptr<const GranthamAAChemicalDistance> gacd_;
  std::unique_ptr<CodonDistanceFrequenciesSubstitutionModel> pmodel_;

public:

  GY94(
      std::shared_ptr<const GeneticCode> gc,
      std::unique_ptr<CodonFrequencySetInterface> codonFreqs);

  virtual ~GY94();

  GY94(const GY94& gy94);

  GY94& operator=(const GY94& gy94);

  GY94* clone() const override { return new GY94(*this); }

public:

  std::string getName() const override { return "GY94"; }

  const SubstitutionModelInterface& substitutionModel() const override { return *pmodel_; }

  std::shared_ptr<const GeneticCode> getGeneticCode() const { return pmodel_->getGeneticCode(); }

  double getCodonsMulRate(size_t i, size_t j) const { return pmodel_->getCodonsMulRate(i, j); }

protected:
  
  SubstitutionModelInterface& substitutionModel_() override { return *pmodel_; }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_GY94_H
