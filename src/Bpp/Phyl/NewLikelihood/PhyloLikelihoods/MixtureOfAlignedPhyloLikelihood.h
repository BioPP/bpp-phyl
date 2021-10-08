//
// File: MixtureOfAlignedPhyloLikelihood.h
// Authors:
//   Laurent GuÃ©guen
// Created: mercredi 7 octobre 2015, Ã  14h 03
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

#ifndef BPP_PHYL_NEWLIKELIHOOD_PHYLOLIKELIHOODS_MIXTUREOFALIGNEDPHYLOLIKELIHOOD_H
#define BPP_PHYL_NEWLIKELIHOOD_PHYLOLIKELIHOODS_MIXTUREOFALIGNEDPHYLOLIKELIHOOD_H


#include "SetOfAlignedPhyloLikelihood.h"

// From SeqLib:
#include <Bpp/Seq/Container/AlignedValuesContainer.h>

#include "../DataFlow/Simplex_DF.h"

namespace bpp
{
/**
 * @brief Likelihood framework based on a mixture of aligned likelihoods
 *
 * The resulting likelihood is the mean value of
 * the AlignedPhyloLikelihoods, ponderated with parametrized probabilities
 * (through a Simplex).
 *
 */

class MixtureOfAlignedPhyloLikelihood :
  public SetOfAlignedPhyloLikelihood
{
private:
  /**
   * DF simplex
   *
   */

  std::shared_ptr<ConfiguredSimplex> simplex_;

  /**
   * Aligned LikelihoodCalculation to store DF nodes
   */

  mutable std::shared_ptr<AlignedLikelihoodCalculation> likCal_;

public:
  MixtureOfAlignedPhyloLikelihood(Context& context, std::shared_ptr<PhyloLikelihoodContainer> pC, const std::vector<size_t>& nPhylo, bool inCollection = true);

  MixtureOfAlignedPhyloLikelihood(const MixtureOfAlignedPhyloLikelihood& mlc);

  ~MixtureOfAlignedPhyloLikelihood() {}

  MixtureOfAlignedPhyloLikelihood* clone() const
  {
    return new MixtureOfAlignedPhyloLikelihood(*this);
  }

protected:
  void fireParameterChanged(const ParameterList& parameters);

public:
  std::shared_ptr<LikelihoodCalculation> getLikelihoodCalculation () const
  {
    return likCal_;
  }

  std::shared_ptr<AlignedLikelihoodCalculation> getAlignedLikelihoodCalculation () const
  {
    return likCal_;
  }

  /**
   * @brief Get the probabilities of the simplex
   *
   */

  Vdouble getPhyloProbabilities() const;

  /**
   * @brief Get the probability of a phylolikelihood
   *
   */

  double getPhyloProb(size_t index) const;

  /**
   * @brief Set the probabilities of the simplex
   *
   */

  void setPhyloProb(Simplex const& simplex);
};
} // end of namespace bpp.
#endif // BPP_PHYL_NEWLIKELIHOOD_PHYLOLIKELIHOODS_MIXTUREOFALIGNEDPHYLOLIKELIHOOD_H
