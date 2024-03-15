// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_PRODUCTOFALIGNEDPHYLOLIKELIHOOD_H
#define BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_PRODUCTOFALIGNEDPHYLOLIKELIHOOD_H



// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/AlignmentData.h>

#include "AlignedPhyloLikelihoodSet.h"

namespace bpp
{
/**
 * @brief The AlignedPhyloLikelihoodProduct class, for phylogenetic
 * likelihood on several independent data.
 */
class AlignedPhyloLikelihoodProduct :
  public AbstractAlignedPhyloLikelihoodSet
{
  /**
   * Aligned LikelihoodCalculation to store DF nodes
   */
  mutable std::shared_ptr<AlignedLikelihoodCalculation> likCal_;

public:
  AlignedPhyloLikelihoodProduct(
      Context& context,
      std::shared_ptr<PhyloLikelihoodContainer> pC,
      bool inCollection = true);

  AlignedPhyloLikelihoodProduct(
      Context& context,
      std::shared_ptr<PhyloLikelihoodContainer> pC,
      const std::vector<size_t>& nPhylo,
      bool inCollection = true);

  virtual ~AlignedPhyloLikelihoodProduct() {}

protected:

  AlignedPhyloLikelihoodProduct(const AlignedPhyloLikelihoodProduct& sd) :
    AbstractPhyloLikelihood(sd),
    AbstractParametrizable(""),
    AbstractPhyloLikelihoodSet(sd),
    AbstractAlignedPhyloLikelihood(sd),
    AbstractAlignedPhyloLikelihoodSet(sd),
    likCal_(sd.likCal_)
  {}

  AlignedPhyloLikelihoodProduct& operator=(const AlignedPhyloLikelihoodProduct& sd)
  {
    AbstractAlignedPhyloLikelihoodSet::operator=(sd);
    likCal_ = sd.likCal_;
    return *this;
  }

  AlignedPhyloLikelihoodProduct* clone() const
  {
    return new AlignedPhyloLikelihoodProduct(*this);
  }

public:
  LikelihoodCalculation& likelihoodCalculation () const
  {
    return *likCal_;
  }

  std::shared_ptr<LikelihoodCalculation> getLikelihoodCalculation () const
  {
    return likCal_;
  }

  AlignedLikelihoodCalculation& alignedLikelihoodCalculation () const
  {
    return *likCal_;
  }

  std::shared_ptr<AlignedLikelihoodCalculation> getAlignedLikelihoodCalculation () const
  {
    return likCal_;
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_PRODUCTOFALIGNEDPHYLOLIKELIHOOD_H
