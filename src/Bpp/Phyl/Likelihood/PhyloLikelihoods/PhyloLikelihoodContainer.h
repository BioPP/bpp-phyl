//
// File: PhyloLikelihoodContainer.h
// Authors:
//   Laurent Guéguen
// Created: mercredi 7 octobre 2015, ÃÂ  22h 34
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

#ifndef BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_PHYLOLIKELIHOODCONTAINER_H
#define BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_PHYLOLIKELIHOODCONTAINER_H



// From bpp-seq:
#include <Bpp/Seq/Container/AlignmentData.h>

#include "PhyloLikelihood.h"
#include "SingleDataPhyloLikelihood.h"

#include "../DataFlow/CollectionNodes.h"

namespace bpp
{
/**
 * @brief The PhyloLikelihoodContainer, owns and assigns numbers to
 * Phylolikelihoods.
 *
 * It owns the PhyloLikelihoods
 */
class PhyloLikelihoodContainer :
  virtual public Clonable
{
private:
  Context& context_;

  std::shared_ptr<CollectionNodes> collectionNodes_;

protected:
  std::map<size_t, std::shared_ptr<PhyloLikelihoodInterface> >  mPhylo_;

public:
  PhyloLikelihoodContainer(
      Context& context, 
      std::shared_ptr<SubstitutionProcessCollection> sColl) :
    context_(context),
    collectionNodes_(std::make_shared<CollectionNodes>(context_, sColl)),
    mPhylo_()
  {}

  PhyloLikelihoodContainer(Context& context, std::shared_ptr<CollectionNodes> sColl) :
    context_(context),
    collectionNodes_(sColl),
    mPhylo_()
  {}

  PhyloLikelihoodContainer* clone() const override
  {
    throw Exception("PhyloLikelihoodContainer::clone should not be called.");
  }

  /**
   * @brief Abstract class destructor
   *
   */
  virtual ~PhyloLikelihoodContainer()
  {}

public:
  /*
   * @brief add a PhyloLikelihood in the container, at a given
   * position.
   *
   * Beware! Takes possession of the PhyloLikelihood through
   *
   */
  void addPhyloLikelihood(size_t pos, std::shared_ptr<PhyloLikelihoodInterface> Ap)
  {
    if (mPhylo_.find(pos) != mPhylo_.end())
      throw Exception("PhyloLikelihoodContainer::addPhylolikelihood: map number already used : " + TextTools::toString(pos));
    mPhylo_[pos] = Ap;
  }

  bool hasPhyloLikelihood(size_t pos) const
  {
    return mPhylo_.find(pos) != mPhylo_.end();
  }

  std::shared_ptr<const PhyloLikelihoodInterface> operator[](size_t pos) const
  {
    auto it = mPhylo_.find(pos);
    return it != mPhylo_.end() ? it->second : nullptr;
  }

  std::shared_ptr<PhyloLikelihoodInterface> operator[](size_t pos)
  {
    auto it = mPhylo_.find(pos);
    return it != mPhylo_.end() ? it->second : nullptr;
  }

  std::shared_ptr<const PhyloLikelihoodInterface> getPhyloLikelihood(size_t pos) const
  {
    auto it = mPhylo_.find(pos);
    return it != mPhylo_.end() ? it->second : nullptr;
  }

  std::shared_ptr<PhyloLikelihoodInterface> getPhyloLikelihood(size_t pos)
  {
    auto it = mPhylo_.find(pos);
    return it != mPhylo_.end() ? it->second : nullptr;
  }

  size_t getSize() const
  {
    return mPhylo_.size();
  }

  std::vector<size_t> getNumbersOfPhyloLikelihoods() const
  {
    std::vector<size_t> vnum;

    for (const auto& it : mPhylo_)
    {
      vnum.push_back(it.first);
    }

    return vnum;
  }

  /**
   *@brief Manage Collection Nodes
   *
   */
  std::shared_ptr<const CollectionNodes> getCollectionNodes() const
  {
    return collectionNodes_;
  }

  std::shared_ptr<CollectionNodes> getCollectionNodes()
  {
    return collectionNodes_;
  }

  /**
   * @brief Set the dataset for which the likelihood must be
   * evaluated, iff the pointed PhyloLikelihood is a SingleDataPhyloLikelihood
   *
   * @param nPhyl The number of the Likelihood.
   * @param sites The data set to use.
   */
  void setData(std::shared_ptr<const AlignmentDataInterface> sites, size_t nPhyl)
  {
    auto it = mPhylo_.find(nPhyl);
    if (it != mPhylo_.end())
    {
      auto sdp = std::dynamic_pointer_cast<SingleDataPhyloLikelihoodInterface>(it->second);
      if (sdp)
        sdp->setData(sites);
    }
  }


  /**
   * @brief Get the dataset for which the likelihood must be evaluated.
   *
   * @return A pointer toward the site container where the sequences are stored.
   */
  std::shared_ptr<const AlignmentDataInterface> getData(size_t nPhyl) const
  {
    const auto it = mPhylo_.find(nPhyl);
    if (it != mPhylo_.end())
    {
      auto sdp = std::dynamic_pointer_cast<SingleDataPhyloLikelihoodInterface>(it->second);
      if (sdp)
        return sdp->getData();
    }
    return nullptr;
  }
  
};
} // end of namespace bpp.
#endif // BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_PHYLOLIKELIHOODCONTAINER_H
