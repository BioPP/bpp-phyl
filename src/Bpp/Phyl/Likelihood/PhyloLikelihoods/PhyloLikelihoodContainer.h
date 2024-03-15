// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

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
