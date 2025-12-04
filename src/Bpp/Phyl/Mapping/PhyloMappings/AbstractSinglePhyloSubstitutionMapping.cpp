// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "AbstractSinglePhyloSubstitutionMapping.h"

using namespace bpp;
using namespace std;


AbstractSinglePhyloSubstitutionMapping::AbstractSinglePhyloSubstitutionMapping(const AbstractSinglePhyloSubstitutionMapping& sppm) :
  modelTree(sppm),
  pReg_(sppm.pReg_),
  weights_(sppm.weights_),
  distances_(sppm.distances_),
  counts_(sppm.counts_ ? sppm.counts_->clone() : 0),
  factors_(sppm.factors_ ? sppm.factors_->clone() : 0),
  modelColl_(sppm.modelColl_),
  mModBrid_(sppm.mModBrid_)
{
  // assign correct addresses
  unique_ptr<modelTree::EdgeIterator> nIT = allEdgesIterator();

  vector<size_t> keys = modelColl_.keys();

  for ( ; !nIT->end(); nIT->next())
  {
    (**nIT)->pMod_ = 0;
    unsigned int brid = getEdgeIndex(**nIT);

    auto sm = sppm.getEdge(brid)->pMod_;

    if (sm)
    {
      for (auto& k:keys)
      {
        if (sppm.modelColl_[k] == sm)
        {
          (**nIT)->pMod_ = modelColl_[k];
          break;
        }
      }

      if (!(**nIT)->pMod_)
        throw Exception("AbstractSinglePhyloSubstitutionMapping::AbstractSinglePhyloSubstitutionMapping: unable to find model for branch " + TextTools::toString(brid));
    }
  }
}


AbstractSinglePhyloSubstitutionMapping& AbstractSinglePhyloSubstitutionMapping::operator=(const AbstractSinglePhyloSubstitutionMapping& sppm)
{
  modelTree::operator=(sppm);
  pReg_ = sppm.pReg_;
  weights_ = sppm.weights_;
  distances_ = sppm.distances_;

  counts_.reset(sppm.counts_ ? sppm.counts_->clone() : 0);
  factors_.reset(sppm.factors_ ? sppm.factors_->clone() : 0);
  modelColl_ = sppm.modelColl_;

  mModBrid_ = sppm.mModBrid_;

  // assign correct addresses
  unique_ptr<modelTree::EdgeIterator> nIT = allEdgesIterator();
  vector<size_t> keys = modelColl_.keys();

  for ( ; !nIT->end(); nIT->next())
  {
    (**nIT)->pMod_ = 0;

    unsigned int brid = getEdgeIndex(**nIT);
    auto sm = sppm.getEdge(brid)->pMod_;

    if (sm)
    {
      for (auto& k:keys)
      {
        if (sppm.modelColl_[k] == sm)
        {
          (**nIT)->pMod_ = modelColl_[k];
          break;
        }
      }

      if (!(**nIT)->pMod_)
        throw Exception("AbstractSinglePhyloSubstitutionMapping::AbstractSinglePhyloSubstitutionMapping: unable to find model for branch " + TextTools::toString(brid));
    }
  }

  return *this;
}


void AbstractSinglePhyloSubstitutionMapping::addModel(size_t index, const TransitionModelInterface& model, Vuint brIds)
{
  modelColl_.addObject(std::shared_ptr<TransitionModelInterface>(model.clone()), index);

  auto tm = modelColl_[index];
  mModBrid_[index] = brIds;

  for (auto& id : brIds)
  {
    auto mb = make_shared<ModelBranch>();
    mb->pMod_ = tm;

    associateEdge(mb, id);
  }
}
