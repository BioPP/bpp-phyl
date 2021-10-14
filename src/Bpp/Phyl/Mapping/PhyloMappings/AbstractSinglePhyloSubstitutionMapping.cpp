//
// File: AbstractSinglePhyloSubstitutionMapping.cpp
// Authors:
//   Laurent GuÃÂ©guen
// Created: mercredi 6 dÃÂ©cembre 2017, ÃÂ  22h 34
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
  // assign correct adresses
  unique_ptr<modelTree::EdgeIterator> nIT = allEdgesIterator();

  vector<size_t> keys = modelColl_.keys();

  for ( ; !nIT->end(); nIT->next())
  {
    (**nIT)->pMod_ = 0;
    uint brid = getEdgeIndex(**nIT);

    TransitionModel* sm = sppm.getEdge(brid)->pMod_;

    if (sm)
    {
      for (auto& k:keys)
      {
        if (sppm.modelColl_[k].get() == sm)
        {
          (**nIT)->pMod_ = modelColl_[k].get();
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

  // assign correct adresses
  unique_ptr<modelTree::EdgeIterator> nIT = allEdgesIterator();
  vector<size_t> keys = modelColl_.keys();

  for ( ; !nIT->end(); nIT->next())
  {
    (**nIT)->pMod_ = 0;

    uint brid = getEdgeIndex(**nIT);
    TransitionModel* sm = sppm.getEdge(brid)->pMod_;

    if (sm)
    {
      for (auto& k:keys)
      {
        if (sppm.modelColl_[k].get() == sm)
        {
          (**nIT)->pMod_ = modelColl_[k].get();
          break;
        }
      }

      if (!(**nIT)->pMod_)
        throw Exception("AbstractSinglePhyloSubstitutionMapping::AbstractSinglePhyloSubstitutionMapping: unable to find model for branch " + TextTools::toString(brid));
    }
  }

  return *this;
}


void AbstractSinglePhyloSubstitutionMapping::addModel(size_t index, const TransitionModel& model, Vuint brIds)
{
  modelColl_.addObject(std::shared_ptr<TransitionModel>(model.clone()), index);

  TransitionModel* tm = modelColl_[index].get();
  mModBrid_[index] = brIds;

  for (auto& id : brIds)
  {
    shared_ptr<ModelBranch> mb(new ModelBranch);
    mb->pMod_ = tm;

    associateEdge(mb, id);
  }
}
