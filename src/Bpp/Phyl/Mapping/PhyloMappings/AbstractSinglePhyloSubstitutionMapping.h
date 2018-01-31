//
// File: AbstractSinglePhyloSubstitutionMapping.h
// Created by: Laurent Guéguen
// Created on: dimanche 3 décembre 2017, à 13h 56
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

   This software is a computer program whose purpose is to provide classes
   for phylogenetic data analysis.

   This software is governed by the CeCILL  license under French law and
   abiding by the rules of distribution of free software.  You can  use,
   modify and/ or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL
   "http://www.cecill.info".

   As a counterpart to the access to the source code and  rights to copy,
   modify and redistribute granted by the license, users are provided only
   with a limited warranty  and the software's author,  the holder of the
   economic rights,  and the successive licensors  have only  limited
   liability.

   In this respect, the user's attention is drawn to the risks associated
   with loading,  using,  modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean  that it is complicated to manipulate,  and  that  also
   therefore means  that it is reserved for developers  and  experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards their
   requirements in conditions enabling the security of their systems and/or
   data to be ensured and,  more generally, to use and operate it in the
   same conditions as regards security.

   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
 */

#ifndef _ABSTRACT_SINGLE_PHYLO_SUBSTITUTION_MAPPING_H_
#define _ABSTRACT_SINGLE_PHYLO_SUBSTITUTION_MAPPING_H_

#include "../../Model/BranchedModelSet.h"
#include "../ProbabilisticSubstitutionMapping.h"
#include "PhyloSubstitutionMapping.h"
#include "../../Tree/PhyloTree.h"

#include <Bpp/Numeric/ParametrizableCollection.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>

namespace bpp
{
/**
 * @brief The AbstractSinglePhyloSubstitutionMapping class: substitution
 * mapping linked with a Single Process PhyloLikelihood
 *
 */
  
  struct ModelBranch 
  {
    TransitionModel* pMod_;
  };


  class AbstractSinglePhyloSubstitutionMapping :
    virtual public BranchedModelSet,
    virtual public PhyloSubstitutionMapping,
    public AssociationTreeGlobalGraphObserver<uint, ModelBranch>
  {
  public:

    typedef  AssociationTreeGlobalGraphObserver<uint, ModelBranch> modelTree;
    
  private:

    const SubstitutionRegister* pReg_;

  protected:
    
    std::unique_ptr<ProbabilisticSubstitutionMapping> counts_;
    std::unique_ptr<ProbabilisticSubstitutionMapping> factors_;

  private:
    
    /**
     * @brief A collection of Transition Models
     */

    ParametrizableCollection<TransitionModel> modelColl_;

    /**
     *
     * @brief a map <model index, vector of branch ids>
     *
     */
    
    std::map<size_t, std::vector<uint> > mModBrid_;
    
  public:
    AbstractSinglePhyloSubstitutionMapping(std::shared_ptr<TreeGlobalGraph> graph, const SubstitutionRegister& reg) :
      modelTree(graph), pReg_(&reg), counts_(), factors_(), modelColl_(), mModBrid_()
    {}

    AbstractSinglePhyloSubstitutionMapping(const AbstractSinglePhyloSubstitutionMapping& sppm);
    
    AbstractSinglePhyloSubstitutionMapping& operator=(const AbstractSinglePhyloSubstitutionMapping& sppm);

    virtual ~AbstractSinglePhyloSubstitutionMapping() {}
    
    // AbstractSinglePhyloSubstitutionMapping* clone() const { return new AbstractSinglePhyloSubstitutionMapping(*this); }
    
    /*
     * @brief From BranchedModelSet
     *
     * @{
     *
     */
    
    TransitionModel* getModelForBranch(uint branchId)
    {
      return (*getEdge(branchId)).pMod_;
    }

    const TransitionModel* getModelForBranch(uint branchId) const
    {
      return (*getEdge(branchId)).pMod_;
    }

    TransitionModel* getModel(unsigned int branchId, size_t classIndex) const
    {
      return (*getEdge(branchId)).pMod_;
    }

    TransitionModel* getModel(size_t index)
    {
      return modelColl_[index];
    }

    const TransitionModel* getModel(size_t index) const
    {
      return modelColl_[index];
    }

    std::vector<uint> getBranchesWithModel(size_t index) const
    {
      return mModBrid_.at(index);
    }
    
    /*
     * @}
     *
     */

    /*
     * @brief From PhyloSubstitutionMapping.
     *
     * @{
     */
    
    void setRegister(const SubstitutionRegister& reg)
    {
      pReg_=&reg;
    }

    const SubstitutionRegister& getRegister() const
    {
      return *pReg_;
    }
    
    bool matchParametersValues(const ParameterList& nullParams)
    {
      return modelColl_.matchParametersValues(nullParams);
    }

    const ParameterList& getParameters() const
    {
      return modelColl_.getParameters();
    }

    /* 
     *
     * @brief add a Substitition Model in the map, and on all branches
     * with given Ids.
     *
     * @param index the index of the model
     * @param model the model that will be COPIED.
     * @param brIds the Ids of the branches that will carry this model.
     *
     *
     */
    
    void addModel(size_t index, const TransitionModel& model, Vuint brIds);

  };
  
} // end of namespace bpp.

#endif  // _ABSTRACT_SINGLE_PHYLO_SUBSTITUTION_MAPPING_H_
