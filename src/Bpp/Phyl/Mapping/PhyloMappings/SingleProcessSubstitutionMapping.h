//
// File: SingleProcessSubstitutionMapping.h
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

#ifndef _SINGLE_PROCESS_SUBSTITUTION_MAPPING_H_
#define _SINGLE_PROCESS_SUBSTITUTION_MAPPING_H_

#include "../../NewLikelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h"
#include "../ProbabilisticSubstitutionMapping.h"
#include "../SubstitutionMappingTools.h"

#include "AbstractSinglePhyloSubstitutionMapping.h"

#include <Bpp/Numeric/Matrix/MatrixTools.h>

namespace bpp
{
/**
 * @brief The SingleProcessSubstitutionMapping class: substitution
 * mapping linked with a SingleProcessPhyloLikelihood
 *
 */

  class SingleProcessSubstitutionMapping :
    public AbstractSinglePhyloSubstitutionMapping
  {
  private:
    SingleProcessPhyloLikelihood* pSPP_;

    /*
     * @brief Set the models of the BranchedModelSet to the adhoc
     * branches, for normalization.
     *
     */
    
    void setBranchedModelSet_();
    
  public:
    SingleProcessSubstitutionMapping(SingleProcessPhyloLikelihood& spp, SubstitutionRegister& reg, std::shared_ptr<const AlphabetIndex2> weights, std::shared_ptr<const AlphabetIndex2> distances, double threshold = -1) :
      AbstractSinglePhyloSubstitutionMapping(spp.getTree().getGraph(), reg, weights, distances),
      pSPP_(&spp)
    {
      computeCounts(threshold);
      setBranchedModelSet_();
    }

    SingleProcessSubstitutionMapping(const SingleProcessSubstitutionMapping& sppm) :
      AbstractSinglePhyloSubstitutionMapping(sppm),
      pSPP_(sppm.pSPP_)
    {
    }

    SingleProcessSubstitutionMapping& operator=(const SingleProcessSubstitutionMapping& sppm)
    {
      AbstractSinglePhyloSubstitutionMapping::operator=(sppm);
      pSPP_ = sppm.pSPP_;
      return *this;
    }

    virtual ~SingleProcessSubstitutionMapping() {}
    
    SingleProcessSubstitutionMapping* clone() const { return new SingleProcessSubstitutionMapping(*this); }

    /*
     * @brief ComputeCounts
     *
     */

    void computeCounts(double threshold = -1)
    {
      counts_.reset(SubstitutionMappingTools::computeCounts(getLikelihoodCalculation(),
                                                            getRegister(),
                                                            getWeights(),
                                                            getDistances(),
                                                            threshold));
    }
    
    void computeNormalizations(const ParameterList& nullParams);

    /*
     * @brief Return the tree of counts
     *
     */

    ProbabilisticSubstitutionMapping& getCounts()
    {
      return *counts_;
    }

    const ProbabilisticSubstitutionMapping& getCounts() const
    {
      return *counts_;
    }

    size_t getNumberOfModels() const
    {
      return pSPP_->getSubstitutionProcess().getNumberOfModels();
    }

    std::vector<size_t> getModelNumbers() const
    {
      return pSPP_->getSubstitutionProcess().getModelNumbers();
    }

        
    /*
     * @brief Return the tree of factors
     *
     */
    
    bool normalizationsPerformed() const
    {
      return factors_!=0;
    }
    
    ProbabilisticSubstitutionMapping& getNormalizations()
    {
      return *factors_;
    }

    const ProbabilisticSubstitutionMapping& getNormalizations() const
    {
      return *factors_;
    }

    RecursiveLikelihoodTreeCalculation& getLikelihoodCalculation()
    {
      return *(dynamic_cast<RecursiveLikelihoodTreeCalculation*>(pSPP_->getLikelihoodCalculation()));
    }

    const RecursiveLikelihoodTreeCalculation& getLikelihoodCalculation() const
    {
      return *(dynamic_cast<const RecursiveLikelihoodTreeCalculation*>(pSPP_->getLikelihoodCalculation()));
    }

  };
  
} // end of namespace bpp.

#endif  // _SINGLE_PROCESS_SUBSTITUTION_MAPPING_H_
