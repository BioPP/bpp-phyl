//
// File: FormulaOfPhyloLikelihood.h
// Created by: Laurent Guéguen
// Created on: jeudi 8 décembre 2016, à 10h 35
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

#ifndef _FORMULA_OFPHYLOLIKELIHOOD_H_
#define _FORMULA_OFPHYLOLIKELIHOOD_H_


#include "SetOfAbstractPhyloLikelihood.h"

#include <Bpp/Numeric/Function/Operators/ComputationTree.h>

namespace bpp
{

  /**
   * @brief The FormulaOfPhyloLikelihood class, for phylogenetic
   * likelihood on several independent data.
   *
   * WARNING: This formula applies on the log-likelihoods (ie getValues())
   *
   */
  
  class FormulaOfPhyloLikelihood:
    public SetOfAbstractPhyloLikelihood
  {
  private:
    std::unique_ptr<ComputationTree> compTree_;

    std::shared_ptr<LikelihoodCalculation> likCal_;
    
  public:
    
    FormulaOfPhyloLikelihood(Context& context, std::shared_ptr<PhyloLikelihoodContainer> pC, const std::string& formula, bool inCollection = true);

    ~FormulaOfPhyloLikelihood() {}
    
    FormulaOfPhyloLikelihood* clone() const
    {
      return new FormulaOfPhyloLikelihood(*this);
    }
    
    FormulaOfPhyloLikelihood(const FormulaOfPhyloLikelihood& sd);

  public:

    /**
     * @ input
     *
     */
    
    void readFormula(const std::string& formula, bool inCollection = true);

    /**
     * @ output
     *
     */
    
    std::string output() const;

    /**
     * @name The likelihood functions.
     *
     */
    
    /**
     * @brief Get the logarithm of the likelihood for the whole dataset.
     *
     * @return The logarithm of the likelihood of the dataset.
     */

    
    std::shared_ptr<LikelihoodCalculation> getLikelihoodCalculation() const
    {
      return likCal_;
    }
    
  private:
    
    /**
     * @brief Build the LikelihoodNode from the computation Tree
     *
     */

    ValueRef<DataLik> makeLikelihoods()
    {
      return makeLikelihoodsFromOperator(compTree_->getRoot());
    }

    /**
     * @brief Build the LikelihoodNode from a node of the computation Tree
     *
     */

    ValueRef<DataLik> makeLikelihoodsFromOperator(std::shared_ptr<Operator> op);
      
      
    /** @} */
    
  };

} //end of namespace bpp.

#endif  //_PRODUCT_OFPHYLOLIKELIHOOD_H_

