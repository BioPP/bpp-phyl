//
// File: SetOfAbstractPhyloLikelihood.h
// Created by: Laurent Guéguen
// Created on: jeudi 8 octobre 2015, à 14h 33
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

#ifndef _SET_OF_ABSTRACT_PHYLOLIKELIHOOD_H_
#define _SET_OF_ABSTRACT_PHYLOLIKELIHOOD_H_


#include "PhyloLikelihood.h"
#include "AbstractPhyloLikelihood.h"
#include "PhyloLikelihoodContainer.h"

namespace bpp
{

    /**
     * @brief The SetOfAbstractPhyloLikelihood class, to manage a
     * subset of AbstractPhyloLikelihoods from a given
     * PhyloLikelihoodContainer
     *
     */
    
  class SetOfAbstractPhyloLikelihood:
    virtual public AbstractPhyloLikelihood,
    public AbstractParametrizable
    {
    protected:

      /**
       * @brief pointer to a  PhyloLikelihoodContainer
       *
       */
      
      PhyloLikelihoodContainer* pPhyloCont_;
      
      /**
       * @brief vector of AbstractPhyloLikelihood numbers
       *
       */

      std::vector<size_t> nPhylo_;

    public:
      SetOfAbstractPhyloLikelihood(Context& context, PhyloLikelihoodContainer* pC, const std::string& prefix = "");

      ~SetOfAbstractPhyloLikelihood() {}

      SetOfAbstractPhyloLikelihood(const SetOfAbstractPhyloLikelihood& sd);
        
    public:

      PhyloLikelihoodContainer* getPhyloContainer()
      {
        return pPhyloCont_;
      }

      const PhyloLikelihoodContainer* getPhyloContainer() const
      {
        return pPhyloCont_;
      }
      
      const std::vector<size_t>& getNumbersOfPhyloLikelihoods() const
      {
        return nPhylo_;
      }

      /**
       *
       * @brief adds a PhyloLikelihood already stored in the
       * PhyloLikelihoodContainer, iff it is an
       * AbstractPhyloLikelihood.
       *
       * @return if the PhyloLikelihood has been added.
       */

      virtual bool addPhyloLikelihood(size_t nPhyl);

      /**
       *
       * @brief adds all PhyloLikelihoods already stored in the
       * PhyloLikelihoodContainer, iff their type fit.
       *
       */

      void addAllPhyloLikelihoods();

      /**
       *
       * @name The AbstractPhyloLikelihood storage.
       *
       * @{
       */
      
      virtual const AbstractPhyloLikelihood* getAbstractPhyloLikelihood(size_t nPhyl) const
      {
        return dynamic_cast<const AbstractPhyloLikelihood*>((*pPhyloCont_)[nPhyl]);
      }
      
      
      virtual AbstractPhyloLikelihood* getAbstractPhyloLikelihood(size_t nPhyl)
      {
        return dynamic_cast<AbstractPhyloLikelihood*>((*pPhyloCont_)[nPhyl]);
      }

      virtual const PhyloLikelihood* getPhyloLikelihood(size_t nPhyl) const
      {
        return (*pPhyloCont_)[nPhyl];
      }
      
      
      virtual PhyloLikelihood* getPhyloLikelihood(size_t nPhyl)
      {
        return (*pPhyloCont_)[nPhyl];
      }

      /**
       *
       * @}
       *
       */

      /**
       *
       * @name Inherited from PhyloLikelihood
       *
       * @{
       */
      

      /**
       * @return 'true' is the likelihood function has been initialized.
       */
      
      bool isInitialized() const 
      {
        for (size_t i=0; i<nPhylo_.size(); i++)
          if (!getPhyloLikelihood(nPhylo_[i])->isInitialized())
            return false;
        return true;
      }
      
    public:
      
      virtual void fireParameterChanged(const ParameterList& params)
      {
        for (size_t i=0; i<nPhylo_.size(); i++){
          getAbstractPhyloLikelihood(nPhylo_[i])->matchParametersValues(params);

          // to ensure each phylolikelihood is recomputed, such as in
          // case of total aliasing
          // getAbstractPhyloLikelihood(nPhylo_[i])->update();
        }
        
      }

      double getFirstOrderDerivative(const std::string& variable) const 
      {
        // if (dValues_.find(variable)==dValues_.end())
        //   computeDLogLikelihood_(variable);

        // if (dValues_.find(variable)==dValues_.end() || std::isnan(dValues_[variable]))
        //   dValues_[variable]=-getDLogLikelihood(variable);
        
        // return dValues_[variable];
        return 0;
      }

      double getSecondOrderDerivative(const std::string& variable) const
      {
        // if (d2Values_.find(variable)==d2Values_.end())
        //   computeD2LogLikelihood_(variable);

        // if (d2Values_.find(variable)==d2Values_.end() || std::isnan(d2Values_[variable]))
        //   d2Values_[variable]=-getD2LogLikelihood(variable);
        
        // return d2Values_[variable];
        return 0;
      }

      double getSecondOrderDerivative(const std::string& variable1, const std::string& variable2) const
      { return 0; } // Not implemented for now.
      
      /**
       * @name Retrieve some particular parameters subsets.
       *
       * @{
       */
    
      /**
       * @brief Get the branch lengths parameters.
       *
       * @return A ParameterList with all branch lengths.
       */

      ParameterList getBranchLengthParameters() const;
      
      /**
      * @brief Get the parameters associated to substitution model(s).
      * 
      * @return A ParameterList.
      */

      ParameterList getSubstitutionModelParameters() const;

      /**
       * @brief Get the parameters associated to the rate distribution(s).
       *
       * @return A ParameterList.
       */

      ParameterList getRateDistributionParameters() const;

      /**
       * @brief Get the parameters associated to the root frequencies(s).
       *
       * @return A ParameterList.
       */

      ParameterList getRootFrequenciesParameters() const;

      /**
       * @brief All non derivable parameters.
       *
       * It is the union of the non derivable parameters from all members.
       *
       * @return A ParameterList.
       */

      virtual ParameterList getNonDerivableParameters() const;


      /** @} */

      
    };

} //end of namespace bpp.

#endif  //_SET_OF_ABSTRACT_PHYLOLIKELIHOOD_H_

