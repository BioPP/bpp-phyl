//
// File: SequencePhyloLikelihood.h
// Created by: Laurent Guéguen
// Created on: mardi 28 avril 2015, à 11h 41
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

#ifndef _SEQUENCEPHYLOLIKELIHOOD_H_
#define _SEQUENCEPHYLOLIKELIHOOD_H_

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/AlignedValuesContainer.h>

#include "PhyloLikelihood.h"
#include "AbstractPhyloLikelihood.h"
#include "SingleDataPhyloLikelihood.h"
#include "../SequenceEvolution.h"

namespace bpp
{
/**
 * @brief PhyloLikelihoods based on Sequence Evolution class, ie for
 * which there is a (set of) processes able to build a sequence.
 *
 */
  
  /**
   * @brief Interface
   *
   */
  
  class SequencePhyloLikelihood :
    public AbstractSingleDataPhyloLikelihood
  {
  protected:
    /**
     * @brief the Sequence Evolution
     *
     **/
    
    SequenceEvolution* seqEvol_;
    
    /**
     * @brief the Sequence Evolution number
     *
     **/

    size_t nSeqEvol_;
      
  public:
    SequencePhyloLikelihood(Context& context, SequenceEvolution& se, size_t nSE = 0, size_t nData = 0) :
      AbstractPhyloLikelihood(context),
      AbstractAlignedPhyloLikelihood(context, 0),
      AbstractSingleDataPhyloLikelihood(context, 0, (se.getSubstitutionProcessNumbers().size()!=0)?se.getSubstitutionProcess(se.getSubstitutionProcessNumbers()[0]).getNumberOfStates():0, nData),
      seqEvol_(&se),
      nSeqEvol_(nSE)
    {
    }

    SequencePhyloLikelihood(const SequencePhyloLikelihood& asd) :
      AbstractPhyloLikelihood(asd),
      AbstractAlignedPhyloLikelihood(asd),
      AbstractSingleDataPhyloLikelihood(asd),
      seqEvol_(asd.seqEvol_),
      nSeqEvol_(asd.nSeqEvol_)
    {
    }
      
    virtual ~SequencePhyloLikelihood() {}

    SequencePhyloLikelihood* clone() const = 0;

  public:
    /**
     * @brief The Sequence Evolution
     *
     */
      
    const SequenceEvolution& getSequenceEvolution() const
    {
      return *seqEvol_;
    }

    const size_t getSequenceEvolutionNumber() const
    {
      return nSeqEvol_;
    }

    /**
     * @name the Likelihood interface
     *
     */
      
    const Alphabet* getAlphabet() const {
      if (seqEvol_->getSubstitutionProcessNumbers().size()==0)
        return NULL;
      else 
        return seqEvol_->getSubstitutionProcess(seqEvol_->getSubstitutionProcessNumbers()[0]).getModel(0,0)->getAlphabet();
    }


  };

  class AbstractSequencePhyloLikelihood :
    public SequencePhyloLikelihood,
    public AbstractParametrizable
  {
  public:
    AbstractSequencePhyloLikelihood(Context& context, SequenceEvolution& se, size_t nSE = 0, size_t nData = 0) :
      AbstractPhyloLikelihood(context),
      AbstractAlignedPhyloLikelihood(context, 0),
      SequencePhyloLikelihood(context, se, nSE, nData),
      AbstractParametrizable("")
    {
      // initialize INDEPENDENT parameters:
      addParameters_(se.getIndependentParameters());
    }

    AbstractSequencePhyloLikelihood(const AbstractSequencePhyloLikelihood& asd) :
      AbstractPhyloLikelihood(asd),
      AbstractAlignedPhyloLikelihood(asd),
      SequencePhyloLikelihood(asd),
      AbstractParametrizable(asd)
    {
    }
      
    virtual ~AbstractSequencePhyloLikelihood() {}

    AbstractSequencePhyloLikelihood* clone() const = 0;

  protected:
    virtual void fireParameterChanged(const ParameterList& parameters)
    {
      seqEvol_->matchParametersValues(parameters);
      update();
    }

  public:

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

    void setData(const AlignedValuesContainer& sites, size_t nData = 0)
    {
      AbstractSingleDataPhyloLikelihood::setData(sites, nData);
      update();
    }
      
    void setNamespace(const std::string& nameSpace)
    {
      seqEvol_->setNamespace(nameSpace);
    }

    ParameterList getNonDerivableParameters() const
    {
      return seqEvol_->getNonDerivableParameters();
    }

    ParameterList getBranchLengthParameters() const
    {
      return seqEvol_->getBranchLengthParameters(true);
    }

    ParameterList getRootFrequenciesParameters() const
    {
      return seqEvol_->getRootFrequenciesParameters(true);
    }

    ParameterList getRateDistributionParameters() const
    {
      return seqEvol_->getRateDistributionParameters(true);
    }

    ParameterList getSubstitutionModelParameters() const
    {
      return seqEvol_->getSubstitutionModelParameters(true);
    }

    ParameterList getSubstitutionProcessParameters() const
    {
      return seqEvol_->getSubstitutionProcessParameters(true);
    }


  };
      
} //end of namespace bpp.

#endif  //_SEQUENCEPHYLOLIKELIHOOD_H_

