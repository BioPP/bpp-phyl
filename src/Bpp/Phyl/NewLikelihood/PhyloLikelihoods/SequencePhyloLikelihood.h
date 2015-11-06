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
#include <Bpp/Seq/Container/SiteContainer.h>

#include "PhyloLikelihood.h"
#include "SingleDataPhyloLikelihood.h"
#include "../SequenceEvolution.h"

namespace bpp
{
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
      SequencePhyloLikelihood(SequenceEvolution& se, size_t nSE = 0, size_t nData = 0) :
        AbstractPhyloLikelihood(),
        AbstractAlignedPhyloLikelihood(0),
        AbstractSingleDataPhyloLikelihood(0, (se.getSubstitutionProcessNumbers().size()!=0)?se.getSubstitutionProcess(se.getSubstitutionProcessNumbers()[0]).getNumberOfStates():0, nData),
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
      
      SequencePhyloLikelihood& operator=(const SequencePhyloLikelihood& asd)
      {
        AbstractSingleDataPhyloLikelihood::operator=(asd);

        seqEvol_=asd.seqEvol_;
        nSeqEvol_=asd.nSeqEvol_;

        return *this;
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
          return seqEvol_->getSubstitutionProcess(seqEvol_->getSubstitutionProcessNumbers()[0]).getSubstitutionModel(0,0).getAlphabet();
      }


    };

    class AbstractSequencePhyloLikelihood :
      public SequencePhyloLikelihood,
      public AbstractParametrizable
    {
    public:
      AbstractSequencePhyloLikelihood(SequenceEvolution& se, size_t nSE = 0, size_t nData = 0) :
        AbstractPhyloLikelihood(),
        AbstractAlignedPhyloLikelihood(0),
        SequencePhyloLikelihood(se, nSE, nData),
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
      
      AbstractSequencePhyloLikelihood& operator=(const AbstractSequencePhyloLikelihood& asd)
      {
        SequencePhyloLikelihood::operator=(asd);

        AbstractParametrizable::operator=(asd);
        
        return *this;
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
      double getFirstOrderDerivative(const std::string& variable) const throw (Exception)
      {
        if (!hasParameter(variable))
          throw ParameterNotFoundException("SequencePhyloLikelihood::getFirstOrderDerivative().", variable);
//        if (!hasDerivableParameter(variable))
        {
          throw Exception("SequencePhyloLikelihood::Derivative is not implemented for " + variable + " parameter.");
        }

        computeDLogLikelihood_(variable);
        return -getDLogLikelihood();
      }

      double getSecondOrderDerivative(const std::string& variable) const throw (Exception)
      {
        if (!hasParameter(variable))
          throw ParameterNotFoundException("AbstractPhyloLikelihood::getSecondOrderDerivative().", variable);
//        if (!hasDerivableParameter(variable))
        {
          throw Exception("Derivative is not implemented for " + variable + " parameter.");
        }
        
        computeD2LogLikelihood_(variable);
        return -getD2LogLikelihood();
      }

      double getSecondOrderDerivative(const std::string& variable1, const std::string& variable2) const throw (Exception) { return 0; } // Not implemented for now.


      void setData(const SiteContainer& sites, size_t nData = 0)
      {
        AbstractSingleDataPhyloLikelihood::setData(sites, nData);
        update();
      }
      
      void setNamespace(const std::string& nameSpace)
      {
        seqEvol_->setNamespace(nameSpace);
      }

      bool hasDerivableParameter(const std::string& variable) const
      {
        return seqEvol_->hasDerivableParameter(variable);
      }
      
      ParameterList getDerivableParameters() const
      {
        return seqEvol_->getDerivableParameters();
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

