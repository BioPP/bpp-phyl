//
// File: MultiProcessSequenceEvolution.h
// Created by: Laurent Guéguen
// Created on: jeudi 11 juillet 2013, à 21h 51
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

#ifndef _MULTI_PROCESS_SEQUENCE_EVOLUTION_H_
#define _MULTI_PROCESS_SEQUENCE_EVOLUTION_H_

#include "SubstitutionProcess.h"
#include "SubstitutionProcessCollection.h"
#include "SequenceEvolution.h"

#include <Bpp/Numeric/AbstractParametrizable.h>

// From SeqLib:
#include <Bpp/Seq/Container/SiteContainer.h>

using namespace std;

namespace bpp
{
/**
 * @brief Partial implementation of multiple processes of sequences.
 *
 */
  
    class MultiProcessSequenceEvolution :
      virtual public SequenceEvolution,
      public AbstractParameterAliasable
    {
    protected:
      SubstitutionProcessCollection* processColl_;

      /*
       * @brief the vector of the substitution process numbers, as
       * they are used in this order.
       *
       */
    
      std::vector<size_t> nProc_;

    public:
      MultiProcessSequenceEvolution(
        SubstitutionProcessCollection* processColl,
        std::vector<size_t> nProc,
        const std::string& prefix = "");

      MultiProcessSequenceEvolution(const MultiProcessSequenceEvolution& lik) :
        AbstractParameterAliasable(lik),
        processColl_(lik.processColl_),
        nProc_(lik.nProc_)
      {
      }

      MultiProcessSequenceEvolution& operator=(const MultiProcessSequenceEvolution& lik)
      {
        AbstractParameterAliasable::operator=(lik);
        
        processColl_=lik.processColl_;
        nProc_=lik.nProc_;
        
        return *this;
      }

    public:
      
      /**
       * @brief The collection
       *
       */

      const SubstitutionProcessCollection& getCollection() const { return *processColl_; }

      SubstitutionProcessCollection& getCollection() { return *processColl_; }

      /**
       * @brief Return the number of process used for computation.
       */
  

      size_t getNumberOfSubstitutionProcess() const { return nProc_.size(); }

      /**
       * @brief Return the SubstitutionProcess of a given index
       * position (in nProc_ vector).
       *
       */
            
      const SubstitutionProcess& getSubstitutionProcess(size_t number) const
      {
        return processColl_->getSubstitutionProcess(number);
      }

      const std::vector<size_t>& getSubstitutionProcessNumbers() const
      {
        return nProc_;
      }  

      bool hasDerivableParameter(const std::string& name) const;

      ParameterList getSubstitutionProcessParameters(bool independent) const;

      ParameterList getSubstitutionModelParameters(bool independent) const;

      ParameterList getRateDistributionParameters(bool independent) const;

      ParameterList getRootFrequenciesParameters(bool independent) const;

      ParameterList getBranchLengthParameters(bool independent) const;

      virtual ParameterList getDerivableParameters() const
      {
        // patch, to be fixed properly later
        return ParameterList();

        return getBranchLengthParameters(true);
      }

      virtual ParameterList getNonDerivableParameters() const;
  

      virtual void fireParameterChanged(const ParameterList& parameters);

      void setParameters(const ParameterList& parameters)   throw (ParameterNotFoundException, ConstraintException);

      /**
       * @brief test if data fits this model
       *
       */
      
      virtual bool isCompatibleWith(const SiteContainer& data) const;

    };
} // end of namespace bpp.

#endif  // _MULTI_PROCESS_SEQUENCE_EVOLUTION_H_

