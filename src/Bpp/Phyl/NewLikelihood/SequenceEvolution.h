//
// File: SequenceEvolution.h
// Created by: Laurent Guéguen
// Created on: mardi 28 avril 2015, à 10h 51
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

#ifndef _SEQUENCE_EVOLUTION_H_
#define _SEQUENCE_EVOLUTION_H_


#include "SubstitutionProcess.h"

//From bpp-core:
#include <Bpp/Numeric/ParameterAliasable.h>
#include <Bpp/Numeric/AbstractParameterAliasable.h>

//From bpp-seq:

#include <Bpp/Seq/Container/SiteContainer.h>

//From the STL:
#include <memory>

namespace bpp
{

/**
 * @brief This interface describes the evolution process of a sequence.
 *
 * It main purpose is to provide the necessary calculus for each
 * management of site specific substitution process, such as partition
 * models and HMM.
 *
 * This object has the INDEPENDENT parameters of the processes. 
 * 
 */
  
  class SequenceEvolution :
    public virtual ParameterAliasable
  {
  public:
    virtual SequenceEvolution* clone() const = 0;

  public:
    virtual bool isCompatibleWith(const SiteContainer& data) const = 0;

    virtual bool hasDerivableParameter(const std::string& name) const = 0;

    virtual const std::vector<size_t>& getSubstitutionProcessNumbers() const = 0;
    
    virtual const SubstitutionProcess& getSubstitutionProcess(size_t number) const = 0;
    
    /**
     * @brief Get the branch lengths parameters.
     *
     * @return A ParameterList with all branch lengths.
     */

    virtual ParameterList getBranchLengthParameters(bool independent) const = 0;
    
    /**
     * @brief Get the parameters associated to substitution model(s).
     *
     * @return A ParameterList.
     */

    virtual ParameterList getSubstitutionModelParameters(bool independent) const = 0;

    /**
     * @brief Get the parameters associated to substitution processes(s).
     *
     * @return A ParameterList.
     */

    virtual ParameterList getSubstitutionProcessParameters(bool independent) const = 0;

    /**
     * @brief Get the parameters associated to the rate distribution(s).
     *
     * @return A ParameterList.
     */

    virtual ParameterList getRateDistributionParameters(bool independent) const = 0;

    /**
     * @brief Get the parameters associated to the root frequencies(s).
     *
     * @return A ParameterList.
     */

    virtual ParameterList getRootFrequenciesParameters(bool independent) const = 0;


    /**
     * @brief All independent derivable parameters.
     *
     * Usually, this contains all branch lengths parameters.
     *
     * @return A ParameterList.
     */

    virtual ParameterList getDerivableParameters() const = 0;

    /**
     * @brief All independent non derivable parameters.
     *
     * Usually, this contains all substitution model parameters and rate distribution.
     *
     * @return A ParameterList.
     */

    virtual ParameterList getNonDerivableParameters() const = 0;

  };

  // class AbstractSequenceEvolution :
  //   public virtual SequenceEvolution,
  //   public virtual AbstractParameterAliasable
  // {
  // protected:
  //   size_t seqLength_;

  // public:
  //   AbstractSequenceEvolution() :
  //     AbstractParameterAliasable(""),
  //     seqLength_(0)
  //   {
  //   }

  //   AbstractSequenceEvolution(const AbstractSequenceEvolution& evol) :
  //     AbstractParameterAliasable(evol),
  //     seqLength_(evol.seqLength_)
  //   {
  //   }

  //   AbstractSequenceEvolution& operator=(const AbstractSequenceEvolution& evol)
  //   {
  //     AbstractParameterAliasable::operator=(*this);
  //     seqLength_=evol.seqLength_;

  //     return *this;
  //   }
    
  //   size_t getNumberOfPositions() const
  //   {
  //     return seqLength_;
  //   }

  //   void setNumberOfPositions(size_t length) 
  //   {
  //     seqLength_=length;
  //   }

//  };

} // end namespace bpp

#endif // _SEQUENCE_EVOLUTION_H_

