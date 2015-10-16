//
// File: HmmSequenceEvolution.h
// Created by: Laurent Guéguen
// Created on: vendredi 1 mai 2015, à 21h 05
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

#ifndef _HMM_SEQUENCE_EVOLUTION_H_
#define _HMM_SEQUENCE_EVOLUTION_H_


#include "MultiProcessSequenceEvolution.h"
#include "HmmProcessAlphabet.h"

// From Numeric
#include <Bpp/Numeric/Hmm/FullHmmTransitionMatrix.h>

#include <memory>

namespace bpp
{
/**
 * @brief Sequence evolution framework based on a hmm
 *
 */

    class HmmSequenceEvolution :
      public MultiProcessSequenceEvolution
    {
    private:
      std::auto_ptr<HmmProcessAlphabet> hmmAlph_;
      std::auto_ptr<FullHmmTransitionMatrix> hmmTransMat_;
  
    public:
      HmmSequenceEvolution(
        SubstitutionProcessCollection* processColl,
        std::vector<size_t>& nProc);

      HmmSequenceEvolution(const HmmSequenceEvolution& mlc) :
        MultiProcessSequenceEvolution(mlc),
        hmmAlph_(mlc.hmmAlph_->clone()),
        hmmTransMat_(mlc.hmmTransMat_->clone()) {}

      HmmSequenceEvolution& operator=(const HmmSequenceEvolution& mlc)
      {
        MultiProcessSequenceEvolution::operator=(mlc);
        hmmAlph_=std::auto_ptr<HmmProcessAlphabet>(new HmmProcessAlphabet(*mlc.hmmAlph_.get()));
        hmmTransMat_=std::auto_ptr<FullHmmTransitionMatrix>(new FullHmmTransitionMatrix(*mlc.hmmTransMat_.get()));
    
        return *this;
      }

      virtual ~HmmSequenceEvolution() {}

      HmmSequenceEvolution* clone() const { return new HmmSequenceEvolution(*this); }

    public:

      void setNamespace(const std::string& nameSpace);

      void fireParameterChanged(const ParameterList& parameters);

      const FullHmmTransitionMatrix& getHmmTransitionMatrix() const
      {
        return *hmmTransMat_.get();
      }

      FullHmmTransitionMatrix& getHmmTransitionMatrix()
      {
        return *hmmTransMat_.get();
      }

      const HmmProcessAlphabet& getHmmProcessAlphabet() const
      {
        return *hmmAlph_.get();
      }

      HmmProcessAlphabet& getHmmProcessAlphabet()
      {
        return *hmmAlph_.get();
      }

    };
} // end of namespace bpp.

#endif  // _HMMPHYLOLIKELIHOOD_H_

