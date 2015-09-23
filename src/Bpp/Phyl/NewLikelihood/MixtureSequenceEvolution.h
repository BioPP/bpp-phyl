//
// File: MixtureSequenceEvolution.h
// Created by: Laurent Guéguen
// Created on: jeudi 30 avril 2015, à 17h 16
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

#ifndef _MIXTURE_SEQUENCE_EVOLUTION_H_
#define _MIXTURE_SEQUENCE_EVOLUTION_H_


#include "MultiProcessSequenceEvolution.h"

#include <Bpp/Numeric/Prob/Simplex.h>

namespace bpp
{
/**
 * @brief Sequence evolution framework based on a mixture of
 * substitution processes
 *
 * @see MultiProcessSequencePhyloLikelihood
 */

    class MixtureSequenceEvolution :
      public MultiProcessSequenceEvolution
    {
    private:
      Simplex simplex_;

    public:
      MixtureSequenceEvolution(
        SubstitutionProcessCollection* processColl,
        std::vector<size_t>& nProc);
      
      MixtureSequenceEvolution(const MixtureSequenceEvolution& mlc) :
        MultiProcessSequenceEvolution(mlc),
        simplex_(mlc.simplex_) {}

      MixtureSequenceEvolution& operator=(const MixtureSequenceEvolution& mlc)
      {
        MultiProcessSequenceEvolution::operator=(mlc);
        simplex_ = mlc.simplex_;
        return *this;
      }

      virtual ~MixtureSequenceEvolution() {}

      MixtureSequenceEvolution* clone() const { return new MixtureSequenceEvolution(*this); }

    public:
      void setNamespace(const std::string& nameSpace);

      void fireParameterChanged(const ParameterList& parameters);

      ParameterList getNonDerivableParameters() const;

      const std::vector<double>& getSubProcessProbabilities() const
      {
        return simplex_.getFrequencies();
      }
      
        
      /**
       * @brief return the probability of a  subprocess
       *
       * @param i the index of the subprocess
       */
  
      double getSubProcessProb(size_t i) const
      {
        return simplex_.prob(i);
      }

      /**
       * @brief Set the probabilities of the subprocess
       *
       */
  
      void setSubProcessProb(const Simplex& si);
  
    };
} // end of namespace bpp.

#endif  // _MIXTURE_SEQUENCE_EVOLUTION_H_

