//
// File: PartitionSequenceEvolution.h
// Created by: Laurent Guéguen
// Created on: vendredi 15 mai 2015, à 18h 29
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

#ifndef _PARTITION_SEQUENCE_EVOLUTION_H_
#define _PARTITION_SEQUENCE_EVOLUTION_H_


#include "MultiProcessSequenceEvolution.h"

namespace bpp
{
/**
 * @brief Sequence evolution framework based on a mixture of
 * substitution processes
 *
 * @see MultiProcessSequenceEvolution
 */

    class PartitionSequenceEvolution :
      public MultiProcessSequenceEvolution
    {
    private:

      /**
       *@ brief vector of the substitution process numbers along the sequence.
       *
       */
      
      std::vector<size_t> vProc_;

      size_t vSize_;
      
      /**
       * @brief On the reverse, for each process number, the vector
       * of the sites where it is used.
       *
       * Convenient for process specific site patterns.
       */

      std::map<size_t, std::vector<size_t> > mProcPos_;
    
    public:

      /*
       * @brief constructor
       *
       * @param the used SubstitutionProcessCollection
       * @param A vector of the number of the processes along the sequence.
       *
       */
      
      PartitionSequenceEvolution(
        SubstitutionProcessCollection* processColl,
        std::vector<size_t>& posProc);
      
      PartitionSequenceEvolution(const PartitionSequenceEvolution& mlc) :
        MultiProcessSequenceEvolution(mlc),
        vProc_(mlc.vProc_),
        vSize_(mlc.vSize_),
        mProcPos_(mlc.mProcPos_) {}

      PartitionSequenceEvolution& operator=(const PartitionSequenceEvolution& mlc)
      {
        MultiProcessSequenceEvolution::operator=(mlc);
        vProc_ = mlc.vProc_;
        vSize_ = mlc.vSize_;
        
        mProcPos_ = mlc.mProcPos_;
        
        return *this;
      }

      virtual ~PartitionSequenceEvolution() {}

      PartitionSequenceEvolution* clone() const { return new PartitionSequenceEvolution(*this); }

    public:
      size_t getNumberOfSites() const
      {
        return vSize_;
      }
      
      size_t getSubstitutionProcessNumber(size_t i) const
      {
        if (i>=vSize_)
          throw IndexOutOfBoundsException("PartitionSequenceEvolution::getSubstitutionProcess",i,0,vSize_);
        return vProc_[i];
      }

      std::map<size_t, std::vector<size_t> >& getMapOfProcessSites()
      {
        return mProcPos_;
      }

      const std::map<size_t, std::vector<size_t> >& getMapOfProcessSites() const
      {
        return mProcPos_;
      }

    };
} // end of namespace bpp.

#endif  // _PARTITION_SEQUENCE_EVOLUTION_H_

