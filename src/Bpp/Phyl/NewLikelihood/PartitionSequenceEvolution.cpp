//
// File: PartitionSequenceEvolution.cpp
// Created by: Laurent Guéguen
// Created on: vendredi 15 mai 2015, à 22h 40
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

#include "PartitionSequenceEvolution.h"

#include <Bpp/Numeric/VectorTools.h>

using namespace std;
using namespace bpp;

/******************************************************************************/

PartitionSequenceEvolution::PartitionSequenceEvolution(
  SubstitutionProcessCollection* processColl,
  std::vector<size_t>& posProc) :
  MultiProcessSequenceEvolution(processColl, vector<size_t>()),
  vProc_(),
  vSize_(0),
  mProcPos_()
{
  for (size_t i=0;i<posProc.size();i++){
    size_t nproc=posProc[i];
    
    if (!processColl_->hasSubstitutionProcessNumber(nproc))
      throw BadIntegerException("PartitionSequenceEvolution::PartitionSequenceEvolution : unknown process number ", int(nproc));
    else
    {
      vProc_.push_back(nproc);
      if (mProcPos_.find(nproc)==mProcPos_.end())
        mProcPos_[nproc]=vector<size_t>();
      mProcPos_[nproc].push_back(i);
    }
  }
  
  vSize_=posProc.size();
  
  for (std::map<size_t, std::vector<size_t> >::const_iterator it=mProcPos_.begin(); it!=mProcPos_.end(); it++)
    nProc_.push_back(it->first);
  
}



