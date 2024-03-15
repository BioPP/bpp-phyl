// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/VectorTools.h>

#include "PartitionSequenceEvolution.h"

using namespace std;
using namespace bpp;

/******************************************************************************/

PartitionSequenceEvolution::PartitionSequenceEvolution(
  shared_ptr<SubstitutionProcessCollection> processColl,
  std::vector<size_t>& posProc) :
  MultiProcessSequenceEvolution(processColl, vector<size_t>()),
  vProc_(),
  vSize_(0),
  mProcPos_()
{
  size_t i = 0;
  for (auto nproc : posProc)
  {
    if (!processColl_->hasSubstitutionProcessNumber(nproc))
      throw BadIntegerException("PartitionSequenceEvolution::PartitionSequenceEvolution : unknown process number ", int(nproc));
    else
    {
      vProc_.push_back(nproc);
      if (mProcPos_.find(nproc) == mProcPos_.end())
        mProcPos_[nproc] = vector<size_t>();
      mProcPos_[nproc].push_back(i++);
    }
  }

  vSize_ = posProc.size();

  for (auto& it : mProcPos_)
  {
    nProc_.push_back(it.first);
  }
}

