// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "SubstitutionRegister.h"

using namespace bpp;
using namespace std;

void GeneralSubstitutionRegister::updateTypes_()
{
  types_.clear();
  for (size_t i = 0; i < size_; i++)
  {
    for (size_t j = 0; j < size_; j++)
    {
      size_t type = matrix_(i, j);
      map<size_t, vector<size_t>> reg = types_[type];
      reg[i].push_back(j);
    }
  }
}
