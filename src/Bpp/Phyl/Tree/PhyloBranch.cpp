// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "PhyloBranch.h"
#include "PhyloBranchParam.h"

using namespace bpp;
using namespace std;

/** Copy constructor: *********************************************************/

PhyloBranch::PhyloBranch(const PhyloBranch& branch) :
  isLengthDefined_(branch.isLengthDefined_),
  length_(branch.length_),
  properties_()
{
  for (map<string, Clonable*>::iterator i = branch.properties_.begin(); i != branch.properties_.end(); i++)
  {
    properties_[i->first] = i->second->clone();
  }
}


PhyloBranch::PhyloBranch(const PhyloBranchParam& branch) :
  isLengthDefined_(true),
  length_(branch.getLength()),
  properties_()
{}

/** Assignation operator: *****************************************************/

PhyloBranch& PhyloBranch::operator=(const PhyloBranch& branch)
{
  isLengthDefined_ = branch.isLengthDefined_;
  length_ = branch.length_;

  for (map<string, Clonable*>::iterator i = branch.properties_.begin(); i != branch.properties_.end(); i++)
  {
    Clonable* p = properties_[i->first];
    if (p)
      delete p;
    properties_[i->first] = i->second->clone();
  }
  return *this;
}
