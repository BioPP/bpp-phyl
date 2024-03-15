// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "PhyloBranch.h"
#include "PhyloBranchParam.h"

using namespace std;
using namespace bpp;

PhyloBranchParam::PhyloBranchParam(const PhyloBranch& branch) :
  AbstractParametrizable("")
{
  double brLen = NumConstants::SMALL();
  if (branch.hasLength() && branch.getLength() >= NumConstants::SMALL())
    brLen = branch.getLength();
  addParameter_(new Parameter("BrLen", brLen, Parameter::R_PLUS_STAR));
}
