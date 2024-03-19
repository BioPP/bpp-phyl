// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "AbstractAutonomousSubstitutionProcess.h"

using namespace bpp;
using namespace std;

AbstractAutonomousSubstitutionProcess::AbstractAutonomousSubstitutionProcess(
    shared_ptr<const PhyloTree> tree,
    shared_ptr<FrequencySetInterface> rootFrequencies,
    const string& prefix) :
  AbstractParameterAliasable(prefix),
  pTree_(0),
  rootFrequencies_(rootFrequencies),
  modelScenario_(0)
{
  if (tree)
    setPhyloTree(*tree);
  if (rootFrequencies_)
    addParameters_(rootFrequencies_->getParameters());
}

AbstractAutonomousSubstitutionProcess::AbstractAutonomousSubstitutionProcess(
    shared_ptr<ParametrizablePhyloTree> tree,
    shared_ptr<FrequencySetInterface> rootFrequencies, const string& prefix) :
  AbstractParameterAliasable(prefix),
  pTree_(tree),
  rootFrequencies_(rootFrequencies),
  modelScenario_(0)
{
  if (tree)
    addParameters_(tree->getParameters()); // Branch lengths
  if (rootFrequencies_)
    addParameters_(rootFrequencies_->getParameters());
}

AbstractAutonomousSubstitutionProcess::AbstractAutonomousSubstitutionProcess(const AbstractAutonomousSubstitutionProcess& asp) :
  AbstractParameterAliasable(asp),
  pTree_(asp.pTree_ ? asp.pTree_->clone() : 0),
  rootFrequencies_(asp.rootFrequencies_ ? asp.rootFrequencies_->clone() : 0),
  modelScenario_(asp.modelScenario_) // this has to be specified by inheriting class to follow model links
{}

AbstractAutonomousSubstitutionProcess& AbstractAutonomousSubstitutionProcess::operator=(const AbstractAutonomousSubstitutionProcess& asp)
{
  AbstractParameterAliasable::operator=(*this);

  pTree_.reset(asp.pTree_ ? asp.pTree_->clone() : 0);
  rootFrequencies_.reset(asp.rootFrequencies_ ? asp.rootFrequencies_->clone() : 0);
  modelScenario_ = asp.modelScenario_; // this has to be specified by inheriting class to follow model links
  return *this;
}

void AbstractAutonomousSubstitutionProcess::fireParameterChanged(const ParameterList& pl)
{
  ParameterList gAP = getAliasedParameters(pl);
  gAP.addParameters(pl);

  if (pTree_)
    pTree_->matchParametersValues(gAP);
  if (rootFrequencies_)
    rootFrequencies_->matchParametersValues(gAP);
}

void AbstractAutonomousSubstitutionProcess::setPhyloTree(const PhyloTree& phyloTree)
{
  if (pTree_)
    getParameters_().deleteParameters(pTree_->getParameters().getParameterNames(), false);

  pTree_.reset(new ParametrizablePhyloTree(phyloTree));
  addParameters_(pTree_->getParameters());
}
