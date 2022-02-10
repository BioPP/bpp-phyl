//
// File: AbstractSubstitutionProcess.cpp
// Authors:
//   Julien Dutheil
// Created: Tue Marc 22 21:17 2013
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/


#include "AbstractAutonomousSubstitutionProcess.h"

using namespace bpp;
using namespace std;

AbstractAutonomousSubstitutionProcess::AbstractAutonomousSubstitutionProcess(std::shared_ptr<const PhyloTree> tree, std::shared_ptr<FrequencySet> rootFrequencies, const string& prefix) :
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

AbstractAutonomousSubstitutionProcess::AbstractAutonomousSubstitutionProcess(std::shared_ptr<ParametrizablePhyloTree> tree, std::shared_ptr<FrequencySet> rootFrequencies, const string& prefix) :
  AbstractParameterAliasable(prefix),
  pTree_(tree),
  rootFrequencies_(rootFrequencies),
  modelScenario_(0)
{
  if (tree)
    addParameters_(tree->getParameters());  // Branch lengths
  if (rootFrequencies_)
    addParameters_(rootFrequencies_->getParameters());
}

AbstractAutonomousSubstitutionProcess::AbstractAutonomousSubstitutionProcess(const AbstractAutonomousSubstitutionProcess& asp) :
  AbstractParameterAliasable(asp),
  pTree_(asp.pTree_?asp.pTree_->clone():0),
  rootFrequencies_(asp.rootFrequencies_?asp.rootFrequencies_->clone():0),
  modelScenario_(asp.modelScenario_) // this has to be specified by inheriting class to follow model links
{}

AbstractAutonomousSubstitutionProcess& AbstractAutonomousSubstitutionProcess::operator=(const AbstractAutonomousSubstitutionProcess& asp)
{
  AbstractParameterAliasable::operator=(*this);

  pTree_.reset(asp.pTree_?asp.pTree_->clone():0);
  rootFrequencies_.reset(asp.rootFrequencies_?asp.rootFrequencies_->clone():0);
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
  
