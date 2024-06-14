// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "GlobalClockTreeLikelihoodFunctionWrapper.h"

using namespace bpp;

void GlobalClockTreeLikelihoodFunctionWrapper::fireParameterChanged(const bpp::ParameterList& pl)
{
  // filter parameters:
  ParameterList pl2;
  bool recomputeHeights = false;
  for (unsigned int i = 0; i < pl.size(); ++i)
  {
    if (pl[i].getName().substr(0, 7) == "HeightP" || pl[i].getName() == "TotalHeight")
      recomputeHeights = true;
    else
      pl2.addParameter(pl[i]);
  }
  if (recomputeHeights)
  {
    TreeTemplate<Node> tree(tl_->tree());
    computeBranchLengthsFromHeights_(tree.getRootNode(), parameter("TotalHeight").getValue(), pl2);
  }
  tl_->setParameters(pl2);
}

ParameterList GlobalClockTreeLikelihoodFunctionWrapper::getHeightParameters() const
{
  ParameterList pl;

  for (unsigned int i = 0; i < getNumberOfParameters(); ++i)
  {
    Parameter p = getParameter_(i);
    if (p.getName().substr(0, 7) == "HeightP" || p.getName() == "TotalHeight")
      pl.addParameter(p);
  }
  return pl;
}

void GlobalClockTreeLikelihoodFunctionWrapper::initParameters_()
{
  // Check if the tree is rooted:
  TreeTemplate<Node> tree(tl_->tree());
  if (!tree.isRooted())
    throw Exception("GlobalClockTreeLikelihoodFunctionWrapper::initParameters_(). Tree is unrooted!");
  if (TreeTemplateTools::isMultifurcating(*(tree.getRootNode())))
    throw Exception("GlobalClockTreeLikelihoodFunctionWrapper::initParameters_(). Tree is multifurcating.");
  std::map<const Node*, double> heights;
  TreeTemplateTools::getHeights(*(tree.getRootNode()), heights);
  double totalHeight = heights[tree.getRootNode()];
  addParameter_(new Parameter("TotalHeight", totalHeight, Parameter::R_PLUS_STAR));
  for (std::map<const Node*, double>::iterator it = heights.begin(); it != heights.end(); it++)
  {
    if (!it->first->isLeaf() && it->first->hasFather())
    {
      double fatherHeight = heights[it->first->getFather()];
      addParameter_(new Parameter("HeightP" + TextTools::toString(it->first->getId()), it->second / fatherHeight, Parameter::PROP_CONSTRAINT_IN));
    }
  }
  // We add other parameters:
  ParameterList pl = tl_->getParameters();
  for (unsigned int i = 0; i < pl.size(); ++i)
  {
    if (pl[i].getName().substr(0, 5) != "BrLen")
      addParameter_(pl[i].clone());
  }
  // Compute everything:
  fireParameterChanged(getParameters());
}

void GlobalClockTreeLikelihoodFunctionWrapper::computeBranchLengthsFromHeights_(const Node* node, double height, ParameterList& brlenPl)
{
  for (unsigned int i = 0; i < node->getNumberOfSons(); i++)
  {
    const Node* son = node->getSon(i);
    if (son->isLeaf())
    {
      brlenPl.addParameter(Parameter("BrLen" + TextTools::toString(son->getId()), std::max(0.0000011, height), std::make_shared<IntervalConstraint>(1, 0.000001, false)));
    }
    else
    {
      double sonHeightP = parameter("HeightP" + TextTools::toString(son->getId())).getValue();
      double sonHeight = sonHeightP * height;
      brlenPl.addParameter(Parameter("BrLen" + TextTools::toString(son->getId()), std::max(0.0000011, height - sonHeight), std::make_shared<IntervalConstraint>(1, 0.000001, false)));
      computeBranchLengthsFromHeights_(son, sonHeight, brlenPl);
    }
  }
}
