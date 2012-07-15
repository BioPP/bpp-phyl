//
// File: GlobalClockTreeLikelihoodFunctionWrapper.cpp
// Created by: Julien Dutheil
// Created on: Thu Jul 14 10:53 2011
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
    TreeTemplate<Node> tree(tl_->getTree());
    computeBranchLengthsFromHeights_(tree.getRootNode(), getParameter("TotalHeight").getValue(), pl2);
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
  TreeTemplate<Node> tree(tl_->getTree());
  if (!tree.isRooted()) throw Exception("GlobalClockTreeLikelihoodFunctionWrapper::initParameters_(). Tree is unrooted!");
  if (TreeTemplateTools::isMultifurcating(*(tree.getRootNode()))) throw Exception("GlobalClockTreeLikelihoodFunctionWrapper::initParameters_(). Tree is multifurcating.");
  std::map<const Node*, double> heights;
  TreeTemplateTools::getHeights(*(tree.getRootNode()), heights);
  double totalHeight = heights[tree.getRootNode()];
  addParameter_(new Parameter("TotalHeight", totalHeight, &Parameter::R_PLUS_STAR));
  for (std::map<const Node*, double>::iterator it = heights.begin(); it != heights.end(); it++)
  {
    if (!it->first->isLeaf() && it->first->hasFather())
    {
      double fatherHeight = heights[it->first->getFather()];
      addParameter_(new Parameter("HeightP" + TextTools::toString(it->first->getId()), it->second / fatherHeight, &Parameter::PROP_CONSTRAINT_IN));
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

void GlobalClockTreeLikelihoodFunctionWrapper::computeBranchLengthsFromHeights_(const Node* node, double height, ParameterList& brlenPl) throw (Exception)
{
  for (unsigned int i = 0; i < node->getNumberOfSons(); i++)
  {
    const Node* son = node->getSon(i);
    if (son->isLeaf())
    {
      brlenPl.addParameter(Parameter("BrLen" + TextTools::toString(son->getId()), std::max(0.0000011, height), new IntervalConstraint(1, 0.000001, false), true));
    }
    else
    {
      double sonHeightP = getParameter("HeightP" + TextTools::toString(son->getId())).getValue();
      double sonHeight = sonHeightP * height;
      brlenPl.addParameter(Parameter("BrLen" + TextTools::toString(son->getId()), std::max(0.0000011, height - sonHeight), new IntervalConstraint(1, 0.000001, false), true));
      computeBranchLengthsFromHeights_(son, sonHeight, brlenPl);
    }
  }
}

