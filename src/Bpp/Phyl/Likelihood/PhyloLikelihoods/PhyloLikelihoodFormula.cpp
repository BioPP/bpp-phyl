//
// File: PhyloLikelihoodFormula.cpp
// Authors:
//   Laurent GuÃÂ©guen
// Created: jeudi 8 dÃÂ©cembre 2016, ÃÂ  15h 23
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

#include <Bpp/Numeric/Function/Operators/BinaryOperator.h>
#include <Bpp/Numeric/Function/Operators/ConstantOperator.h>
#include <Bpp/Numeric/Function/Operators/FunctionOperator.h>
#include <Bpp/Numeric/Function/Operators/MathOperator.h>
#include <Bpp/Numeric/Function/Operators/NegativeOperator.h>

#include "PhyloLikelihoodFormula.h"

using namespace bpp;
using namespace std;

PhyloLikelihoodFormula::PhyloLikelihoodFormula(
    Context& context,
    std::shared_ptr<PhyloLikelihoodContainer> pC,
    const std::string& formula,
    bool inCollection) :
  AbstractPhyloLikelihood(context),
  AbstractParametrizable(""),
  AbstractPhyloLikelihoodSet(context, pC, {}, inCollection),
  compTree_(),
  likCal_(make_shared<LikelihoodCalculation>(context))
{
  readFormula(formula, inCollection);
  likCal_->setLikelihoodNode(makeLikelihoods());

  // using bpp::DotOptions;
  // writeGraphToDot(
  //   "formula.dot", {likCal_->getLikelihoodNode().get()});//, DotOptions::DetailedNodeInfo | DotOp
}


void PhyloLikelihoodFormula::readFormula(const std::string& formula, bool inCollection)
{
  map<string, shared_ptr<FunctionInterface>> functionNames;

  const vector<size_t>& nPhyl = getPhyloContainer()->getNumbersOfPhyloLikelihoods();

  for (size_t i = 0; i < nPhyl.size(); i++)
  {
    functionNames["phylo" + TextTools::toString(nPhyl[i])] = getPhyloLikelihood(nPhyl[i]);
  }

  compTree_ = make_unique<ComputationTree>(formula, functionNames);

  // add used Phylolikelihoods

  string popout = output();

  StringTokenizer st(popout, "phylo", true, true);
  st.nextToken();

  vector<size_t> phyldep;

  while (st.hasMoreToken())
  {
    string ex = st.nextToken();
    size_t np = (size_t)(atoi(ex.c_str()));
    if (!hasPhyloLikelihood(np))
      addPhyloLikelihood(np, inCollection ? "" : "_" + ex);
  }
}

std::string PhyloLikelihoodFormula::output() const
{
  return compTree_->output();
}

ValueRef<DataLik> PhyloLikelihoodFormula::makeLikelihoodsFromOperator(std::shared_ptr<Operator> op)
{
  auto cst = dynamic_pointer_cast<ConstantOperator>(op);
  if (cst)
    return NumericConstant<DataLik>::create(context_, cst->getValue());

  auto neg = dynamic_pointer_cast<NegativeOperator>(op);
  if (neg)
  {
    auto sonDF = makeLikelihoodsFromOperator(neg->getSon());
    return CWiseNegate<DataLik>::create(context_, {sonDF}, Dimension<DataLik> ());
  }

  auto bin = dynamic_pointer_cast<BinaryOperator>(op);
  if (bin)
  {
    auto left = makeLikelihoodsFromOperator(bin->getLeftSon());
    auto right = makeLikelihoodsFromOperator(bin->getRightSon());

    switch (bin->getSymbol())
    {
    case '+':
      return CWiseAdd<DataLik, std::tuple<DataLik, DataLik> >::create(context_, {left, right}, Dimension<DataLik> ());
    case '-':
      return CWiseSub<DataLik, std::tuple<DataLik, DataLik> >::create(context_, {left, right}, Dimension<DataLik> ());
    case '/':
    {
      auto inv = CWiseInverse<DataLik>::create(context_, {right}, Dimension<DataLik> ());
      return CWiseMul<DataLik, std::tuple<DataLik, DataLik> >::create(context_, {left, inv}, Dimension<DataLik> ());
    }
    case '*':
      return CWiseMul<DataLik, std::tuple<DataLik, DataLik> >::create(context_, {left, right}, Dimension<DataLik> ());
    default:
      return NumericConstant<DataLik>::create(context_, 0);
    }
  }

  auto mat = dynamic_pointer_cast<MathOperator>(op);
  if (mat)
  {
    auto sonDF = makeLikelihoodsFromOperator(mat->getSon());
    auto name = mat->getName();
    if (name == "exp")
      return CWiseExp<DataLik>::create(context_, {sonDF}, Dimension<DataLik> ());
    else if (name == "log")
      return CWiseLog<DataLik>::create(context_, {sonDF}, Dimension<DataLik> ());
    else
      throw Exception("PhyloLikelihoodFormula::Makelikelihoodsfromoperator : unknown function " + name + ". Ask developpers.");
  }


  auto func = dynamic_pointer_cast<FunctionOperator<SecondOrderDerivable> >(op);
  if (func)
  {
    auto name = func->getName();
    if (name.substr(0, 5) == "phylo")
    {
      auto phyl = getPhyloLikelihood(TextTools::to<size_t>(name.substr(5, string::npos)));
      shareParameters_(phyl->getParameters());

      return phyl->getLikelihoodNode();
    }
  }

  throw Exception("PhyloLikelihoodFormula::makeLikelihoodsFromOperator : Unknown operator: " + op->output());
}
