// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

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
      return CWiseAdd<DataLik, std::tuple<DataLik, DataLik>>::create(context_, {left, right}, Dimension<DataLik> ());
    case '-':
      return CWiseSub<DataLik, std::tuple<DataLik, DataLik>>::create(context_, {left, right}, Dimension<DataLik> ());
    case '/':
    {
      auto inv = CWiseInverse<DataLik>::create(context_, {right}, Dimension<DataLik> ());
      return CWiseMul<DataLik, std::tuple<DataLik, DataLik>>::create(context_, {left, inv}, Dimension<DataLik> ());
    }
    case '*':
      return CWiseMul<DataLik, std::tuple<DataLik, DataLik>>::create(context_, {left, right}, Dimension<DataLik> ());
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


  auto func = dynamic_pointer_cast<FunctionOperator<SecondOrderDerivable>>(op);
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
