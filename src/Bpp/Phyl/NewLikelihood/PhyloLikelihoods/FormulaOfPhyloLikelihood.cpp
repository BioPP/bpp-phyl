//
// File: FormulaOfPhyloLikelihood.cpp
// Created by: Laurent Guéguen
// Created on: jeudi 8 décembre 2016, à 15h 23
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

#include "FormulaOfPhyloLikelihood.h"

#include <Bpp/Numeric/Function/Operators/BinaryOperator.h>
#include <Bpp/Numeric/Function/Operators/FunctionOperator.h>
#include <Bpp/Numeric/Function/Operators/ConstantOperator.h>
#include <Bpp/Numeric/Function/Operators/NegativeOperator.h>
#include <Bpp/Numeric/Function/Operators/MathOperator.h>

using namespace bpp;
using namespace std;

FormulaOfPhyloLikelihood::FormulaOfPhyloLikelihood(Context& context, PhyloLikelihoodContainer* pC, const std::string& formula, bool inCollection) :
  AbstractPhyloLikelihood(context),
  SetOfAbstractPhyloLikelihood(context, pC, {}, inCollection),
  compTree_(),
  likCal_(new AlignedLikelihoodCalculation(context))
{
  readFormula(formula, inCollection);
  likCal_->setLikelihoodNode(makeLikelihoods());
}


FormulaOfPhyloLikelihood::FormulaOfPhyloLikelihood(const FormulaOfPhyloLikelihood& sd) :
  AbstractPhyloLikelihood(sd),
  SetOfAbstractPhyloLikelihood(sd),
  compTree_(sd.compTree_->clone()),
  likCal_(sd.likCal_)  
{
}

void FormulaOfPhyloLikelihood::readFormula(const std::string& formula, bool inCollection)
{
  std::map<std::string, Function*> functionNames;

  const vector<size_t>& nPhyl = getPhyloContainer()->getNumbersOfPhyloLikelihoods();
  
  for (size_t i = 0; i < nPhyl.size(); i++)
  {
    functionNames["phylo"+TextTools::toString(nPhyl[i])]=getPhyloLikelihood(nPhyl[i]);
  }
  
  compTree_=unique_ptr<ComputationTree>(new ComputationTree(formula, functionNames));

  // add used Phylolikelihoods
  
  string popout=output();
  
  StringTokenizer st(popout,"phylo",true, true);
  st.nextToken();

  vector<size_t> phyldep;
  
  while (st.hasMoreToken())
  {
    string ex=st.nextToken();
    size_t np =(size_t)(atoi(ex.c_str()));
    if (!hasPhyloLikelihood(np))
      addPhyloLikelihood(np, inCollection?"":"_"+ex);
  }
}

std::string FormulaOfPhyloLikelihood::output() const
{
  return compTree_->output();
}

ValueRef<double> FormulaOfPhyloLikelihood::makeLikelihoodsFromOperator(std::shared_ptr<Operator> op)
{
  auto cst = dynamic_pointer_cast<ConstantOperator>(op);
  if (cst)
    return NumericConstant<double>::create(context_, cst->getValue());

  auto neg = dynamic_pointer_cast<NegativeOperator>(op);
  if (neg)
  {
    auto sonDF = makeLikelihoodsFromOperator(neg->getSon());
    return CWiseNegate<double>::create(context_, {sonDF}, Dimension<double> ());
  }

  auto bin = dynamic_pointer_cast<BinaryOperator>(op);
  if (bin)
  {
    auto left = makeLikelihoodsFromOperator(bin->getLeftSon());
    auto right = makeLikelihoodsFromOperator(bin->getRightSon());
    
    switch(bin->getSymbol())
    {
    case '+':
      return CWiseAdd<double, std::tuple<double, double>>::create(context_, {left, right}, Dimension<double> ());
    case '-':
      return CWiseSub<double, std::tuple<double, double>>::create(context_, {left, right}, Dimension<double> ());
    case '/':
      {
        auto inv = CWiseInverse<double>::create(context_, {right}, Dimension<double> ());
        return CWiseMul<double, std::tuple<double, double>>::create(context_, {left, inv}, Dimension<double> ());
      }
    case '*':
      return CWiseMul<double, std::tuple<double, double>>::create(context_, {left, right}, Dimension<double> ());
    default:
      return NumericConstant<double>::create(context_, 0);
    }
  }
  
  auto mat = dynamic_pointer_cast<MathOperator>(op);
  if (mat)
  {
    auto sonDF = makeLikelihoodsFromOperator(mat->getSon());
    auto name = mat->getName();
    if (name=="exp")
      return CWiseExp<double>::create(context_, {sonDF}, Dimension<double> ());
    else if (name=="log")
      return CWiseLog<double>::create(context_, {sonDF}, Dimension<double> ());
    else
      throw Exception("FormulaOfPhyloLikelihood::Makelikelihoodsfromoperator : unknown function " + name + ". Ask developpers.");
  }


  auto func = dynamic_pointer_cast<FunctionOperator<DerivableSecondOrder>>(op);
  if (func)
  {
    auto name = func->getName();
    if (name.substr(0,5)=="phylo")
    {
      auto phyl  = getAbstractPhyloLikelihood((size_t)(TextTools::toInt(name.substr(5, string::npos))));
      return phyl->getLikelihoodNode();
    }
  }

  throw Exception("FormulaOfPhyloLikelihood::makeLikelihoodsFromOperator : Unknown operator: " + op->output());
}





