// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

// From bpp-core:
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Numeric/AbstractParametrizable.h>
#include <Bpp/Numeric/Matrix/Matrix.h>
#include <Bpp/Numeric/NumTools.h>

// From bpp-seq:
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/SequenceTools.h>

// From bpp-phyl
#include <Bpp/Phyl/Model/TwoParameterBinarySubstitutionModel.h>
#include <Bpp/Phyl/Tree/TreeTemplate.h>
#include <Bpp/Phyl/Tree/TreeTemplateTools.h>

// from std
#include <string>
#include <vector>
#include <iostream>

using namespace bpp;
using namespace std;


int main()
{
  try
  {
    // create a binary model
    auto alphabet = make_shared<BinaryAlphabet>();
    double mu = 42.;
    double pi0 = 0.1;
    auto model = make_unique<TwoParameterBinarySubstitutionModel>(alphabet, mu, pi0); // second arguent stands for mu

    // make sure the parameters of the model were assihned correctly
    double muVal = model->getParameterValue("mu");
    if (muVal != mu)
    {
      cerr << "Setting of initial mu parameter failed. Value is " << muVal << " instead of " << mu << endl;
      return 1;
    }
    double pi0Val = model->getParameterValue("pi0");
    if (pi0Val != pi0)
    {
      cerr << "Setting of initial mu parameter failed. Value is " << pi0Val << " instead of " << pi0 << endl;
      return 1;
    }

    // make sure the generator matrix is set correctly
    RowMatrix<double> refGenerator;
    refGenerator.resize(alphabet->getSize(), alphabet->getSize());
    refGenerator(0, 0) = -1 * muVal * (1 - pi0Val);
    refGenerator(0, 1) = muVal * (1 - pi0Val);
    refGenerator(1, 0) = muVal * pi0Val;
    refGenerator(1, 1) = -1 * muVal * pi0Val;
    for (size_t i = 0; i < alphabet->getSize(); i++)
    {
      for (size_t j = 0; j < alphabet->getSize(); j++)
      {
        double rate = model->Qij(i, j);
        if (rate != refGenerator(i, j))
        {
          cerr << "Q matrix definition is faulty. Entry at (" << i << ", " << j << ") is: " << rate << " instead of " << refGenerator(i, j) << endl;
          return 1;
        }
      }
    }

    // make sure the transition matrix is set correctly
    RowMatrix<double> refTransition;
    refTransition.resize(alphabet->getSize(), alphabet->getSize());
    double branchLen = 1;

    // case 1: branch of length 1
    double expVal = exp(-1 * muVal * branchLen);
    refTransition(0, 0) = (1 - pi0Val) + pi0Val * expVal;
    refTransition(0, 1) = pi0Val * (1 - expVal);
    refTransition(1, 0) =  (1 - pi0Val) * (1 - expVal);
    refTransition(1, 1) = pi0Val + (1 - pi0Val) * expVal;
    for (size_t i = 0; i < alphabet->getSize(); i++)
    {
      for (size_t j = 0; j < alphabet->getSize(); j++)
      {
        double transProb = model->Pij_t(i, j, branchLen);
        if (transProb != refTransition(i, j))
        {
          cerr << "Transition matrix definition is faulty for branch length " << branchLen << ". Entry (" << i << ", " << j << ") is: " << transProb << " instead of " << refTransition(i, j) << endl;
          return 1;
        }
      }
    }

    // case 2: branch of length 2
    branchLen = 2;
    expVal = exp(-1 * muVal * branchLen);
    refTransition(0, 0) = (1 - pi0Val) + pi0Val * expVal;
    refTransition(0, 1) = pi0Val * (1 - expVal);
    refTransition(1, 0) =  (1 - pi0Val) * (1 - expVal);
    refTransition(1, 1) = pi0Val + (1 - pi0Val) * expVal;
    for (size_t i = 0; i < alphabet->getSize(); i++)
    {
      for (size_t j = 0; j < alphabet->getSize(); j++)
      {
        double transProb = model->Pij_t(i, j, branchLen);
        if (transProb != refTransition(i, j))
        {
          cerr << "Transition matrix definition is faulty for branch length " << branchLen << ". Entry (" << i << ", " << j << ") is: " << transProb << " instead of " << refTransition(i, j) << endl;
          return 1;
        }
      }
    }
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }
  return 0;
}
