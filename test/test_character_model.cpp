//
// File: test_models.cpp
// Created by: Julien Dutheil
// Created on: Tue Jan 22 22:18 2013
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for numerical calculus. This file is part of the Bio++ project.

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
    refGenerator.resize(alphabet->getSize(),alphabet->getSize());
    refGenerator(0,0) = -1 * muVal * (1-pi0Val);
    refGenerator(0,1) = muVal * (1-pi0Val);
    refGenerator(1,0) = muVal * pi0Val;
    refGenerator(1,1) = -1* muVal * pi0Val;
    for (size_t i=0; i < alphabet->getSize(); i++)
    {
      for (size_t j=0; j < alphabet->getSize(); j++)
      {
        double rate = model->Qij(i,j);
        if (rate != refGenerator(i,j))
        {
          cerr << "Q matrix definition is faulty. Entry at (" << i << ", " << j << ") is: " << rate << " instead of " << refGenerator(i,j) << endl;
          return 1;
        }
      }
    }

    // make sure the transtition matrix is set correctly
    RowMatrix<double> refTransition;
    refTransition.resize(alphabet->getSize(),alphabet->getSize());
    double branchLen = 1;

    // case 1: branch of length 1
    double expVal = exp(-1 * muVal * branchLen);
    refTransition(0,0) = (1 - pi0Val) + pi0Val * expVal;
    refTransition(0,1) = pi0Val * (1 - expVal);
    refTransition(1,0) =  (1 - pi0Val) * (1 - expVal);
    refTransition(1,1) = pi0Val + (1 - pi0Val) * expVal;
    for (size_t i=0; i < alphabet->getSize(); i++)
    {
      for (size_t j=0; j < alphabet->getSize(); j++)
      {
        double transProb = model->Pij_t(i,j, branchLen);
        if (transProb != refTransition(i,j))
        {
          cerr << "Transition matrix definition is faulty for branch length " << branchLen << ". Entry (" << i << ", " << j << ") is: " << transProb << " instead of " << refTransition(i,j) << endl;
          return 1;
        }
      }
    }

    // case 2: branch of length 2
    branchLen = 2;
    expVal = exp(-1 * muVal * branchLen);
    refTransition(0,0) = (1 - pi0Val) + pi0Val * expVal;
    refTransition(0,1) = pi0Val * (1 - expVal);
    refTransition(1,0) =  (1 - pi0Val) * (1 - expVal);
    refTransition(1,1) = pi0Val + (1 - pi0Val) * expVal;
    for (size_t i=0; i < alphabet->getSize(); i++)
    {
      for (size_t j=0; j < alphabet->getSize(); j++)
      {
        double transProb = model->Pij_t(i,j, branchLen);
        if (transProb != refTransition(i,j))
        {
          cerr << "Transition matrix definition is faulty for branch length " << branchLen << ". Entry (" << i << ", " << j << ") is: " << transProb << " instead of " << refTransition(i,j) << endl;
          return 1;
        }
      }
    }
  }
  catch (exception & e)
  {
    cout << e.what() << endl;
    return 1;
  }
  return 0;
}
