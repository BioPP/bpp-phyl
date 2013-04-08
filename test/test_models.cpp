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

#include <Bpp/Phyl/Model/Nucleotide/GTR.h>
#include <Bpp/Phyl/Model/Codon/YN98.h>
#include <Bpp/Phyl/Model/FrequenciesSet/CodonFrequenciesSet.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Alphabet/StandardCodonAlphabet.h>
#include <Bpp/Seq/GeneticCode/StandardGeneticCode.h>
#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/Function/ReparametrizationFunctionWrapper.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Numeric/AbstractParametrizable.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <iostream>

using namespace bpp;
using namespace std;

class DummyFunction:
  public virtual Function,
  public AbstractParametrizable
{
  public:
    DummyFunction(const SubstitutionModel& model):
      AbstractParametrizable("")
    {
      addParameters_(model.getParameters());
    }

    DummyFunction* clone() const { return new DummyFunction(*this); }

    void setParameters(const ParameterList& pl) throw (bpp::ParameterNotFoundException
, bpp::ConstraintException, bpp::Exception) {
      matchParametersValues(pl);
    }

    double getValue() const throw (Exception) { return 0; }

    void fireParameterChanged(const bpp::ParameterList&) {}

};

bool testModel(SubstitutionModel& model) {
  ParameterList pl = model.getParameters();
  DummyFunction df(model);
  ReparametrizationFunctionWrapper redf(&df, pl, false);

  //Try to modify randomly each parameter and check that the new parameter apply correctly:
  for (unsigned int i = 0; i < 10; ++i) {
    //Get random parameters (unconstrained):
    ParameterList pl2 = redf.getParameters();
    for (size_t j = 0; j < pl.size(); ++j) {
      double value = (RandomTools::flipCoin() ? 1. : -1) * RandomTools::giveRandomNumberBetweenZeroAndEntry(1);
      pl2[j].setValue(value);
    }
    //Apply unconstrained parameters:
    redf.setParameters(pl2);
    //Retrieve transformed parameters:
    pl2 = df.getParameters();
    //pl2.printParameters(cout);
    //Now apply the new parameters and retrieve them again:
    model.matchParametersValues(pl2);
    ParameterList pl3 = model.getParameters();
    //Compare the two lists:
    for (size_t j = 0; j < pl.size(); ++j) {
      if (abs(pl2[j].getValue() - pl3[j].getValue()) > 0.0000001) {
        cerr << "ERROR for parameter " << pl2[j].getName() << ": " << pl2[j].getValue() << "<>" << pl3[j].getValue() << endl;
        return false;
      }
    }
  }

  return true;
}

int main() {
  //Nucleotide models:
  GTR gtr(&AlphabetTools::DNA_ALPHABET);
  if (!testModel(gtr)) return 1;

  //Codon models:
  StandardGeneticCode gc(&AlphabetTools::DNA_ALPHABET);
  const CodonAlphabet* codonAlphabet = new StandardCodonAlphabet(&AlphabetTools::DNA_ALPHABET);
  FrequenciesSet* fset = CodonFrequenciesSet::getFrequenciesSetForCodons(CodonFrequenciesSet::F3X4, *codonAlphabet);
  YN98 yn98(&gc, fset);
  if (!testModel(yn98)) return 1;

  delete codonAlphabet;

  return 0;
}
