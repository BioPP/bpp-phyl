// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Phyl/Model/Nucleotide/GTR.h>
#include <Bpp/Phyl/Model/Codon/YN98.h>
#include <Bpp/Phyl/Model/Protein/JTT92.h>
#include <Bpp/Phyl/Model/POMO.h>
#include <Bpp/Phyl/Model/FrequencySet/CodonFrequencySet.h>
#include <Bpp/Phyl/Model/FrequencySet/NucleotideFrequencySet.h>
#include <Bpp/Phyl/Model/FrequencySet/ProteinFrequencySet.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Alphabet/CodonAlphabet.h>
#include <Bpp/Seq/Alphabet/ProteicAlphabet.h>
#include <Bpp/Seq/Alphabet/AllelicAlphabet.h>
#include <Bpp/Seq/GeneticCode/StandardGeneticCode.h>
#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/Function/ReparametrizationFunctionWrapper.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Numeric/AbstractParametrizable.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <iostream>

using namespace bpp;
using namespace std;

class DummyFunction :
  public virtual FunctionInterface,
  public AbstractParametrizable
{
public:
  DummyFunction(const SubstitutionModelInterface& model) :
    AbstractParametrizable("")
  {
    addParameters_(model.getParameters());
  }

  DummyFunction* clone() const { return new DummyFunction(*this); }

  void setParameters(const ParameterList& pl)
  {
    matchParametersValues(pl);
  }

  double getValue() const { return 0; }

  void fireParameterChanged(const bpp::ParameterList&) {}
};

bool testModel(SubstitutionModelInterface& model)
{
  ParameterList pl = model.getParameters();
  auto df = make_shared<DummyFunction>(model);
  ReparametrizationFunctionWrapper redf(df, pl, false);

  // Try to modify randomly each parameter and check that the new parameter apply correctly:
  for (unsigned int i = 0; i < 10; ++i)
  {
    // Get random parameters (unconstrained):
    ParameterList pl2 = redf.getParameters();
    for (size_t j = 0; j < pl.size(); ++j)
    {
      double value = (RandomTools::flipCoin() ? 1. : -1) * RandomTools::giveRandomNumberBetweenZeroAndEntry(1);
      pl2[j].setValue(value);
    }
    // Apply unconstrained parameters:
    redf.setParameters(pl2);
    // Retrieve transformed parameters:
    pl2 = df->getParameters();
    // pl2.printParameters(cout);
    // Now apply the new parameters and retrieve them again:
    model.matchParametersValues(pl2);
    ParameterList pl3 = model.getParameters();
    // Compare the two lists:
    for (size_t j = 0; j < pl.size(); ++j)
    {
      if (abs(pl2[j].getValue() - pl3[j].getValue()) > 0.0000001)
      {
        cerr << "ERROR for parameter " << pl2[j].getName() << ": " << pl2[j].getValue() << "<>" << pl3[j].getValue() << endl;
        return false;
      }
    }
  }

  return true;
}

int main()
{
  // Nucleotide models:
  auto gtr = std::make_unique<GTR>(AlphabetTools::DNA_ALPHABET);

  if (!testModel(*gtr))
    return 1;

  // Codon models:
  auto gc = make_shared<StandardGeneticCode>(AlphabetTools::DNA_ALPHABET);
  auto fset = CodonFrequencySetInterface::getFrequencySetForCodons(CodonFrequencySetInterface::F3X4, gc);
  YN98 yn98(gc, std::move(fset));

  if (!testModel(yn98))
    return 1;

  // Allelic models

  auto allalph = make_shared<AllelicAlphabet>(AlphabetTools::DNA_ALPHABET, 4);
  auto fit = std::make_unique<FullNucleotideFrequencySet>(AlphabetTools::DNA_ALPHABET);

  auto statemod = std::move(gtr);

  // AllelicAlphabet allalph(AlphabetTools::PROTEIN_ALPHABET, 4);
  // auto fit = std::make_shared<FullProteinFrequencySet>(&AlphabetTools::PROTEIN_ALPHABET);
  // auto freq = std::make_shared<FullProteinFrequencySet>(&AlphabetTools::PROTEIN_ALPHABET);

  // auto statemod = std::make_shared<JTT92>(&AlphabetTools::PROTEIN_ALPHABET, freq);

  POMO pomo(allalph, std::move(statemod), std::move(fit));

  auto& Q = pomo.generator();

  MatrixTools::printForR(Q, "Q", cerr);

  VectorTools::printForR(pomo.getFrequencies(), "freq", cerr);

  return 0;
}
