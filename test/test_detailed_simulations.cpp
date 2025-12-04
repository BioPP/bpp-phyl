// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/Matrix/Matrix.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Seq/Alphabet/DNA.h>

#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Model/Nucleotide/GTR.h>

#include <Bpp/Phyl/Simulation/SimpleSubstitutionProcessSequenceSimulator.h>

#include <Bpp/Phyl/Likelihood/RateAcrossSitesSubstitutionProcess.h>

#include <iostream>

using namespace bpp;
using namespace std;

int main()
{
  Newick reader;
  auto phyloTree = std::shared_ptr<PhyloTree>(reader.parenthesisToPhyloTree("((A:0.001, B:0.002):0.003,C:0.01,D:0.1);"));

  auto ids = phyloTree->getAllEdgesIndexes();

  // -------------

  shared_ptr<NucleicAlphabet> alphabet(new DNA());
  auto model = make_shared<GTR>(alphabet, 1, 0.2, 0.3, 0.4, 0.4, 0.1, 0.35, 0.35, 0.2);
  // DiscreteDistribution* rdist = new GammaDiscreteDistribution(4, 0.4, 0.4);
  auto rdist = std::make_shared<ConstantDistribution>(1.0);
  auto process = make_shared<RateAcrossSitesSubstitutionProcess>(model, rdist, phyloTree);

  process->setPhyloTree(*phyloTree);

  SimpleSubstitutionProcessSiteSimulator simulatorS(process);

  unsigned int n = 200000;
  map<unsigned int, RowMatrix<unsigned int>> counts;
  for (size_t j = 0; j < ids.size(); ++j)
  {
    counts[ids[j]].resize(4, 4);
  }
  for (unsigned int i = 0; i < n; ++i)
  {
    auto result = simulatorS.dSimulateSite();
    for (size_t j = 0; j < ids.size(); ++j)
    {
      result->getMutationPath(ids[j]).getEventCounts(counts[ids[j]]);
    }
  }

  const auto& Q = model->generator();

  map<unsigned int, RowMatrix<double>> freqs;
  map<unsigned int, double> sums;
  for (size_t k = 0; k < ids.size(); ++k)
  {
    RowMatrix<double>* freqsP = &freqs[ids[k]];
    RowMatrix<unsigned int>* countsP = &counts[ids[k]];
    freqsP->resize(4, 4);
    for (unsigned int i = 0; i < 4; ++i)
    {
      for (unsigned int j = 0; j < 4; ++j)
      {
        (*freqsP)(i, j) = static_cast<double>((*countsP)(i, j)) / (static_cast<double>(n));
      }
    }

    // For now we simply compare the total number of substitutions:
    sums[ids[k]] = MatrixTools::sumElements(*freqsP);

    cout << "Br" << ids[k] << " BrLen = " << phyloTree->getEdge(ids[k])->getLength() << " counts = " << sums[ids[k]] << endl;
    MatrixTools::print(*freqsP);

    cout << " Comparison with generator (more or less same non-diagonal values on each line):" << endl;
    RowMatrix<double> comp(4, 4);

    for (unsigned int i = 0; i < 4; ++i)
    {
      for (unsigned int j = 0; j < 4; ++j)
      {
        comp(i, j) = (*freqsP)(i, j) / Q(i, j);
      }
    }

    MatrixTools::print(comp);
  }


  for (size_t k = 0; k < ids.size(); ++k)
  {
    if (abs(sums[ids[k]] - phyloTree->getEdge(ids[k])->getLength()) > 0.01)
    {
      return 1;
    }
  }

  // -------------

  // return (abs(obs - 0.001) < 0.001 ? 0 : 1);
  return 0;
}
