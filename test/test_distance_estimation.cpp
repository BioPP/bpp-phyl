// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Matrix/MatrixTools.h>

#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>

#include <Bpp/Phyl/Model/Nucleotide/K80.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Distance/DistanceEstimation.h>

#include <Bpp/Phyl/Likelihood/RateAcrossSitesSubstitutionProcess.h>

#include <iostream>

using namespace bpp;
using namespace std;

int main()
{
  // Note: this only tests instantiation and exceptions. The correctness of calculations is not assessed.

  try
  {
    // Creates a sequence alignment, with missing data:

    shared_ptr<const NucleicAlphabet> nucAlphabet = AlphabetTools::DNA_ALPHABET;
    shared_ptr<const Alphabet> alphabet = AlphabetTools::DNA_ALPHABET;

    
    auto sites = make_shared<VectorSiteContainer>(alphabet);
    auto seqA = make_unique<Sequence>("A", "ATCCAGACATGCCGGGACTTTGCAGAGAAGGAGTTGTTTCCCATTGCAGCCCAGGTGGATAAGGAACAGC", alphabet);
    sites->addSequence("A", seqA);
    auto seqB = make_unique<Sequence>("B", "CGTCAGACATGCCGTGACTTTGCCGAGAAGGAGTTGGTCCCCATTGCGGCCCAGCTGGACAGGGAGCATC", alphabet);
    sites->addSequence("B", seqB);
    auto seqC = make_unique<Sequence>("C", "GGTCAGACATGAAGGGAATTTGCTGGTCAGGAGCTGGTTCCCATTGCAGCCCAGGTAGACAAGGAGCATC", alphabet);
    sites->addSequence("C", seqC);
    auto seqD = make_unique<Sequence>("D", "TTCCAGACATGCCGGGACTTTACCGAGAAGGAGTTGTTTTCCATTGCAGCCCAGGTGGATAAGGAACATC", alphabet);
    sites->addSequence("D", seqD);

    // SiteContainerTools::changeGapsToUnknownCharacters(*sites);

    // Test with a JC model:
    auto model1 = make_shared<K80>(nucAlphabet);
    auto rdist = make_shared<ConstantRateDistribution>();
    shared_ptr<ParametrizablePhyloTree> pt(nullptr);
    auto process = std::make_shared<RateAcrossSitesSubstitutionProcess>(model1, rdist, pt);
    DistanceEstimation de(process, sites, 0, true);
    cout << endl;
    VectorTools::printForR(de.getMatrix()->getNames(),"Names",cout);
    MatrixTools::printForR(de.getMatrix()->asMatrix(),"Dist",cout);
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}
