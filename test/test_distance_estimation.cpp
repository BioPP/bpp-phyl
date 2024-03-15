// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>

#include <Bpp/Phyl/Model/Nucleotide/JCnuc.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Distance/DistanceEstimation.h>

#include <iostream>

using namespace bpp;
using namespace std;

int main() {

  //Note: this only tests instanciation and exceptions. The correctness of calculations is not assessed.
  
  try {
  
    // Creates a sequence alignment, with missing data:
    
    shared_ptr<const NucleicAlphabet> nucAlphabet = AlphabetTools::DNA_ALPHABET;
    shared_ptr<const Alphabet> alphabet = AlphabetTools::DNA_ALPHABET;

    auto sites = make_shared<VectorSiteContainer>(alphabet);
    auto seqA = make_unique<Sequence>("A", "ATCCAGACATGCCGGGACTTTGCAGAGAAGGAGTTGTTTCCCATTGCAGCCCAGGTGGATAAGGAACAGC", alphabet);
    sites->addSequence("A", seqA);
    auto seqB = make_unique<Sequence>("B", "CGTCAGACATGCCGTGACTTTGCCGAGAAGGAGTTGGTCCCCATTGCGGCCCAGCTGGACAGGGAGCATC", alphabet);
    sites->addSequence("B", seqB);
    auto seqC = make_unique<Sequence>("C", "GGTCAGACATGCCGGGAATTTGCTGAAAAGGAGCTGGTTCCCATTGCAGCCCAGGTAGACAAGGAGCATC", alphabet);
    sites->addSequence("C", seqC);
    auto seqD = make_unique<Sequence>("D", "TTCCAGACATGCCGGGACTTTACCGAGAAGGAGTTGTTTTCCATTGCAGCCCAGGTGGATAAGGAACATC", alphabet);
    sites->addSequence("D", seqD);

    //SiteContainerTools::changeGapsToUnknownCharacters(*sites);

    // Test with a JC model:
    auto model1 = make_shared<JCnuc>(nucAlphabet);
    auto rdist = make_shared<ConstantRateDistribution>();
    DistanceEstimation de(model1, rdist, sites, 1, true);

  } catch (exception& e) {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}
