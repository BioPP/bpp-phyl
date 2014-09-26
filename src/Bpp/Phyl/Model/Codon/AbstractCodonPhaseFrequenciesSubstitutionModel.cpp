//
// File: AbstractCodonPhaseFrequenciesSubstitutionModel.cpp
// Created by:  Laurent Gueguen
// Created on: vendredi 23 septembre 2011, à 16h 29
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

#include "AbstractCodonPhaseFrequenciesSubstitutionModel.h"
#include "../FrequenciesSet/NucleotideFrequenciesSet.h"

using namespace bpp;

using namespace std;

/******************************************************************************/

AbstractCodonPhaseFrequenciesSubstitutionModel::AbstractCodonPhaseFrequenciesSubstitutionModel(
  FrequenciesSet* pfreq,
  const std::string& prefix) :
  CodonSubstitutionModel(),
  AbstractParameterAliasable(prefix),
  posfreqset_(),
  freqName_("")
{
  CodonFrequenciesSet* pCFS = dynamic_cast<CodonFrequenciesSet*>(pfreq);
  if (!pCFS)
    throw Exception("Bad type for equilibrium frequencies " + pfreq->getName());

  if (dynamic_cast<CodonFromUniqueFrequenciesSet*>(pCFS)
   || dynamic_cast<CodonFromIndependentFrequenciesSet*>(pCFS))
    posfreqset_ = dynamic_cast<WordFrequenciesSet*>(pfreq)->clone();
  else
  {
    vector<FrequenciesSet*> vFS;
    if (dynamic_cast<FixedCodonFrequenciesSet*>(pCFS)) {
      for (unsigned int i = 0; i < 3; i++)
      {
        vFS.push_back(new FixedNucleotideFrequenciesSet(pCFS->getAlphabet()->getNucleicAlphabet()));
      }
    } else {
      for (unsigned int i = 0; i < 3; i++)
      {
        vFS.push_back(new FullNucleotideFrequenciesSet(pCFS->getAlphabet()->getNucleicAlphabet()));
      }
    }
    posfreqset_ = new CodonFromIndependentFrequenciesSet(
        pCFS->getGeneticCode(),
        vFS, "");

    posfreqset_->setFrequencies(pfreq->getFrequencies());
  }

  freqName_ = pfreq->getNamespace();
  posfreqset_->setNamespace(prefix + pfreq->getNamespace());
  //  if (dynamic_cast<FixedCodonFrequenciesSet*>(pCFS)!=NULL)
  addParameters_(posfreqset_->getParameters());
}

AbstractCodonPhaseFrequenciesSubstitutionModel::~AbstractCodonPhaseFrequenciesSubstitutionModel()
{
  if (posfreqset_)
    delete posfreqset_;
}

void AbstractCodonPhaseFrequenciesSubstitutionModel::fireParameterChanged(const ParameterList& parameters)
{
  posfreqset_->matchParametersValues(parameters);
}


void AbstractCodonPhaseFrequenciesSubstitutionModel::setFreq(map<int, double>& frequencies)
{
  posfreqset_->setFrequenciesFromAlphabetStatesFrequencies(frequencies);
  matchParametersValues(posfreqset_->getParameters());
}

double AbstractCodonPhaseFrequenciesSubstitutionModel::getCodonsMulRate(size_t i, size_t j) const
{
  size_t i2(i), j2(j);

  double x = 1.;
  for (size_t k = 0; k < 3; k++)
  {
    if ((i2 % 4) != (j2 % 4))
      x *= posfreqset_->getFrequenciesSetForLetter(2 - k).getFrequencies()[j2 % 4];
    i2 /= 4;
    j2 /= 4;
  }
  return x;
}

