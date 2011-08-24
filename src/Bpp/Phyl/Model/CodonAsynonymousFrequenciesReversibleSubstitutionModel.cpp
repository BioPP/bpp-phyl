//
// File: CodonAsynonymousFrequenciesReversibleSubstitutionModel.cpp
// Created by:  Laurent Gueguen
// Created on: Feb 2009
//

/*
   Copyright or © or Copr. CNRS, (November 16, 2004)
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

#include "CodonAsynonymousFrequenciesReversibleSubstitutionModel.h"

using namespace bpp;

using namespace std;

/******************************************************************************/

CodonAsynonymousFrequenciesReversibleSubstitutionModel::CodonAsynonymousFrequenciesReversibleSubstitutionModel(
  const GeneticCode* palph,
  FrequenciesSet* pfreq,
  const AlphabetIndex2<double>* pdist) throw (Exception) :
  AbstractCodonFrequenciesReversibleSubstitutionModel(
    dynamic_cast<const CodonAlphabet*>(palph->getSourceAlphabet()),
    pfreq,
    "CodonAsynonymousFrequencies."),
  geneticCode_(palph),
  pdistance_(pdist)
{
  if (pdistance_)
    addParameter_(Parameter("CodonAsynonymousFrequencies.alpha", 10000, &Parameter::R_PLUS_STAR));

  addParameter_(Parameter("CodonAsynonymousFrequencies.beta", 1, new IncludingInterval(NumConstants::TINY, 999), true));
  updateMatrices();
}

string CodonAsynonymousFrequenciesReversibleSubstitutionModel::getName() const
{
  return "CodonAsynonymousFrequenciesReversibleSubstitutionModel model : " + pfreqset_->getName();
}

void CodonAsynonymousFrequenciesReversibleSubstitutionModel::completeMatrices()
{
  unsigned int i, j;
  unsigned int salph = getNumberOfStates();
  double alpha = pdistance_ ? getParameterValue("alpha") : 1;
  double beta = getParameterValue("beta");
  const CodonAlphabet* ca = dynamic_cast<const CodonAlphabet*>(geneticCode_->getSourceAlphabet());

  for (i = 0; i < salph; i++)
  {
    for (j = 0; j < salph; j++)
    {
      if (i != j)
      {
        if (ca->isStop(j) || ca->isStop(i))
        {
          generator_(i,j) = 0;
          exchangeability_(i,j) = 0;
        }
        else
        {
          if (!geneticCode_->areSynonymous(i,j))
          {
            generator_(i,j) *= beta * (pdistance_ ? exp(-pdistance_->getIndex(geneticCode_->translate(i), geneticCode_->translate(j)) / alpha) : 1);
            exchangeability_(i,j) *= beta * (pdistance_ ? exp(-pdistance_->getIndex(geneticCode_->translate(i), geneticCode_->translate(j)) / alpha) : 1);
          }
        }
      }
    }
  }

  AbstractCodonFrequenciesReversibleSubstitutionModel::completeMatrices();
}

