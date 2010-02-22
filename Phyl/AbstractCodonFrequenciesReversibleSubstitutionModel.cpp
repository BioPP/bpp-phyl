//
// File: AbstractCodonFrequenciesReversibleSubstitutionModel.cpp
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

#include "AbstractCodonFrequenciesReversibleSubstitutionModel.h"
#include "K80.h"

#include <Seq/AlphabetTools.h>

using namespace bpp;

using namespace std;

/******************************************************************************/

AbstractCodonFrequenciesReversibleSubstitutionModel::AbstractCodonFrequenciesReversibleSubstitutionModel(
  const CodonAlphabet* palph,
  FrequenciesSet* pfreq,
  const std::string& st) throw (Exception) :
  AbstractWordReversibleSubstitutionModel(palph,st),
  pfreqset_(pfreq)
{
  enableEigenDecomposition(1);

  rate_.resize(3);

  SubstitutionModel* pmodel = new K80(palph->getNucleicAlphabet());

  for (unsigned i = 0; i < 3; i++)
  {
   VSubMod_.push_back(pmodel);
   VnestedPrefix_.push_back(pmodel->getNamespace());
    rate_[i] = 1.0 / 3;
  }

  pmodel->setNamespace(st + "123_" + VnestedPrefix_[0]);
  addParameters_(pmodel->getParameters());

  if (pfreqset_->getAlphabet()->getSize() != 64)
    throw Exception("Bad Alphabet for equilibrium frequencies " + pfreqset_->getAlphabet()->getAlphabetType());

  pfreqset_->setNamespace(st + pfreqset_->getNamespace());
  addParameters_(pfreqset_->getParameters());
}

AbstractCodonFrequenciesReversibleSubstitutionModel::~AbstractCodonFrequenciesReversibleSubstitutionModel()
{
  if (pfreqset_)
    delete pfreqset_;
  pfreqset_ = 0;
}

void AbstractCodonFrequenciesReversibleSubstitutionModel::fireParameterChanged(const ParameterList& parameters)
{
   pfreqset_->matchParametersValues(parameters);
   AbstractWordReversibleSubstitutionModel::fireParameterChanged(parameters);
}


void AbstractCodonFrequenciesReversibleSubstitutionModel::setFreq(map<int,double>& frequencies)
{
   pfreqset_->setFrequenciesFromMap(frequencies);
  updateMatrices();
}

void AbstractCodonFrequenciesReversibleSubstitutionModel::completeMatrices()
{
   unsigned int i, j;
   unsigned int salph = getNumberOfStates();

  freq_ = pfreqset_->getFrequencies();

  for (i = 0; i < salph; i++)
  {
    for (j = 0; j < salph; j++)
    {
      generator_(i, j) = exchangeability_(i, j) * freq_[j];
    }
  }
}

