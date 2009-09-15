//
// File: CodonReversibleSubstitutionModel.cpp
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

#include "CodonNeutralReversibleSubstitutionModel.h"


using namespace bpp;

using namespace std;

/******************************************************************************/

CodonNeutralReversibleSubstitutionModel::CodonNeutralReversibleSubstitutionModel(const CodonAlphabet* palph,
                                                                                 NucleotideSubstitutionModel* pmod1) : AbstractCodonReversibleSubstitutionModel(palph, pmod1,"CodonNeutral.") 
{
  int i;
  
  // relative rates
  for (i=0; i< 2; i++){
    addParameter_(Parameter("CodonNeutral.relrate"+TextTools::toString(i) , 1.0/(3-i),&Parameter::PROP_CONSTRAINT_EX));
  }
  
  updateMatrices();
}

CodonNeutralReversibleSubstitutionModel::CodonNeutralReversibleSubstitutionModel(const CodonAlphabet* palph,
                                                                                 NucleotideSubstitutionModel* pmod1,
                                                                                 NucleotideSubstitutionModel* pmod2,
                                                                                 NucleotideSubstitutionModel* pmod3) : AbstractCodonReversibleSubstitutionModel(palph, pmod1, pmod2, pmod3, "CodonNeutral.") 
{
  int i;
  
  // relative rates
  for (i=0; i< 2; i++){
    addParameter_(Parameter("CodonNeutral.relrate"+TextTools::toString(i) , 1.0/(3-i),&Parameter::PROP_CONSTRAINT_EX));
  }
  
  updateMatrices();
}

string CodonNeutralReversibleSubstitutionModel::getName() const
{
  string s = "CodonNeutralReversibleSubstitutionModel model:";
  for (unsigned int i = 0; i < _VSubMod.size(); i++)
    s += " "+ _VSubMod[i]->getName();
  
  return s;
}

void CodonNeutralReversibleSubstitutionModel::completeMatrices()
{
  unsigned int i, j;
  unsigned int salph = getNumberOfStates();

  CodonAlphabet* ca = (CodonAlphabet*)(alphabet_);
  
  for (i = 0; i < salph; i++)
    for (j = 0; j < salph; j++)
      if (ca->isStop(i) || ca->isStop(j))
      {
        generator_(i,j) = 0;
      }
}


void CodonNeutralReversibleSubstitutionModel::updateMatrices()
{
  unsigned int i, k, nbmod = _VSubMod.size();
  double x;
  for (k = nbmod - 1; k >= 0; k--)
  {
    x = 1.0;
    for (i = 0; i < k; i++)
      x *= 1 - getParameterValue("relrate" + TextTools::toString(i));
    if (k != nbmod - 1)
      x *= getParameterValue("relrate"+TextTools::toString(k));
    _rate[k] = x;
  }

  AbstractCodonReversibleSubstitutionModel::updateMatrices();
}
