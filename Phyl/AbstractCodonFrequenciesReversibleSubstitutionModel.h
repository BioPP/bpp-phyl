//
// File: AbstractCodonFrequenciesReversibleSubstitutionModel.h
// Created by: Laurent Gueguen
// Created on: Tue Dec 24 11:03:53 2003
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

#ifndef _ABSTRACTCODONFREQUENCIESREVERSIBLESUBSTITUTIONMODEL_H_
#define _ABSTRACTCODONFREQUENCIESREVERSIBLESUBSTITUTIONMODEL_H_

#include "AbstractWordReversibleSubstitutionModel.h"
#include "FrequenciesSet.h"

// From SeqLib:
#include <Seq/CodonAlphabet.h>

namespace bpp
{

/**
 * @brief Class for reversible substitution models on codons
 *  parametrized by the equilibrium frequencies. The basic
 *  substitution model is the K80 model. </p>
 *
 * Objects of this class are built from three reversible substitution
 *  models of NucleicAlphabets. No model is directly accessible. </p>
 *
 * Only substitutions with one letter changed are accepted. </p>
 *
 * There is one substitution per word per unit of time
 * on the equilibrium frequency, and each position has its specific rate.</p>
 *
 */
  
class AbstractCodonFrequenciesReversibleSubstitutionModel :
  public AbstractWordReversibleSubstitutionModel
{

protected:

  AbstractFrequenciesSet* AbsFreq;
  
public:

  /**
   *@brief Build a new
   *AbstractCodonFrequenciesReversibleSubstitutionModel object from
   *three instances of the K80 model.
   *
   *@param pointer to a CodonAlphabet
   *@param pointer to the AbstractFrequenciesSet* equilibrium frequencies
   */
  
  AbstractCodonFrequenciesReversibleSubstitutionModel(const CodonAlphabet*,
                                                      AbstractFrequenciesSet*,
                                                      const std::string&) throw(Exception);
  
  AbstractCodonFrequenciesReversibleSubstitutionModel(const AbstractCodonFrequenciesReversibleSubstitutionModel&);

  virtual ~AbstractCodonFrequenciesReversibleSubstitutionModel();

  virtual void fireParameterChanged(const ParameterList&);

public:

  virtual void completeMatrices();
  
};

} //end of namespace bpp.

#endif	

