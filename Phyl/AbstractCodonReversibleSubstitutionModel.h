//
// File: AbstractCodonReversibleSubstitutionModel.h
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

#ifndef _ABSTRACTCODONREVERSIBLESUBSTITUTIONMODEL_H_
#define _ABSTRACTCODONREVERSIBLESUBSTITUTIONMODEL_H_

#include "AbstractWordReversibleSubstitutionModel.h"
#include "NucleotideSubstitutionModel.h"

// From SeqLib:
#include <Seq/CodonAlphabet.h>

namespace bpp
{

/**
 * @brief Abstract class for reversible substitution models on codons.
 * @author Laurent Guéguen
 *
 * Objects of this class are built from three reversible substitution
 * models of NucleicAlphabets. No model is directly accessible. </p>
 *
 * Only substitutions with one letter changed are accepted. </p>
 *
 * There is one substitution per word per unit of time
 * on the equilibrium frequency, and each position has its specific rate.
 *
 */
  
class AbstractCodonReversibleSubstitutionModel :
  public AbstractWordReversibleSubstitutionModel
{
public:

  /**
   *@brief Build a new AbstractCodonReversibleSubstitutionModel object from
   *a pointer to an AbstractReversibleSubstitutionModel. 
   *
   *@param palph pointer to a CodonAlphabet
   *@param pmod pointer to the AbstractReversibleSubstitutionModel to use in the three positions.
   *@param st string of the Namespace
   */
  
  AbstractCodonReversibleSubstitutionModel(
      const CodonAlphabet* palph,
      NucleotideSubstitutionModel* pmod,
      const std::string& st);
  
  /**
   *@brief Build a new AbstractReversibleSubstitutionModel object
   *from three pointers to AbstractReversibleSubstitutionModels. 
   *
   *@param palph pointer to a CodonAlphabet
   *@param pmod1, pmod2, pmod3 are pointers to the
   *AbstractReversibleSubstitutionModel to use in the three positions.
   * All the models must be different objects to avoid redondant
   * parameters.
   *@param st string of the Namespace
   */

  AbstractCodonReversibleSubstitutionModel(
      const CodonAlphabet*,
      NucleotideSubstitutionModel* pmod1,
      NucleotideSubstitutionModel* pmod2,  
      NucleotideSubstitutionModel* pmod3,
      const std::string& st);

  virtual ~AbstractCodonReversibleSubstitutionModel() {};  
};

} //end of namespace bpp.

#endif	

