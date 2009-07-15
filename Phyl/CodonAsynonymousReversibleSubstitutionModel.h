//
// File: CodonAsynonymousReversibleSubstitutionModel.h
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

#ifndef _CODONASYNONYMOUSREVERSIBLESUBSTITUTIONMODEL_H_
#define _CODONASYNONYMOUSREVERSIBLESUBSTITUTIONMODEL_H_

#include "AbstractCodonReversibleSubstitutionModel.h"
#include "NucleotideSubstitutionModel.h"

// From SeqLib:
#include <Seq/CodonAlphabet.h>
#include <Seq/GeneticCode.h>
#include <Seq/AlphabetIndex2.h>

namespace bpp
{

/**
 * @brief Class for asynonymous substitution models on codons.
 *
 * Objects of this class are built from three reversible substitution
 * models of NucleicAlphabets. No model is directly accessible. </p>
 *
 * Only substitutions with one letter changed are accepted. </p>
 *
 *
 * If a distance $d$ between amino-acids is defined, the ratio between
 * non-synonymous and synonymous substitutions rates is, if the codied
 * amino-acids are $x$ and $y$, $\beta*\exp(-\alpha.d(x,y))$ with
 * non-negative parameter $\alpha$ and positive parameter $\beta$.
 *
 * If such a distance is not defined, the ratio between non-synonymous
 * and synonymous substitutions rates is $\beta$ with positive
 * parameter $\beta$.
 */
  
class CodonAsynonymousReversibleSubstitutionModel :
  public AbstractCodonReversibleSubstitutionModel
{
private:

  const GeneticCode* _geneticCode;
  const AlphabetIndex2<double>*   _pdistance;
  
public:

  /**
   *@brief Build a new CodonNeutralReversibleSubstitutionModel object from
   *a pointer to AbstractReversibleSubstitutionModels. 
   *
   *@param pointer to a CodonAlphabet
   *@param pmodel is a pointer to the
   *AbstractReversibleSubstitutionModel to use in the three positions.
   *@param optional pointer to a distance between amino-acids
   */
  
  CodonAsynonymousReversibleSubstitutionModel(const GeneticCode*,
                                              NucleotideSubstitutionModel*,
                                              const AlphabetIndex2<double>* = 0);
  
  /**
   *@brief Build a new CodonNeutralReversibleSubstitutionModel object
   *from three pointers to AbstractReversibleSubstitutionModels. 
   *
   *@param pointer to a CodonAlphabet @param pmodel1, pmodel2 and
   * pmodel3 are pointers to the AbstractReversibleSubstitutionModel
   * to use in the three positions. All the models must be different
   * objects to avoid parameters redondancy, otherwise only the first
   * model is used.
   *@param optional pointer to a distance between *amino-acids
   */

  CodonAsynonymousReversibleSubstitutionModel(const GeneticCode*,
                                              NucleotideSubstitutionModel*,
                                              NucleotideSubstitutionModel*, 
                                              NucleotideSubstitutionModel*,
                                              const AlphabetIndex2<double>* = 0);

  CodonAsynonymousReversibleSubstitutionModel(const CodonAsynonymousReversibleSubstitutionModel&);


  ~CodonAsynonymousReversibleSubstitutionModel() {};
  
#ifndef NO_VIRTUAL_COV
  CodonAsynonymousReversibleSubstitutionModel*
#else
  Clonable*
#endif
  clone() const { return new CodonAsynonymousReversibleSubstitutionModel(*this);}
  
public:
  void completeMatrices();

  string getName() const;
};

} //end of namespace bpp.

#endif	

