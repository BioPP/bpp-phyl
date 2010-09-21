//
// File: CodonNeutralReversibleSubstitutionModel.h
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

#ifndef _CODONNEUTRALREVERSIBLESUBSTITUTIONMODEL_H_
#define _CODONNEUTRALREVERSIBLESUBSTITUTIONMODEL_H_

#include "AbstractCodonReversibleSubstitutionModel.h"
#include "NucleotideSubstitutionModel.h"

// From SeqLib:
#include <Seq/CodonAlphabet.h>

namespace bpp
{
/**
 * @brief Class for reversible substitution models on non stop codons.
 *
 * Objects of this class are built from three reversible substitution
 * models of NucleicAlphabets. No model is directly accessible. </p>
 *
 * Only substitutions with one letter changed are accepted. </p>
 *
 * There is one substitution per word per unit of time
 * on the equilibrium frequency, and each position has its specific rate.
 *
 * The generator is constructed in two steps:
 * First, if there are @f$n@f$ models and @f$\rho_i@f$ is the rate of
 * model i (@f$\sum_{i=1}^{n} \rho_i = 1@f$):
 * @f[
 * Q_{abc \rightarrow abd} = \rho_3 Q^{(2)}_{c \rightarrow d}
 * Q_{abc \rightarrow aed} = 0
 * Q_{abc \rightarrow abc} = \rho_1 Q^{(0)}_{a \rightarrow a} + \rho_2 Q^{(1)}_{b \rightarrow b} + \rho_3 Q^{(2)}_{c \rightarrow c})
 * @f]
 *
 * Second, the subsitution between STOP codons and non-STOP codons are
 * nullified, and the lines of Q are normalized in accordance.
 *
 * The parameters of this word model are the same as the ones of the
 * models used. Their names have a new suffix, "_phi" where i stands
 * for the position (i.e. the phase) in the word.
 *
 * The rates are defined by relative rates parameters @f$r_i@f$
 * (called "relrate_i") with:
 * @f[
 * 1 <= i < n, \rho_i = (1-r_1).(1-r_2)...(1-r_{i-1}).r_{i}
 * \rho_n = (1-r_1).(1-r_2)...(1-r_{n-1})
 * @f]
 * and
 * @f[
 * \forall 1 <= i < n, r_i = \frac{\rho_i}{1-(\rho_0+...\rho_{i-1})}
 * @f]
 *
 * The parameters of this codon model are the same as the ones of the
 * models used. Their names have a new suffix, "_phi" where i stands
 * for the position (i.e. the phase) in the word.
 */

class CodonNeutralReversibleSubstitutionModel :
  public AbstractCodonReversibleSubstitutionModel
{
public:
  /**
   *@brief Build a new CodonNeutralReversibleSubstitutionModel object from
   *a pointer to NucleotideSubstitutionModels.
   * @author Laurent Guéguen
   *
   *@param palph pointer to a CodonAlphabet
   *@param pmod1 pointer to the NucleotideSubstitutionModel to use in the
   *       three positions. It is owned by the instabce.
   */

  CodonNeutralReversibleSubstitutionModel(const CodonAlphabet* palph,
                                          NucleotideSubstitutionModel* pmod1);

  /**
   *@brief Build a new CodonNeutralReversibleSubstitutionModel object
   *from three pointers to NucleotideSubstitutionModels.
   *
   *@param palph pointer to a CodonAlphabet
   *@param pmod1, pmod2, pmod3 pointers to the
   *   NucleotideSubstitutionModel to use in the three positions.
   *   All the models must be different objects to avoid parameters
   *   redondancy, otherwise only the first model is used. The used models
   *   are owned by the instance.
   */

  CodonNeutralReversibleSubstitutionModel(const CodonAlphabet* palph,
                                          NucleotideSubstitutionModel* pmod1,
                                          NucleotideSubstitutionModel* pmod2,
                                          NucleotideSubstitutionModel* pmod3);

  ~CodonNeutralReversibleSubstitutionModel(){}

#ifndef NO_VIRTUAL_COV
  CodonNeutralReversibleSubstitutionModel*
#else
  Clonable*
#endif
  clone() const { return new CodonNeutralReversibleSubstitutionModel(*this); }

public:
  void completeMatrices();
  void updateMatrices();
  std::string getName() const;
};
} // end of namespace bpp.

#endif

