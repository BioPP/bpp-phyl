//
// File: CodonRateSubstitutionModel.h
// Created by: Laurent Gueguen
// Created on: Tue Dec 24 11:03:53 2003
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

#ifndef _CODONRATEFREQUENCIESSUBSTITUTIONMODEL_H_
#define _CODONRATEFREQUENCIESSUBSTITUTIONMODEL_H_

#include "AbstractCodonSubstitutionModel.h"
#include "AbstractCodonFrequenciesSubstitutionModel.h"

namespace bpp
{
/**
 * @brief Class for substitution models on non stop codons, with
 * different parametrized rates on the models, depending on their
 * phase and defined through a specific equilibrium distribution.
 *
 * @author Laurent Guéguen
 *
 * This class should be used with models which equilibrium
 * distribution is fixed, ans does not depend on the parameters.
 * Otherwise there may be problems of identifiability of the
 * parameters.
 *
 * See description in AbstractCodonRateSubstitutionModel and
 * AbstractCodonFrequenciesSubstitutionModel class.
 *
 */
class CodonRateFrequenciesSubstitutionModel :
    public AbstractCodonSubstitutionModel,
    public AbstractCodonFrequenciesSubstitutionModel
{
public:
  
  /**
   * @brief Build a new CodonRateSubstitutionModel object from
   *  a pointer to NucleotideSubstitutionModels.
   * @author Laurent Guéguen
   *
   * @param gCode pointer to a genetic code, which will be owned by this instance.
   * @param pmod pointer to the NucleotideSubstitutionModel to use in the
   *       three positions. It is owned by the instabce.
   * @param pfreq pointer to the FrequenciesSet* equilibrium frequencies
   */
  CodonRateFrequenciesSubstitutionModel(
      const GeneticCode* gCode,
      NucleotideSubstitutionModel* pmod,
      FrequenciesSet* pfreq);

  /**
   * @brief Build a new CodonRateSubstitutionModel object
   *  from three pointers to NucleotideSubstitutionModels.
   *
   * @param gCode pointer to a genetic code, which will be owned by this instance.
   * @param pmod1, pmod2, pmod3 pointers to the
   *   NucleotideSubstitutionModel to use in the three positions.
   *   All the models must be different objects to avoid parameters
   *   redondancy, otherwise only the first model is used. The used models
   *   are owned by the instance.
   * @param pfreq pointer to the FrequenciesSet* equilibrium frequencies
   */

  CodonRateFrequenciesSubstitutionModel(
      const GeneticCode* gCode,
      NucleotideSubstitutionModel* pmod1,
      NucleotideSubstitutionModel* pmod2,
      NucleotideSubstitutionModel* pmod3,
      FrequenciesSet* pfreq);

  virtual ~CodonRateFrequenciesSubstitutionModel(){}

#ifndef NO_VIRTUAL_COV
  CodonRateFrequenciesSubstitutionModel*
#else
  Clonable*
#endif
  clone() const { return new CodonRateFrequenciesSubstitutionModel(*this); }

public:
  void fireParameterChanged(const ParameterList& parameterlist);
  
  std::string getName() const;

  double getCodonsMulRate(size_t i, size_t j) const;

  void setNamespace(const std::string& st);

  void setFreq(std::map<int,double>& frequencies);
};

} // end of namespace bpp.

#endif

