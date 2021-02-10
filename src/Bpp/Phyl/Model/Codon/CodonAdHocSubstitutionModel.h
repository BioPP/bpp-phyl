//
// File: CodonAdHocSubstitutionModel.h
// Created by: Laurent Gueguen
// Created on: lundi 30 octobre 2017, à 06h 31
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

#ifndef _CODON_ADHOC_SUBSTITUTIONMODEL_H_
#define _CODON_ADHOC_SUBSTITUTIONMODEL_H_

#include "AbstractCodonSubstitutionModel.h"

namespace bpp
{
  /**
   * @brief Class for substitution models of codons with
   * several layers of codon models
   *
   * @author Laurent Guéguen
   *
   * Objects of this class are built from three substitution models of
   * NucleicAlphabets. No model is directly accessible. </p>
   *
   * Only substitutions with one letter changed are accepted. </p>
   *
   *
   */

  class CodonAdHocSubstitutionModel :
    public AbstractCodonSubstitutionModel
  {
  private:

    std::vector<std::unique_ptr<CoreCodonSubstitutionModel> > vModel_;

    std::string name_;

    /*
     * @brief optional FrequencySet if model is defined through a
     * FrequencySet.
     *
     */
     
    std::shared_ptr<FrequencySet> freqSet_;
    
  public:
    /**
     * @brief Build a new CodonAdHocSubstitutionModel object from
     * a pointer to NucleotideSubstitutionModel.
     *
     * @param gCode pointer to a GeneticCode
     * @param pmod  pointer to the NucleotideSubstitutionModel to use
     * in the three positions.
     * The instance will then own this substitution model.
     * @param vpmodel vector of codon models. They will be owned by the model.
     * @param name the name of the model
     */

    CodonAdHocSubstitutionModel(
      const GeneticCode* gCode,
      NucleotideSubstitutionModel* pmod,
      std::vector<CoreCodonSubstitutionModel*>& vpmodel,
      const std::string& name);

    /**
     * @brief Build a new CodonAdHocSubstitutionModel object
     * from three pointers to NucleotideSubstitutionModels.
     *
     * @param gCode pointer to a GeneticCode
     * @param pmod1, pmod2, pmod3 pointers to the
     *   NucleotideSubstitutionModels to use in the three positions.
     *   Either all the models are different objects to avoid parameters
     *   redondancy, or only the first model is used in every position.
     *   The used models are owned by the instance.
     * @param vpmodel vector of codon models. They will be owned by the
     * model.
     * @param name the name of the model
     */

    CodonAdHocSubstitutionModel(
        const GeneticCode* gCode,
        NucleotideSubstitutionModel* pmod1,
        NucleotideSubstitutionModel* pmod2,
        NucleotideSubstitutionModel* pmod3,
        std::vector<CoreCodonSubstitutionModel*>& vpmodel,
        const std::string& name);

    CodonAdHocSubstitutionModel(const CodonAdHocSubstitutionModel& model);

    CodonAdHocSubstitutionModel& operator=(const CodonAdHocSubstitutionModel& model);

    virtual ~CodonAdHocSubstitutionModel() {}

    CodonAdHocSubstitutionModel* clone() const
    {
      return new CodonAdHocSubstitutionModel(*this);
    }

  public:
    void fireParameterChanged(const ParameterList& parameterlist);
  
    std::string getName() const
    {
      return name_;
    }

    void setNamespace(const std::string& prefix){
      AbstractCodonSubstitutionModel::setNamespace(prefix);
      for (auto& model : vModel_)
        model->setNamespace(prefix);
    }

    size_t  getNumberOfModels() const
    {
      return vModel_.size();
    }

    const std::unique_ptr<CoreCodonSubstitutionModel>& getNModel(size_t i) const
    {
      return vModel_[i];
    }
    
    double getCodonsMulRate(size_t i, size_t j) const;

    void setFreq(std::map<int,double>& frequencies);

    const std::shared_ptr<FrequencySet> getFrequencySet() const
    {
      return freqSet_;
    }
    
  };
  
} // end of namespace bpp.

#endif //_CODON_ADHOC_SUBSTITUTIONMODEL_H_

