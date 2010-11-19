//
// File: MixedSubstitutionModel.h
// Created by: Laurent Gueguen
// On: vendredi 19 novembre 2010, à 15h 48
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

#ifndef _ABSTRACTMIXEDSUBSTITUTIONMODEL_H_
#define _ABSTRACTMIXEDSUBSTITUTIONMODEL_H_

#include "MixedSubstitutionModel.h"
#include <Bpp/Seq/Alphabet.all>

#include <vector>
#include <string>
#include <map>
#include <cstring> // C lib for string copy

namespace bpp
{
  /**
   * @brief Partial implementation for Mixed Substitution models,
   defined as a mixture of "simple" substitution models.
   * @author Laurent Guéguen
   *
   */

  class AbstractMixedSubstitutionModel :
    public MixedSubstitutionModel
  {
  protected:

    /*
     * @brief vector of pointers to SubstitutionModels.
     *
     * Beware: these SubstitutionModels are owned by the object, so
     * will be deleted at destruction
     *
     */
    
    std::vector<SubstitutionModel*> modelsContainer_;

    std::vector<double> Vprobas_;
    
  public:

    AbstractMixedSubstitutionModel(const Alphabet*, const std::string& prefix);
    
    AbstractMixedSubstitutionModel(const AbstractMixedSubstitutionModel&);
  
    AbstractMixedSubstitutionModel& operator=(const AbstractMixedSubstitutionModel&);

    virtual ~AbstractMixedSubstitutionModel();

    virtual AbstractMixedSubstitutionModel* clone() const = 0;

  public:

    /**
     * @brief Returns a specific model from the mixture
     */
    const SubstitutionModel* getNModel(unsigned int i) const
    {
      return modelsContainer_[i];
    }
    
    SubstitutionModel* getNModel(unsigned int i)
    {
      return modelsContainer_[i];
    }

    /**
     * @brief returns the number of models in the mixture
     */
    
    unsigned int getNumberOfModels() const
    {
      return modelsContainer_.size();
    }
 
    /**
     * @brief Returns the  probability of a specific model from the mixture
     */
  
    double getNProbability(unsigned int i) const
    {
      return Vprobas_[i];
    }
    
    const std::vector<double>& getProbabilities() const
    {
      return Vprobas_;
    }
    
    /**
     * @brief From SubstitutionModel interface
     *
     */

    unsigned int getNumberOfStates() const;
    
    const Matrix<double>& getPij_t(double t) const;
    const Matrix<double>& getdPij_dt(double t) const;
    const Matrix<double>& getd2Pij_dt2(double t) const;
    const Vdouble& getFrequencies();
    double freq(unsigned int i) const;

  };
} // end of namespace bpp.

#endif  // _ABSTRACTMIXEDSUBSTITUTIONMODEL_H_
