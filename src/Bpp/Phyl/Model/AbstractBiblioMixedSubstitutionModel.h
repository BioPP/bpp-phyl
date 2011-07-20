//
// File: AbstractBiblioMixedSubstitutionModel.h
// Created by: Laurent Gueguen
// Created on: lundi 18 juillet 2011, à 15h 17
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

#ifndef _ABSTRACTBIBLIOMIXEDSUBSTITUTIONMODEL_H_
#define _ABSTRACTBIBLIOMIXEDSUBSTITUTIONMODEL_H_

#include "AbstractBiblioSubstitutionModel.h"
#include "MixedSubstitutionModel.h"

namespace bpp
{

/**
 * @brief Abstract class for mixture models based on the bibliography.
 * @author Laurent Guéguen
 *
 **/
  
  class AbstractBiblioMixedSubstitutionModel:
    public AbstractBiblioSubstitutionModel
  {
  public:
    AbstractBiblioMixedSubstitutionModel(const std::string& prefix);

    AbstractBiblioMixedSubstitutionModel(const AbstractBiblioMixedSubstitutionModel& model);

    AbstractBiblioMixedSubstitutionModel& operator=(const AbstractBiblioMixedSubstitutionModel& model);

    virtual ~AbstractBiblioMixedSubstitutionModel();

#ifndef NO_VIRTUAL_COV
    virtual AbstractBiblioMixedSubstitutionModel* clone() const = 0;
#endif

  public:
    virtual const MixedSubstitutionModel* getMixedModel() const = 0; 

    virtual MixedSubstitutionModel* getMixedModel() = 0;

    /*
     *@brief Returns the submodel from the mixture.
     *
     */
    
    const SubstitutionModel* getNModel(unsigned int i) const {
      return getMixedModel()->getNModel(i);
    }

    SubstitutionModel* getNModel(unsigned int i) {
      return getMixedModel()->getNModel(i);
    }

    /**
     * @brief Returns the  probability of a specific model from the mixture
     */
    
    double getNProbability(unsigned int i) const {
      return getMixedModel()->getNProbability(i);
    }

    /**
     * @brief Returns the vector of the probabilities of the
     * submodels of the mixture.
     *
     */
    
    const std::vector<double>& getProbabilities() const {
      return getMixedModel()->getProbabilities();
    }
    
    /**
     * @brief Sets the probabilities of the submodels of the mixture.
     *
     */
    
    void setNProbability(unsigned int i, double prob)
    {
      getMixedModel()->setNProbability(i, prob);
    }

    /**
     * @brief Returns the number of submodels
     *
     */
    
    unsigned int getNumberOfModels() const {
      return getMixedModel()->getNumberOfModels();
    }
  
    /**
     * @brief inactivated method to prevent out of model manipulations
     *
     **/
    
    //    virtual void setVRates(Vdouble & vd){};

  };
  
} //end of namespace bpp.

#endif	//_AbstractBiblioMixedSubstitutionModel_H_

