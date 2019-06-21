//
// File: YNGP_M.h
// Created by: Laurent Gueguen
// Created on: mardi 26 septembre 2017, à 23h 06
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

#ifndef _YNGP_M_H_
#define _YNGP_M_H_

#include "../AbstractBiblioMixedTransitionModel.h"
#include "../FrequenciesSet/CodonFrequenciesSet.h"
#include "../MixtureOfASubstitutionModel.h"

#include <Bpp/Seq/GeneticCode/GeneticCode.h>

namespace bpp
{

/**
 * @brief Abstract generic class for The Yang et al (2000) M
 * substitution models for codons. al (2004).
 * @author Laurent Guéguen
 *
 * References:
 *
 * Yang, Z., R. Nielsen, N. Goldman, and A.-M. K. Pedersen (2000)
 * Genetics 155:431-449.
 * 
 * Wong, W. S. W., Z. Yang, N. Goldman, and R. Nielsen. (2004)
 * Genetics 168:1041--1051.
 */

  class YNGP_M:
    public AbstractBiblioMixedTransitionModel,
    virtual public TransitionModel
  {
  protected:

    /*
     * redefined mixed model pointer
     *
     */
    
    const MixtureOfASubstitutionModel*  pmixsubmodel_;
  
    /**
     * @brief indexes of 2 codons states between which the substitution is
     * synonymous, to set a basis to the homogeneization of the rates.
     *
     */
    size_t synfrom_, synto_;
  
  public:
    YNGP_M(const std::string& name) :
      AbstractBiblioMixedTransitionModel(name),
      pmixsubmodel_(),
      synfrom_(),
      synto_()
    {
      pmixsubmodel_=dynamic_cast<const MixtureOfASubstitutionModel*>(&getMixedModel());      
    }
    
    YNGP_M(const YNGP_M& mod2) :
      AbstractBiblioMixedTransitionModel(mod2),
      pmixsubmodel_(),
      synfrom_(mod2.synfrom_),
      synto_(mod2.synto_)
    {
      pmixsubmodel_=dynamic_cast<const MixtureOfASubstitutionModel*>(&getMixedModel());      
    }

    virtual YNGP_M* clone() const = 0;
    
    YNGP_M& operator=(const YNGP_M& mod2)
    {
      const auto& eq = AbstractBiblioMixedTransitionModel::operator=(mod2);
      
      synfrom_ = mod2.synfrom_;
      synto_ = mod2.synto_;
      const auto& mm=eq.getMixedModel();
      
      pmixsubmodel_ = &dynamic_cast<const MixtureOfASubstitutionModel&>(mm);
        
      return *this;
    }
    
  };

} //end of namespace bpp.

#endif	//_YNGP_M_H_

