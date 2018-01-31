//
// File: YNGP_M8.h
// Created by: Laurent Gueguen
// Created on: May 2010
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

#ifndef _YNGP_M8_H_
#define _YNGP_M8_H_

#include "YNGP_M.h"

#include <Bpp/Seq/GeneticCode/GeneticCode.h>

namespace bpp
{

/**
 * @brief The Yang et al (2000) M8 substitution model for codons.
 * @author Laurent Guéguen
 *
 * This model is a mixture of models as described in YN98 class, the
 * mixture being defined on the selection parameter oomega to allow it
 * to vary among sites, following a mixture of a Beta distribution and
 * of another value above 1.
 *
 * This model includes 5 parameters (@f$\kappa@f$, @f$ p @f$ and
 * @f$q@f$ of the @f$ Beta(p,q) @f$ distribution, @f$p0@f$ the weight of the
 * Beta distribution and @f$\omega @f$ the selection parameter above 1
 * (with weight @f$ 1-p0 @f$)). The codon frequencies @f$ \pi_j @f$ are
 * either observed or infered.
 *
 * References:
 *
 * Yang, Z., R. Nielsen, N. Goldman, and A.-M. K. Pedersen (2000)
 * Genetics 155:431-449.
 * 
 */
  class YNGP_M8:
    public YNGP_M
  {
  public:
    /*
     *@brief Constructor that requires the number of classes of the
     * BetaDiscreteDistribution.
     *
     */
  
    YNGP_M8(const GeneticCode* gc, FrequenciesSet* codonFreqs, unsigned int nbclass);

    YNGP_M8* clone() const { return new YNGP_M8(*this); }

    YNGP_M8(const YNGP_M8& mod2) :
      YNGP_M(mod2)
    {
    }
  
    YNGP_M8& operator=(const YNGP_M8& mod2)
    {
      YNGP_M::operator=(mod2);
      return *this;
    }
  
  protected:
    void updateMatrices();

  public:
    std::string getName() const { return "YNGP_M8"; }

  };

} //end of namespace bpp.

#endif	//_YNGP_M8_H_

