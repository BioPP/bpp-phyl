//
// File: YNGP_M3.h
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

#ifndef _YNGP_M3_H_
#define _YNGP_M3_H_

#include "YNGP_M.h"
#include <Bpp/Seq/GeneticCode/GeneticCode.h>

namespace bpp
{

/**
 * @brief The Yang et al (2000) M3 substitution model for codons.
 * @author Laurent Guéguen
 *
 * This model is a mixture of models as described in YN98 class, the
 * mixture being defined on the selection parameter to allow it to
 * vary among sites. There are $K$ selection parameters @f$ \omega_0 <
 * ... \omega_{K-1} @f$, with their respective probabilities @f$ p_0,
 * ..., p_{K-1} @f$ with @f$ p_0+p_1+...+p_{K-1}=1@f$. To garantee
 * that the @f$\omega_i@f$ are in increasing order, we define
 * @f$\delta_i=\omega_i - \omega_{i-1}@f$.
 *
 * This model includes 2*K parameters (@f$\kappa@f$, relative
 * probabilities @f$ theta1, theta2, ..., thetaK-1 @f$ and @f$omega0,
 * delta1, deltaK-1@f$). The codon frequencies @f$\pi_j@f$ are either
 * observed or infered.
 *
 * References:
 *
 * Yang, Z., R. Nielsen, N. Goldman, and A.-M. K. Pedersen (2000)
 * Genetics 155:431-449.
 * 
 */
  class YNGP_M3:
    public YNGP_M
  {
  public:
    YNGP_M3(const GeneticCode* gc, std::shared_ptr<FrequencySet> codonFreqs, unsigned int nclass = 3);

    YNGP_M3* clone() const { return new YNGP_M3(*this); }

    YNGP_M3(const YNGP_M3& mod2) :
      YNGP_M(mod2)
    {
    }

    YNGP_M3& operator=(const YNGP_M3& mod2)
    {
      YNGP_M::operator=(mod2);
      return *this;
    }

  protected:
    void updateMatrices();

  public:
    std::string getName() const { return "YNGP_M3"; }

  };

} //end of namespace bpp.

#endif	//_YNGP_M3_H_

