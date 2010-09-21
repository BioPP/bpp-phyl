//
// File: YNGKP_M2.h
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

#ifndef _YNGKP_M2_H_
#define _YNGKP_M2_H_

#include "MixedSubstitutionModel.h"
#include "MixtureOfSubstitutionModels.h"
#include "FrequenciesSet.h"

#include <Bpp/Seq/GeneticCode/GeneticCode.h>

namespace bpp
{

/**
 * @brief The Yang et al (2000) M2 substitution model for codons, with
 * the more realistic modification in Wong & al (2004).
 * @author Laurent Guéguen
 *
 * This model is a mixture of models as described in YN98 class, the
 * mixture being defined on the selection parameter to allow it to
 * vary among sites. A site is either negatively selected @f$ 0 <
 * \omega_0 < 1 @f$ (with probability @f$p_0 @f$), or neutral (@f$
 * \omega_1 = 1 @f$) with probability @f$p_1 @f$, or positively
 * selected @f$ 1 < \omega_2 @f$ (with probability @f$p_2 @f$). This
 * model includes 5 parameters (@f$kappa@f$, @f$ theta1=p0,
 * theta2=\frac{p1}{p1+p2}, omega0 @f$ and @f$ omega2 @f$). The codon
 * frequencies @f$\pi_j@f$ are either observed or infered.
 *
 * References:
 *
 * Yang, Z., R. Nielsen, N. Goldman, and A.-M. K. Pedersen (2000)
 * Genetics 155:431-449.
 * 
 * Wong, W. S. W., Z. Yang, N. Goldman, and R. Nielsen. (2004)
 * Genetics 168:1041--1051.
 */
class YNGKP_M2:
  public MixedSubstitutionModel
{
private:
  MixtureOfSubstitutionModels* pmixmodel_;

private:
  /**
   * @brief Tools to make the link between the Parameters of the
   * object and those of pmixmodel_.
   *
   */

  std::map<std::string,std::string> mapParNamesFromPmodel_;
  
  ParameterList lParPmodel_;
  

public:
  YNGKP_M2(const GeneticCode* gc, FrequenciesSet* codonFreqs);

  ~YNGKP_M2();
  
  YNGKP_M2* clone() const { return new YNGKP_M2(*this); }

  YNGKP_M2(const YNGKP_M2&);

  YNGKP_M2& operator=(const YNGKP_M2&);

  void setFreq(std::map<int, double>& m) { pmixmodel_->setFreq(m);  }

  unsigned int getNumberOfStates() const  { return pmixmodel_->getNumberOfStates();  }

protected:
  void updateMatrices();

public:
  const SubstitutionModel* getNModel(unsigned int i) const {
    return pmixmodel_->getNModel(i);
  }

  SubstitutionModel* getNModel(unsigned int i) {
    return pmixmodel_->getNModel(i);
  }

  /**
   * @brief Returns the  probability of a specific model from the mixture
   */
  
  double getNProbability(unsigned int i) const {
    return pmixmodel_->getNProbability(i);
  }

  const std::vector<double>& getProbabilities() const {
    return pmixmodel_->getProbabilities();
  }
    
  unsigned int getNumberOfModels() const {
    return pmixmodel_->getNumberOfModels();
  }

  std::string getName() const { return "YNGKP_M2"; }

  double Pij_t(unsigned int i, unsigned int j, double t) const {
    return pmixmodel_->Pij_t(i,j,t);
  }
  double dPij_dt(unsigned int i, unsigned int j, double t) const {
    return pmixmodel_->dPij_dt(i,j,t);
  };
  double d2Pij_dt2(unsigned int i, unsigned int j, double t) const {
    return pmixmodel_->dPij_dt(i,j,t);
  };
  const Matrix<double>& getPij_t(double t) const {
    return pmixmodel_->getPij_t(t);
  };
  const Matrix<double>& getdPij_dt(double t) const {
    return pmixmodel_->getdPij_dt(t);
  };
  const Matrix<double>& getd2Pij_dt2(double t) const {
    return pmixmodel_->getd2Pij_dt2(t);
  };
  
  const Vdouble& getFrequencies() {
    return pmixmodel_->getFrequencies();
  };
  
  double freq(unsigned int i) const {
    return pmixmodel_->freq(i);
  };


};

} //end of namespace bpp.

#endif	//_YNGKP_M2_H_

