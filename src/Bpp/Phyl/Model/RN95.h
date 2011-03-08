//
// File: RN95.h
// Created by: Laurent Guéguen
// Created on: jeudi 24 février 2011, à 20h 43
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

#ifndef _RN95_H_
#define _RN95_H_

#include "NucleotideSubstitutionModel.h"
#include "AbstractSubstitutionModel.h"

#include <Bpp/Numeric/Constraints.h>

// From SeqLib:
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>

namespace bpp
{

  /**
   * @brief The model described by Rhetsky \& Ney, where the only
   * hypothesis is that the transversion rates are only dependent of
   * the target nucleotide. This model is not reversible.
   *
   * After normalization (not analytic) this model has 7 parameters:
   * \f[
   * Q= \frac 1K
   * \begin{pmatrix}
   * . & \beta_3 & \alpha_4 & \beta_2 \\
   * \beta_1 & . &  1 & \alpha_2 \\
   * \alpha_1 & \beta_3 & . & \beta_2 \\
   * \beta_1 & \alpha_3 & 1 & .\\
   * \end{pmatrix}
   *\f]
   *
   * The generator of this model is diagonalized numerically.
   * See AbstractSubstitutionModel for details of how the probabilities are computed.
   *
   * The parameters are named \c "alpha1", \c "alpha2", \c "alpha3",
   * \c "alpha4", \c "beta1", \c "beta2", \c "beta3". Their values are
   * positive.
   *
   * Reference:
   * - Rhetsky A. \& Ney M. (1995) MBE 12(1) 131-151.
   */
  
  class RN95:
    public virtual NucleotideSubstitutionModel,
    public AbstractSubstitutionModel
  {
  private:
    double alpha1_, alpha2_, alpha3_, alpha4_, beta1_, beta2_, beta3_;
  
  public:
    RN95(
         const NucleicAlphabet* alphabet,
         double alpha1 = 1,
         double alpha2 = 1,
         double alpha3 = 1,
         double alpha4 = 1,
         double beta1 = 1,
         double beta2 = 1,
         double beta3 = 1);
  
    virtual ~RN95() {}
  
#ifndef NO_VIRTUAL_COV
    RN95*
#else
    Clonable*
#endif
    clone() const { return new RN95(*this); }
  
  public:
    std::string getName() const { return "RN95"; }
  
    void updateMatrices();
  
    /**
     * @brief This method is undefined in this model
     */
    void setFreq(std::map<int, double>&);
  };

} //end of namespace bpp.

#endif	//_RN95_H_

