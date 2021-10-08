//
// File: RELAX.h
// Created by: Keren Halabi
// Created on: May 2018
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

#ifndef _RELAX_H_
#define _RELAX_H_


#include "YNGP_M.h"
#include <Bpp/Seq/GeneticCode/GeneticCode.h>

namespace bpp
{
/**
 * @brief The RELAX (2014) branch-site model for codons
 *
 * @author Keren Halabi
 *
 * This model consists of a mixture of YN98 models, and allows for
 * mixture in two levels: the site level and the branch level.
 *
 * @f$\omega_0 = omega_1 * p < 1 @f$ (with probability @f$p_0 @f$)
 *
 * @f$\omega_1 <= 1 @f$ (with probability @f$p_1 @f$)
 *
 * @f$\omega_2 > 1 @f$ (with probability @f$1-p_1-p_0 @f$)
 *
 * Each omega is raised to the power of a selection intensity
 * parameter k.
 *
 * Used alone, this model is over-parameterized, and should be used in
 * a branch heterogeneous modeling with shared parameters with other
 * RELAX models.
 *
 * References:
 *
 * Joel O., et al. "RELAX: detecting relaxed selection in a phylogenetic framework."
 * Molecular biology and evolution 32.3 (2014): 820-832.‏
 */
class RELAX :
  public YNGP_M
{
public:
  // gc is a pointer to a constant GeneticCode instance
  // unlike GeneticCode const* gc - gc is a constant pointer to a GeneticCode instance
  // see: https://stackoverflow.com/questions/3984989/what-is-the-differnce-between-const-x-a-and-x-const-a-if-x-is-the-class?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
  RELAX(const GeneticCode* gc, std::shared_ptr<FrequencySet> codonFreqs);

  RELAX* clone() const { return new RELAX(*this); }

  RELAX(const RELAX& mod2) :
    YNGP_M(mod2)
  {}

  RELAX& operator=(const RELAX& mod2)
  {
    YNGP_M::operator=(mod2);
    return *this;
  }

protected:
  void updateMatrices();

public:
  std::string getName() const { return "RELAX"; }
};
} // end of namespace bpp.

#endif// _RELAX_H
