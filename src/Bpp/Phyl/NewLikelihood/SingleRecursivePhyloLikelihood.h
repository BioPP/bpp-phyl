//
// File: SingleRecursivePhyloLikelihood.h
// Created by: Julien Dutheil
// Created on: Tue Jul 23 16:01 2013
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

#ifndef _SINGLERECURSIVEPHYLOLIKELIHOOD_H_
#define _SINGLERECURSIVEPHYLOLIKELIHOOD_H_

#include "SinglePhyloLikelihood.h"

namespace bpp
{
namespace newlik
{

/**
 * @brief The SingleRecursivePhyloLikelihood class: phylogenetic likelihood computation with a single process and a single recursion.
 *
 * This class is a convenient wrapper for the more general SinglePhyloLikelihood class.
 */
class SingleRecursivePhyloLikelihood :
  public virtual SinglePhyloLikelihood
{
public:
  SingleRecursivePhyloLikelihood(SubstitutionProcess* process, bool verbose = true, bool usePatterns = true):
    SinglePhyloLikelihood(process, new SingleRecursiveTreeLikelihoodCalculation(process, verbose, usePatterns), verbose)
  {}

  SingleRecursivePhyloLikelihood(const SiteContainer& sites, SubstitutionProcess* process, bool verbose = true, bool usePatterns = true):
    SinglePhyloLikelihood(process, new SingleRecursiveTreeLikelihoodCalculation(sites, process, verbose, usePatterns), verbose)
  {}

  ~SingleRecursivePhyloLikelihood() {}

  SingleRecursivePhyloLikelihood* clone() const { return new SingleRecursivePhyloLikelihood(*this); }
  
};
} // end of namespace newlik.
} // end of namespace bpp.
  
#endif //_SINGLERECURSIVEPHYLOLIKELIHOOD_H_

