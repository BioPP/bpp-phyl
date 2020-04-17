//
// File: LikelihoodCalculationMultiProcess.h
// Authors: Laurent Guéguen (2018)
// Created: vendredi 3 avril 2020, à 20h 58
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.

  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".

  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.

  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef LIKELIHOOD_CALCULATION_MULTI_PROCESS_H
#define LIKELIHOOD_CALCULATION_MULTI_PROCESS_H

#include <Bpp/Numeric/AbstractParametrizable.h>

namespace bpp {

  /*
   * Build LikelihoodCalculation from a set of several
   * LikelihoodCalculation.
   *
   * The computation (ie passage from the set of
   * LikelihoodCalculation to this) will be defined in ad
   * hoc classes (specifically PhyloLikelihood classes).
   *
   *
   */

  template<class LikElem>
  class LikelihoodCalculationMultiProcess :
    public AbstractParametrizable
  {
    using LikRef = std::shared_ptr<LikElem>;

  private:
    std::vector<LikRef> vLikCal_;

  public:
    LikelihoodCalculationMultiProcess() :
      AbstractParametrizable("")
    {}

    LikelihoodCalculationMultiProcess(const LikelihoodCalculationMultiProcess<LikElem>& lik):
      AbstractParametrizable(""),
      vLikCal_(lik.vLikCal_)
    {
    }
    
    LikelihoodCalculationMultiProcess<LikElem>* clone() const
    {
      return new LikelihoodCalculationMultiProcess<LikElem>(*this);
    }

    LikRef getSingleLikelihood(size_t nL)
    {
      if (nL>vLikCal_.size())
        throw Exception("LikelihoodCalculationMultiProcess::getSingleLikelihood : Wrong number " + std::to_string(nL) + "not under " + std::to_string(vLikCal_.size()));
      return vLikCal_[nL];
    }

    void addSingleLikelihood(LikRef lik)
    {
      vLikCal_.push_back(lik);
      shareParameters_(lik->getParameters());
    }

    size_t getNumberOfSingleProcess() const
    {
      return vLikCal_.size();
    }

    bool isInitialized() const {

      for (auto& lik : vLikCal_)
        if (!lik->isInitialized())
          return false;
      
      return true;
    };

    /**************************************************/

    /*
     * @brief Return the loglikehood (see as a function, ie
     * computation node).
     *
     */
      
    void makeLikelihoods()
    {
      for (auto& lik : vLikCal_)
        lik->makeLikelihoods();
    }
  };

} // namespace bpp

#endif // LIKELIHOOD_CALCULATION_SINGLE_PROCESS_H

