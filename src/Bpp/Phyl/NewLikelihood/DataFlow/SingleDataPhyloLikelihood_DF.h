//
// File: SingleDataPhyloLikelihood_DF.h
// Authors: Laurent Guéguen (2019)
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

#ifndef SINGLE_DATA_PHYLOLIKELIHOOD_DF_H
#define SINGLE_DATA_PHYLOLIKELIHOOD_DF_H

#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/SingleDataPhyloLikelihood.h>

#include "AlignedPhyloLikelihood_DF.h"

namespace bpp {

    /* Wraps a dataflow graph as a function: resultNode = f(variableNodes).
     *
     */
  
    class SingleDataPhyloLikelihood_DF :
      public SingleDataPhyloLikelihood,
      virtual public AlignedPhyloLikelihood_DF
    {
    protected:
      size_t nbStates_;

      /**
       * @brief Number of the concerned data.
       *
       **/
    
      size_t nData_;
    
    public:
      SingleDataPhyloLikelihood_DF(size_t nbSites, size_t nbStates, size_t nData = 0) :
        AlignedPhyloLikelihood_DF(nbSites),
        nbStates_(nbStates),
        nData_(nData)
      {}
    
    
      SingleDataPhyloLikelihood_DF(const SingleDataPhyloLikelihood_DF& asd) :
        AlignedPhyloLikelihood_DF(asd),
        nbStates_(asd.nbStates_),
        nData_(asd.nData_)
      {
      }
    
      virtual ~SingleDataPhyloLikelihood_DF() {}
    
      SingleDataPhyloLikelihood_DF* clone() const = 0;
    
      SingleDataPhyloLikelihood_DF& operator=(const SingleDataPhyloLikelihood_DF& asd)
      {
        AlignedPhyloLikelihood_DF::operator=(asd);
        nbStates_=asd.nbStates_;
      
        nData_=asd.nData_;
      
        return *this;
      }
    
      virtual void setData(const AlignedValuesContainer& sites, size_t nData = 0)
      {
        setNumberOfSites(sites.getNumberOfSites());
        nbStates_ = sites.getAlphabet()->getSize();
        nData_=nData;
      }
    
      size_t getNData() const
      {
        return nData_;
      }
    
      void setNData(size_t nData)
      {
        nData_=nData;
      }
    
      size_t getNumberOfStates() const { return nbStates_; }
    
    };
} // namespace bpp

#endif // SINGLE_DATA_PHYLOLIKELIHOOD_DF_H
