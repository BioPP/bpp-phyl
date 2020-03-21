//
// File: AlignedPhyloLikelihood_DF.h
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

#ifndef ALIGNED_PHYLOLIKELIHOOD_DF_H
#define ALIGNED_PHYLOLIKELIHOOD_DF_H

#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/AlignedPhyloLikelihood.h>

namespace bpp {

  namespace dataflow
  {
    
    /* Wraps a dataflow graph as a function: resultNode = f(variableNodes).
     *
     */
  
    class AlignedPhyloLikelihood_DF :
      virtual public AlignedPhyloLikelihood
    {
    protected:

      size_t nbSites_;

    public:      
      AlignedPhyloLikelihood_DF (size_t nbSites) :
        nbSites_(nbSites)
      {}

      AlignedPhyloLikelihood_DF(const AlignedPhyloLikelihood_DF& asd) :
        nbSites_(asd.nbSites_)
      {}
      
      AlignedPhyloLikelihood_DF& operator=(const AlignedPhyloLikelihood_DF& asd)
      {
        nbSites_=asd.nbSites_;        
        return *this;
      }

      size_t getNumberOfSites() const { return nbSites_; }

    protected:
      void setNumberOfSites(size_t nbSites)
      {
        nbSites_ = nbSites;
      }
      
    };
  }
} // namespace bpp

#endif // ALIGNEDPHYLOLIKELIHOOD_DF_H
