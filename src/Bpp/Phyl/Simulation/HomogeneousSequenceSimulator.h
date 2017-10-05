//
// File: HomogeneousSequenceSimulator.h
// Created by: Julien Dutheil
// Created on: Wed Aug  24 15:20 2005
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

#ifndef _HOMOGENEOUSSEQUENCESIMULATOR_H_
#define _HOMOGENEOUSSEQUENCESIMULATOR_H_

#include "NonHomogeneousSequenceSimulator.h"

namespace bpp
{

/**
 * @brief Site and sequences simulation under homogeneous models.
 *
 * This is an alias class, for clarity and backward compatibility.
 */
  class HomogeneousSequenceSimulator:
    public NonHomogeneousSequenceSimulator
  {
  public:
		
    HomogeneousSequenceSimulator(
      const TransitionModel* model,
      const DiscreteDistribution* rate,
      const Tree* tree
      ) : NonHomogeneousSequenceSimulator(model, rate, tree) {}
			
    virtual ~HomogeneousSequenceSimulator() {}

  public:
    const TransitionModel* getModel() const
    {
      return getSubstitutionModelSet()->getModel(0);
    }
	
  };

} //end of namespace bpp.

#endif //_HOMOGENEOUSSEQUENCESIMULATOR_H_

