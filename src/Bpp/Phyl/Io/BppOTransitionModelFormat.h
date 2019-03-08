//
// File: BppOTransitionModelFormat.h
// Created by: Laurent Guéguen
// Created on: mercredi 4 juillet 2012, à 13h 26
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

#ifndef _BPPO_TRANSITION_MODELFORMAT_H_
#define _BPPO_TRANSITION_MODELFORMAT_H_

#include "BppOSubstitutionModelFormat.h"
#include "../Model/MixedTransitionModel.h"

namespace bpp
{
/**
 * @brief Transition model I/O in BppO format.
 *
 * Creates a new transition model object according to model description syntax
 * (see the Bio++ Progam Suite manual for a detailed description of this syntax).
 *
 */

  class BppOTransitionModelFormat :
    public BppOSubstitutionModelFormat
  {
  private:
    MixedTransitionModel* readMixed_(const Alphabet* alphabet, const std::string& modelDescription, const AlignedValuesContainer* data);

  public:
    /**
     * @brief Create a new BppOTransitionModelFormat object.
     *
     * @param alphabetCode     Bit saying which alphabets are allowed in the model specification.
     * @param allowCovarions   Tell is a covarion model can be returned.
     * @param allowMixed       Tell is a mixture model can be returned.
     * @param allowGaps        Tell is a gap model can be returned.
     * @param verbose          Tell if the construction is verbose.
     * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
     */
    BppOTransitionModelFormat(unsigned char alphabetCode, bool allowCovarions, bool allowMixed, bool allowGaps, bool verbose, int warn):
      BppOSubstitutionModelFormat(alphabetCode,allowCovarions,allowMixed,allowGaps,verbose,warn)
    {}

    BppOTransitionModelFormat(const BppOTransitionModelFormat& format):
      BppOSubstitutionModelFormat(format)
    {}

    BppOTransitionModelFormat& operator=(const BppOTransitionModelFormat& format)
    {
      BppOSubstitutionModelFormat::operator=(format);
      return *this;
    }

    virtual ~BppOTransitionModelFormat() {}

  public:
    TransitionModel* readTransitionModel(const Alphabet* alphabet, const std::string& modelDescription, const AlignedValuesContainer* data = 0, bool parseArguments = true);

  };

} // end of namespace bpp.

#endif // _BPPOSUBSTITUTIONMODELFORMAT_H_

