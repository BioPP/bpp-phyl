//
// File: SiteSimulator.h
// Created by: Julien Dutheil
// created on: wed aug  24 15:20 2005
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

#ifndef _SITESIMULATOR_H_
#define _SITESIMULATOR_H_

// From bpp-seq:
#include <Bpp/Seq/Site.h>

namespace bpp
{

/**
 * @brief The SiteSimulator interface.
 * SiteSimulator classes can simulate single sites.
 *
 * @see SequenceSimulator interface for simulating whole sequence sets.
 */
class SiteSimulator
{
  public:
    SiteSimulator() {}
    virtual ~SiteSimulator() {}
    
  public:
    virtual Site* simulateSite() const = 0;
    virtual Site* simulateSite(size_t ancestralStateIndex) const = 0;
    virtual Site* simulateSite(size_t ancestralStateIndex, double rate) const = 0;
    virtual Site* simulateSite(double rate) const = 0;
    virtual std::vector<std::string> getSequencesNames() const = 0;
    virtual const Alphabet* getAlphabet() const = 0;
};

} //end of namespace bpp.

#endif // _SITESIMULATOR_H_

