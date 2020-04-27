//
// File IOFrequencySetFactory.h
// Created by: Laurent Guéguen
// Created on: lundi 9 juillet 2012, à 12h 58
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for sequences analysis.

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

#ifndef _IOFREQUENCYSETFACTORY_H_
#define _IOFREQUENCYSETFACTORY_H_

#include "../Model/SubstitutionModel.h"
#include "IoFrequencySet.h"
#include <Bpp/Exceptions.h>

//From the STL:
#include <string>

namespace bpp
{

/**
 * @brief Utilitary class for creating frequencies set readers and writers.
 *
 * @see IOSequenceFactory
 * @see IOTreeFactory
 */
class IOFrequencySetFactory
{
public:
  static const std::string BPPO_FORMAT;  

public:

  /**
   * @brief Creates a new factory object.
   *
   * Example:
   * @code
   * IFrequencySet * freqReader = IOFrequencySetFactory().createReader(IOFrequencySetFactory::BPP_FORMAT);
   * FrequencySet * freqset = freqReader->read(...);
   * delete freqReader;
   * @endcode
   */
  IOFrequencySetFactory() {}
  virtual ~IOFrequencySetFactory() {}
  
  /**
   * @brief Get a new dynamically created IFrequencySet object.
   *
   * @param format The input file format.
   * @return A pointer toward a new IFrequencySet object.
   * @throw Exception If the format name do not match any available format.
   */
  virtual IFrequencySet* createReader(const std::string& format);
  
  /**
   * @brief Get a new dynamically created OFrequencySet object.
   *
   * @param format The output file format.
   * @return A pointer toward a new OFrequencySet object.
   * @throw Exception If the format name do not match any available format.
   */
  virtual OFrequencySet * createWriter(const std::string& format);
};

} //end of namespace bpp.

#endif //_IOFREQUENCYSETFACTORY_H_

