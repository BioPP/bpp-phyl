// File: FixedVector.h
// Authors:
//   Francois Gindraud (2017)
// Created: 15/03/2017

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

#ifndef _FIXEDVECTOR_H_
#define _FIXEDVECTOR_H_

#include <memory>
#include <iterator>
#include <type_traits>

namespace bpp
{
  /** Non redimensionnable vector.
   *
   * Usable by non copyable non movable classes.
   * Copyable or movable classes can just use vector.
   */
  template <typename T>
  class FixedVector
  {
    private:
      std::unique_ptr<T[]> data_ {};
      std::size_t size_ {};

    public:
    // TODO will be used when optimizing...
  };
} // end namespace bpp
#endif // _FIXEDVECTOR_H_
