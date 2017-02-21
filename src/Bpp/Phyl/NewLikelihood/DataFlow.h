//
// File: DataFlow.h
// Created by: François Gindraud
// Created on: mardi 21 février 2017, à 11h19
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

#ifndef _DATAFLOW_H_
#define _DATAFLOW_H_

#if !(__cplusplus >= 201103L)
#error "Bpp/Phyl/DataFlow.h requires C++11 support"
#endif

#include <cassert>
#include <utility> // move, forward

namespace bpp
{
  /**
   * @brief A class that stores a value that can be valid of invalid
   *
   * This class is temporary, as I will move to a more integrated data flow approach later.
   * It will be used to simplify code at first (grouping values and flags).
   */
  template <typename T>
  class CachedValue
  {
  private:
    T value_;
    bool valid_{false};

    void make_valid(void) { valid_ = true; }

  public:
    bool is_valid(void) const { return valid_; }
    void invalidate(void) { valid_ = false; }

    void set(const T& v)
    {
      value_ = v;
      make_valid();
    }
    void set(T&& v) { value_ = std::move(v); }

    const T& get(void) const
    {
      assert(is_valid());
      return value_;
    }

    // TODO deal with complex objects (release-like func that moves, invalid and access internal storage)
  };
} // end namespace bpp

#endif // _DATAFLOW_H_
