//
// File: Utils.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-07-04 00:00:00
// Last modified: 2018-05-24
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

#ifndef BPP_NEWPHYL_UTILS_H
#define BPP_NEWPHYL_UTILS_H



// FIXME move away ?

#include <functional> // std::hash
#include <utility>
#include <vector>

namespace bpp {
/** @brief Create a new vector filled with results from calling a function on another vector.
 *
 * The type of the new vector is std::vector<T> for T the result type of the function.
 */
template <typename Container, typename Function>
auto mapToVector (const Container & container, Function function)
    -> std::vector<decltype (function (container.front ()))> {
	std::vector<decltype (function (container.front ()))> r;
	r.reserve (container.size ());
	for (const auto & v : container)
		r.emplace_back (function (v));
	return r;
}

/** @brief hash capability for std::vector
 *
 * std::Vector<T> can be used as a key for hash tables (std::unordered_map / set).
 * This requires the map/set to use this custom hasher.
 *
 * The vector itself must be const after insertion in the table.
 * If not, the key will change, which breaks the hash table invariants.
 */
template <typename T> struct VectorHasher {
	size_t operator() (const std::vector<T> & vec) const {
		auto h = vec.size ();
		std::hash<T> hasher{};
		for (const auto & e : vec)
			h ^= hasher (e) + 0x9e3779b9 + (h << 6) + (h >> 2);
		return h;
	}
};
} // namespace bpp
#endif // BPP_NEWPHYL_UTILS_H
