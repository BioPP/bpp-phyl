//
// File: Debug.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-04-28
// Last modified: 2017-04-28
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

#pragma once
#ifndef BPP_NEWPHYL_DEBUG_H
#define BPP_NEWPHYL_DEBUG_H

#include <Bpp/NewPhyl/FrozenPtr.h>
#include <iosfwd>
#include <string>
#include <typeindex>
#include <typeinfo>

namespace bpp {
// Demangle a C++ symbol name
std::string demangle (const char * name);

// Pretty type name
std::string prettyTypeName (const std::type_info & ti);
std::string prettyTypeName (std::type_index ti);
template <typename T> std::string prettyTypeName () {
	return prettyTypeName (typeid (T));
}
template <typename T> std::string prettyTypeName (const T & t) {
	return prettyTypeName (typeid (t));
}

// TODO more general to_string ? (duck::format ?)
inline const std::string & debug_to_string (const std::string & s) {
	return s;
}
inline std::string && debug_to_string (std::string && s) {
	return std::move (s);
}
template <typename T, typename = decltype (std::to_string (std::declval<T> ()))>
std::string debug_to_string (T && t) {
	return std::to_string (std::forward<T> (t));
}

namespace Topology {
	class Tree;

	// Output a dot format graph representing the tree
	void debugTree (std::ostream & os, FrozenPtr<Tree> tree);
}
namespace DF {
	class Node;
	class NodeSpecification;
	class Registry;

	// Output a dot format graph representing the dataflow dag
	void debugDag (std::ostream & os, const Node & entryPoint);

	// Outputs debugDag + registry pointers to dataflow dag
	void debugRegistry (std::ostream & os, const Registry & registry);

	// Instantiate a nodespec, prints nodespec hierarchy & DF graph
	void debugNodeSpecInstantiation (std::ostream & os, const NodeSpecification & nodeSpec);

	// Trace previous instantiation of a NodeSpec in a registry
	// Prints: registry, nodespecs, df graph
	void debugNodeSpecInstantiationInRegistry (std::ostream & os, const NodeSpecification & nodeSpec,
	                                           const Registry & registry,
	                                           bool showRegistryLinks = false);
}
}

#endif // BPP_NEWPHYL_DEBUG_H
