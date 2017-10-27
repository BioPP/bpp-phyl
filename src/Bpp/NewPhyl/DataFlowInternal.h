//
// File: DataFlowInternal.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-10-27
// Last modified: 2017-10-27
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

#ifndef BPP_NEWPHYL_DATAFLOWINTERNAL_H
#define BPP_NEWPHYL_DATAFLOWINTERNAL_H

/* Internal useful primitives for writing DataFlow code.
 *
 * Contains the simplest primitives:
 * - raw value access
 */

#include <Bpp/NewPhyl/DataFlow.h>

namespace bpp {
namespace DF {
	/// Dynamic test that node inherits from Value<T>.
	template <typename T> bool isValueNode (const Node & n) noexcept {
		return dynamic_cast<const Value<T> *> (&n) != nullptr;
	}

	/// Unsafe cast from Node to Value<T> (faster than dynamic cast).
	template <typename T> const Value<T> & nodeValueCast (const Node & node) noexcept {
		assert (isValueNode<T> (node));
		return static_cast<const Value<T> &> (node);
	}

	/** Class used to access protected value_ in Value<T>.
	 * Should only be used for DataFlow internal code.
	 * Friend of node classes, bypasses encapsulation.
	 * Used for wrappers, for manipulating possibly invalid values.
	 */
	class InternalAccessor {
	public:
		template <typename T> static const T & value (const Value<T> & node) noexcept {
			return node.value_;
		}
		template <typename T> static T & value (Value<T> & node) noexcept { return node.value_; }
	};

	/// Access maybe invalid value as const.
	template <typename T> const T & accessValueConst (const Value<T> & node) noexcept {
		return InternalAccessor::value (node);
	}

	/// Access maybe invalid value as mutable.
	template <typename T> T & accessValueMutable (Value<T> & node) noexcept {
		return InternalAccessor::value (node);
	}

	/** Access maybe invalid const from raw Node &.
	 * Typically used to get matrix dimensions during DF graph transformations.
	 * Dimensions are preset in DFNumeric, so they are always valid.
	 */
	template <typename T> const T & accessValueConstCast (const Node & node) {
		return accessValueConst (nodeValueCast<T> (node));
	}

	/// Access const with assert checking validity.
	template <typename T> const T & accessValidValueConst (const Value<T> & node) {
		assert (node.isValid ());
		return accessValueConst (node);
	}

	/// Access valid const while casting from a raw Node &.
	template <typename T> const T & accessValidValueConstCast (const Node & node) {
		return accessValidValueConst (nodeValueCast<T> (node));
	}

	/// Access valid const while casting from NodeRef &.
	template <typename T> const T & accessValidValueConstCast (const NodeRef & nodeRef) {
		assert (nodeRef);
		return accessValidValueConstCast<T> (*nodeRef);
	}
} // namespace DF
} // namespace bpp

#endif // BPP_NEWPHYL_DATAFLOWINTERNAL_H
