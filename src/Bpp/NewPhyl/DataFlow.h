//
// File: DataFlow.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-04-18
// Last modified: 2017-04-18
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
#ifndef BPP_NEWPHYL_DATAFLOW_H
#define BPP_NEWPHYL_DATAFLOW_H

#include <algorithm>
#include <iosfwd>
#include <memory>
#include <typeinfo> // std::bad_cast
#include <utility>
#include <vector>

#include <iostream>

namespace bpp {
namespace DF {

	/* Base Node.
	 */
	class Node {
	public:
		class Impl;
		using Ref = Impl &;

		template <typename U, typename = typename U::Impl>
		explicit Node (const U & u) : pImpl_ (std::dynamic_pointer_cast<Impl> (u.getShared ())) {
			if (!pImpl_)
				throw std::bad_cast ();
		}
		const std::shared_ptr<Impl> & getShared () const noexcept { return pImpl_; }

		template <typename T, typename... Args> static Node create (Args &&... args) {
			return Node (std::make_shared<T> (std::forward<Args> (args)...));
		}

		Ref get () const noexcept { return *pImpl_; }

	private:
		explicit Node (std::shared_ptr<Impl> p) noexcept : pImpl_ (std::move (p)) {}
		std::shared_ptr<Impl> pImpl_;
	};

	class Node::Impl {
	public:
		Impl () = default;
		Impl (const Impl &) = delete;
		Impl (Impl &&) = delete;
		Impl & operator= (const Impl &) = delete;
		Impl & operator= (Impl &&) = delete;

		virtual ~Impl () {
			foreachDependencyNode ([this](Impl * node) { node->unregisterNode (this); });
		}

		bool isValid () const noexcept { return isValid_; }
		void invalidate () noexcept {
			if (isValid ()) {
				isValid_ = false;
				foreachDependentNode ([](Impl * node) { node->invalidate (); });
			}
		}

		void registerNode (Impl * n) { dependentNodes_.emplace_back (n); }
		void unregisterNode (const Impl * n) {
			dependentNodes_.erase (std::remove (dependentNodes_.begin (), dependentNodes_.end (), n),
			                       dependentNodes_.end ());
		}

		virtual void compute () = 0;

		// TODO Replace with ranges
		template <typename F> void foreachDependentNode (F f) const {
			for (auto & p : dependentNodes_)
				f (p);
		}
		template <typename F> void foreachDependencyNode (F f) const {
			for (auto & p : dependencyNodes_)
				f (&p.get ());
		}

	protected:
		void makeValid () noexcept { isValid_ = true; }

		// TODO definitely need a small opt vector
		std::vector<Impl *> dependentNodes_{}; // Nodes that depend on us.
		std::vector<Node> dependencyNodes_{};  // Nodes that we depend on.

	private:
		bool isValid_{false};
	};

	void debugDagStructure (std::ostream & os, const Node & entryPoint);
	void debugDag (std::ostream & os, const Node & entryPoint);

	/* Valued node.
	 */
	template <typename T> class Value {
	public:
		class Impl;
		using Ref = Impl &;

		template <typename U, typename = typename U::Impl>
		explicit Value (const U & u) : pImpl_ (std::dynamic_pointer_cast<Impl> (u.getShared ())) {
			if (!pImpl_)
				throw std::bad_cast ();
		}
		const std::shared_ptr<Impl> & getShared () const noexcept { return pImpl_; }

		template <typename U, typename... Args> static Value create (Args &&... args) {
			return Value (std::make_shared<U> (std::forward<Args> (args)...));
		}

		Ref get () const noexcept { return *pImpl_; }
		const T & getValue () noexcept { return pImpl_->getValue (); }

	private:
		explicit Value (std::shared_ptr<Impl> p) noexcept : pImpl_ (std::move (p)) {}
		std::shared_ptr<Impl> pImpl_;
	};

	template <typename T> class Value<T>::Impl : public Node::Impl {
	public:
		template <typename... Args>
		Impl (Args &&... args) : Node::Impl (), value_ (std::forward<Args> (args)...) {}

		const T & getValue () {
			if (!this->isValid ()) {
				this->compute ();
				this->makeValid ();
			}
			return value_;
		}

	protected:
		T value_;
	};

	/* Param node.
	 */
	template <typename T> class Parameter {
	public:
		class Impl;
		using Ref = Impl &;

		template <typename U, typename = typename U::Impl>
		explicit Parameter (const U & u) : pImpl_ (std::dynamic_pointer_cast<Impl> (u.getShared ())) {
			if (!pImpl_)
				throw std::bad_cast ();
		}
		const std::shared_ptr<Impl> & getShared () const noexcept { return pImpl_; }

		template <typename... Args> static Parameter create (Args &&... args) {
			return Parameter (std::make_shared<Impl> (std::forward<Args> (args)...));
		}

		Ref get () const noexcept { return *pImpl_; }
		const T & getValue () noexcept { return pImpl_->getValue (); }
		void setValue (T t) noexcept { pImpl_->setValue (std::move (t)); }

	private:
		explicit Parameter (std::shared_ptr<Impl> p) noexcept : pImpl_ (std::move (p)) {}
		std::shared_ptr<Impl> pImpl_;
	};

	template <typename T> class Parameter<T>::Impl : public Value<T>::Impl {
	public:
		template <typename... Args>
		Impl (Args &&... args) : Value<T>::Impl (std::forward<Args> (args)...) {
			this->makeValid ();
		}

		void setValue (T t) noexcept {
			this->invalidate ();
			this->value_ = std::move (t);
			this->makeValid ();
		}

	private:
		void compute () override final { throw std::runtime_error ("compute(): parameter node"); }
	};
}
}

#endif // BPP_NEWPHYL_DATAFLOW_H
