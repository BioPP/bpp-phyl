//
// File: DataFlow.cpp
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-05-30
// Last modified: 2017-05-30
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

#include <Bpp/Exceptions.h>
#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/DataFlowTemplateUtils.h>
#include <Bpp/NewPhyl/Debug.h>
#include <algorithm>
#include <stack>
#include <typeinfo>

namespace bpp {
namespace DF {
	// Error functions DataFlow.h
	void failureNodeConversion (const std::type_info & handleType, const Node & node) {
		throw Exception (prettyTypeName (handleType) +
		                 " cannot store: " + prettyTypeName (typeid (node)));
	}

	void failureComputeWasCalled (const std::type_info & nodeType) {
		throw Exception (prettyTypeName (nodeType) + ": compute() was called");
	}

	// Error functions DataFlowTemplateUtils.h
	static void failureDependencyNumberMismatch (const std::type_info & inNodeType,
	                                             SizeType expectedSize, SizeType givenSize) {
		throw Exception (prettyTypeName (inNodeType) + ": expected " + std::to_string (expectedSize) +
		                 " dependencies, got " + std::to_string (givenSize));
	}
	void checkDependencyNumber (const std::type_info & inNodeType, SizeType expectedSize,
	                            SizeType givenSize) {
		if (expectedSize != givenSize)
			failureDependencyNumberMismatch (inNodeType, expectedSize, givenSize);
	}

	void failureEmptyDependency (const std::type_info & inNodeType, IndexType depIndex) {
		throw Exception (prettyTypeName (inNodeType) + ": " + std::to_string (depIndex) +
		                 "-th dependency is empty (nullptr)");
	}

	void failureDependencyTypeMismatch (const std::type_info & inNodeType, IndexType depIndex,
	                                    const std::type_info & expectedType, const Node & givenNode) {
		throw Exception (prettyTypeName (inNodeType) + ": expected class derived from " +
		                 prettyTypeName (expectedType) + " as " + std::to_string (depIndex) +
		                 "-th dependency, got " + prettyTypeName (typeid (givenNode)));
	}

	// Node impls

	Node::Node (const NodeRefVec & dependencies) : dependencyNodes_ (dependencies) {
		for (auto & n : dependencyNodes_)
			n->registerNode (this);
	}
	Node::Node (NodeRefVec && dependencies) : dependencyNodes_ (std::move (dependencies)) {
		for (auto & n : dependencyNodes_)
			n->registerNode (this);
	}

	Node::~Node () {
		for (auto & n : dependencyNodes_)
			n->unregisterNode (this);
	}

	void Node::invalidate () noexcept {
		if (!isValid ())
			return;
		std::stack<Node *> nodesToInvalidate;
		nodesToInvalidate.push (this);
		while (!nodesToInvalidate.empty ()) {
			auto * n = nodesToInvalidate.top ();
			nodesToInvalidate.pop ();
			if (n->isValid ()) {
				n->isValid_ = false;
				for (auto * dependent : n->dependentNodes_)
					nodesToInvalidate.push (dependent);
			}
		}
	}

	void Node::computeRecursively () {
		// Compute the current node (and dependencies recursively) if needed
		if (isValid ())
			return;

		// Discover then recompute needed nodes
		std::stack<Node *> nodesToVisit;
		std::stack<Node *> nodesToRecompute;
		nodesToVisit.push (this);
		while (!nodesToVisit.empty ()) {
			auto * n = nodesToVisit.top ();
			nodesToVisit.pop ();
			if (!n->isValid ())
				nodesToRecompute.push (n);
			for (auto & dep : n->dependencies ())
				nodesToVisit.push (dep.get ());
		}
		while (!nodesToRecompute.empty ()) {
			auto * n = nodesToRecompute.top ();
			nodesToRecompute.pop ();
			n->compute ();
			n->makeValid ();
		}
	}

	NodeRef Node::derive (const Node & variable) {
		throw Exception ("Node does not support derivation: " + description ());
	}

	std::string Node::description () const { return prettyTypeName (typeid (*this)); }

	std::string Node::debugInfo () const { return {}; }

	bool Node::isConstant () const { return false; }

	void Node::appendDependency (NodeRef node) {
		node->registerNode (this);
		dependencyNodes_.emplace_back (std::move (node));
		invalidate ();
	}

	void Node::registerNode (Node * n) { dependentNodes_.emplace_back (n); }
	void Node::unregisterNode (const Node * n) {
		dependentNodes_.erase (std::remove (dependentNodes_.begin (), dependentNodes_.end (), n),
		                       dependentNodes_.end ());
	}

		// TODO use in TemplateUtils for merging
#if 0
#include <Bpp/NewPhyl/Optional.h>
#include <typeindex>
#include <unordered_map>
	class Registry {
	public:
		class Key {
			/* Nodes are indexed by node type and dependencies.
			 * This is sufficient to ensure merging of DF nodes with similar values.
			 *
			 * Key is a temporary class (stores a ref to the dependency vector).
			 * When searching, the ref should point to the built dep vector.
			 * When stored, the ref points to the Node dep vector.
			 * The Key will be destroyed at the same time as the Node.
			 */
		public:
			Key (std::type_index nodeType, const NodeRefVec & dependencies)
			    : nodeType_ (nodeType), dependencies_ (dependencies) {}

			bool operator== (const Key & other) const noexcept {
				return nodeType_ == other.nodeType_ && dependencies_ == other.dependencies_;
			}
			std::size_t hashCode () const noexcept {
				auto nodeTypeHash = std::hash<std::type_index>{}(nodeType_);
				auto vecHash = std::hash<NodeRefVec>{}(dependencies_);
				return vecHash ^ (nodeTypeHash << 1);
			}

			std::type_index operation () const noexcept { return nodeType_; }
			const NodeRefVec & dependencies () const noexcept { return dependencies_; }

		private:
			std::type_index nodeType_;
			const NodeRefVec & dependencies_;
		};

		Optional<NodeRef> get (const Key & key) const {
			return optional_find (nodes_, key).cast<NodeRef> ();
		}

		void set (NodeRef nodeRef) {
			// Node is independent of NodeRef moves, so we can build a key referencing its dependencies.
			auto & node = *nodeRef;
			// By definition of the standard, typeid (node) performs lookup to derived type.
			auto result = nodes_.emplace (Key{typeid (node), node.dependencies ()}, std::move (nodeRef));
			if (!result.second)
				throw std::runtime_error ("Registry::set: key already used");
		}

		template <typename Callable> void foreachKeyValue (Callable callable) const {
			for (auto & it : nodes_)
				callable (it.first, it.second);
		}

	private:
		struct Hash {
			std::size_t operator() (const Key & k) const noexcept { return k.hashCode (); }
		};
		std::unordered_map<Key, NodeRef, Hash> nodes_;
	};
#endif
} // namespace DF
} // namespace bpp
