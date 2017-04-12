#include <algorithm>
#include <memory>
#include <vector>

#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace DF {
///////////////////////////////////////// DF nodes //////////////////////////////////////

class BaseNode {
	// Base node type
public:
	BaseNode () = default;
	BaseNode (const BaseNode &) = delete;
	BaseNode (BaseNode &&) = delete;
	BaseNode & operator= (const BaseNode &) = delete;
	BaseNode & operator= (BaseNode &&) = delete;
	virtual ~BaseNode () = default;

	bool is_valid () const { return is_valid_; }
	void invalidate () {
		if (is_valid_) {
			is_valid_ = false;
			for (auto p : dependent_nodes_)
				p->invalidate ();
		}
	}
	void register_node (BaseNode * n) { dependent_nodes_.emplace_back (n); }
	void unregister_node (const BaseNode * n) {
		dependent_nodes_.erase (std::remove (dependent_nodes_.begin (), dependent_nodes_.end (), n),
		                        dependent_nodes_.end ());
	}

	virtual void compute () = 0;

protected:
	bool is_valid_{false};
	std::vector<BaseNode *> dependent_nodes_{};
};

// Value type for any node
class Node {
public:
	Node (std::shared_ptr<BaseNode> p) : p_node_ (std::move (p)) {}

private:
	std::shared_ptr<BaseNode> p_node_;
};

template <typename T> struct ValuedNode : public BaseNode {
	// Hold a T value
public:
	// Decouple get_value from compute later
	const T & get_value () {
		if (!this->is_valid_) {
			compute ();
			this->is_valid_ = true;
		}
		return value_;
	}

protected:
	T value_;
};

// Value type for ValuedNode
template <typename T> class Value {
public:
	Value (std::shared_ptr<ValuedNode<T>> p) : p_node_ (std::move (p)) {}
	const T & get_value () const { return p_node_->get_value (); }

private:
	std::shared_ptr<ValuedNode<T>> p_node_;
};

template <typename T> struct ParamNode final : public ValuedNode<T> {
	void set_value (T t) {
		this->value_ = t;
		this->invalidate ();
		this->is_valid_ = true;
	}
	void compute () override { throw std::runtime_error ("Param node unset"); }
};

// TODO keep or trash
template <typename T> struct AddNode final : public ValuedNode<T> {
	void add_dep (Value<T> v) {
		auto p = v.get ();
		p->register_node (this);
		deps.emplace_back (std::move (v));
		this->invalidate ();
	}
	void compute () override {
		T sum{};
		for (auto & dep : deps)
			sum += dep.get_value ();
	}
	std::vector<Value<T>> deps;
};
}

//////////////////////////////////// Operation ids ///////////////////////////////////

using OperationId = int;
static OperationId new_operation_id () {
	static OperationId counter{0};
	return counter++;
}

template <typename Operation> OperationId get_operation_id () {
	static OperationId id{new_operation_id ()};
	return id;
}

//////////////////////////////////// Topology //////////////////////////////////////

using NodeIndex = int;
constexpr NodeIndex invalid_index{-1};
static NodeIndex new_global_node_index () {
	static NodeIndex counter{0};
	return counter++;
}

class Node;
using Tree = std::shared_ptr<const Node>;

class Node {
public:
	Node (NodeIndex id, std::vector<Tree> childrens) : id_ (id), childrens_ (std::move (childrens)) {
		std::sort (childrens_.begin (), childrens_.end (),
		           [](const Tree & a, const Tree & b) { return a->id () < b->id (); });
	}
	bool is_leaf () const { return childrens_.empty (); }
	NodeIndex id () const { return id_; }
	const std::vector<Tree> & childrens () const { return childrens_; }
	bool operator== (const Node & o) const {
		return std::tie (id_, childrens_) == std::tie (o.id_, o.childrens_);
	}

private:
	NodeIndex id_;
	std::vector<Tree> childrens_;
};
bool operator== (const Tree & a, const Tree & b) {
	return *a == *b;
}

Tree make_tree (std::vector<Tree> childrens = {}, NodeIndex id = invalid_index) {
	if (id == invalid_index)
		id = new_global_node_index ();
	return std::make_shared<Node> (id, std::move (childrens));
}

void print_tree (std::ostream & os, Tree t, int depth = 0) {
	os << std::string (depth * 4, ' ') << t->id () << "\n";
	for (auto c : t->childrens ())
		print_tree (os, c, depth + 1);
}

//////////////////////////////////// DF registry ///////////////////////////////////

struct Key {
	const OperationId op_id;
	const Tree application_point;

	bool operator== (const Key & k) const {
		return std::tie (op_id, application_point) == std::tie (k.op_id, k.application_point);
	}
	struct Hash {
		std::size_t operator() (const Key & k) const { return std::hash<OperationId>{}(k.op_id); }
	};
};

template <typename OperationType> Key make_key () {
	return Key{get_operation_id<OperationType> ()};
}

class Registry {
public:
private:
	std::unordered_map<Key, DF::Node, typename Key::Hash> registry_;
};

int main () {
	// Building a tree.
	auto a = make_tree ();
	auto b = make_tree ();
	auto c = make_tree ({a, b});
	auto d = make_tree ();
	auto e = make_tree ();
	auto f = make_tree ();
	auto g = make_tree ({d, e, f});
	auto root = make_tree ({c, g});
	auto alt_g = make_tree ({d, e, f});
	auto alt_root = make_tree ({b, alt_g}, root->id ());
	std::cout << "root\n";
	print_tree (std::cout, root);
	std::cout << "alt_root\n";
	print_tree (std::cout, alt_root);
	std::cout << "end tree printing\n";
	// End tree

	std::cout << (g == alt_root) << "\n";
	std::cout << get_operation_id<int> () << " " << get_operation_id<int> () << " "
	          << get_operation_id<bool> () << "\n";

	return 0;
}

#if 0
static constexpr int invalid = -1;

class Tree {
public:
	struct Node {
		int id{invalid};
		int father{invalid};
		std::vector<int> childrens{};
	};

	// copy and move auto

	Node & node (int id) { return nodes_[id]; }
	const Node & node (int id) const { return nodes_[id]; }

	int add_node (std::vector<int> childrens) {
		auto id = int(nodes_.size ());
		for (auto i : childrens)
			nodes_[i].father = id;
		auto node = Node{};
		node.id = id;
		node.childrens = std::move (childrens);
		nodes_.emplace_back (std::move (node));
		return id;
	}

	int root_id () const { return int(nodes_.size ()) - 1; } // invalid if empty

	void print (std::ostream & os, int from_node, int depth = 0) const {
		auto & n = node (from_node);
		os << std::string (depth * 4, ' ') << n.id << "\n";
		for (auto i : n.childrens)
			print (os, i, depth + 1);
	}

private:
	std::vector<Node> nodes_;
};
#endif
