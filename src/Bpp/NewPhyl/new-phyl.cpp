#include <algorithm>
#include <memory>
#include <vector>

#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_map>

template <typename T> std::string basic_type_name ();
#define SET_TYPE_NAME(t)                                                                           \
	template <> std::string basic_type_name<t> () { return #t; }
SET_TYPE_NAME (int);

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
	const std::shared_ptr<BaseNode> & get_pointer () const { return p_node_; }

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

// Value type for ParamNode
template <typename T> class Param {
public:
	Param (std::shared_ptr<ParamNode<T>> p) : p_node_ (std::move (p)) {}

	Param (const Node & node)
	    : p_node_ (std::dynamic_pointer_cast<ParamNode<T>> (node.get_pointer ())) {
		if (!p_node_)
			throw std::runtime_error ("bad dynamic cast");
	}

	Param (T t) : Param (std::make_shared<ParamNode<T>> ()) { set_value (t); }

	void set_value (T t) const { p_node_->set_value (t); }
	const T & get_value () const { return p_node_->get_value (); }

	Node as_node () const { return Node (std::static_pointer_cast<BaseNode> (p_node_)); }

private:
	std::shared_ptr<ParamNode<T>> p_node_;
};

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

using DataSet = int;
class Registry;

struct Key {
	const OperationId op_id;
	const Tree application_point;
	const DataSet data_set;

	bool operator== (const Key & k) const {
		return std::tie (op_id, application_point, data_set) ==
		       std::tie (k.op_id, k.application_point, k.data_set);
	}
	struct Hash {
		std::size_t operator() (const Key & k) const {
			return std::hash<OperationId>{}(k.op_id) ^ (std::hash<Tree>{}(k.application_point) << 1) ^
			       (std::hash<DataSet>{}(k.data_set) << 2);
		}
	};
};

class Registry {
public:
	// Data set and data points
	DataSet create_data_set () { return counter++; }
	template <typename T>
	DF::Param<T> create_parameter (const Tree & tree, const DataSet & data_set, T t = T ()) {
		auto res = registry_.emplace (make_key<DF::ParamNode<T>> (tree, data_set),
		                              DF::Param<T> (t).as_node ());
		if (!res.second)
			throw std::runtime_error ("param node already created !");
		return DF::Param<T> (res.first->second);
	}
	template <typename T> DF::Param<T> get_parameter (const Tree & tree, const DataSet & data_set) {
		auto it = registry_.find (make_key<DF::ParamNode<T>> (tree, data_set));
		if (it == registry_.end ())
			throw std::runtime_error ("param not found");
		return DF::Param<T> (it->second);
	}

  template<typename DFNode> std::shared_ptr<DFNode> get_or_create_node_raw (const Tree & tree, const DataSet & data_set) {
    // TODO
  }

private:
	template <typename OperationType>
	Key make_key (const Tree & tree, const DataSet & data_set) const {
		return Key{get_operation_id<OperationType> (), tree, data_set};
	}

	DataSet counter{0};
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

	Registry registry;

	auto data1 = registry.create_data_set ();
	registry.create_parameter (a, data1, 1);
	registry.create_parameter (b, data1, 2);
	registry.create_parameter (d, data1, 3);
	registry.create_parameter (e, data1, 4);
	registry.create_parameter (f, data1, 5);

	return 0;
}
