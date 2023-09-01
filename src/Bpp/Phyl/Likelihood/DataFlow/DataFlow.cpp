//
// File: DataFlow.cpp
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-05-30 00:00:00
// Last modified: 2017-05-30 00:00:00
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  
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
#include <algorithm>
#include <fstream> // debug
#include <functional> // std::hash
#include <iomanip>
#include <iostream> // debug
#include <ostream> // debug
#include <regex>
#include <stack> // invalidate/compute recursively + debug
#include <type_traits> // DotOptions flags
#include <typeinfo>
#include <unordered_set> // debug
#include <set>

#include "DataFlow.h"
#include "DataFlowNumeric.h"

/* std::type_info::name() returns a "mangled" type name, not very readable.
 * Compilers can optionally provide an ABI header cxxabi.h.
 * This header contain the demangle function, which makes the type name readable.
 * By testing on godbolt, all versions of gcc and clang have this header.
 */
#if defined(__GNUC__) || defined(__clang__)
#include <cstdlib>
#include <cxxabi.h>
static std::string demangle (const char* name)
{
  int status{};
  std::unique_ptr<char, void (*)(void*)> res{abi::__cxa_demangle (name, nullptr, nullptr, &status),
                                             std::free};
  return status == 0 ? res.get () : name;
}
#else
static std::string demangle (const char* name) { return name; }
#endif

namespace bpp
{
std::string prettyTypeName (const std::type_info& ti) { return demangle (ti.name ()); }
} // namespace bpp

namespace bpp
{
/*****************************************************************************
 * Error & dependency check functions.
 */
void failureComputeWasCalled (const std::type_info& nodeType)
{
  throw Exception (prettyTypeName (nodeType) + ": compute() was called");
}

void failureNodeConversion (const std::type_info& handleType, const Node_DF& node)
{
  throw Exception (prettyTypeName (handleType) + " cannot store: " + prettyTypeName (typeid (node)));
}

void failureDependencyNumberMismatch (const std::type_info& contextNodeType, std::size_t expectedSize,
                                      std::size_t givenSize)
{
  throw Exception (prettyTypeName (contextNodeType) + ": expected " + std::to_string (expectedSize) +
                   " dependencies, got " + std::to_string (givenSize));
}

void failureEmptyDependency (const std::type_info& contextNodeType, std::size_t depIndex)
{
  throw Exception (prettyTypeName (contextNodeType) + ": " + std::to_string (depIndex) +
                   "-th dependency is empty (nullptr)");
}

void failureDependencyTypeMismatch (const std::type_info& contextNodeType, std::size_t depIndex,
                                    const std::type_info& expectedType,
                                    const std::type_info& givenNodeType)
{
  throw Exception (prettyTypeName (contextNodeType) + ": expected class derived from " +
                   prettyTypeName (expectedType) + " as " + std::to_string (depIndex) +
                   "-th dependency, got " + prettyTypeName (givenNodeType));
}

void checkDependencyVectorSize (const std::type_info& contextNodeType, const NodeRefVec& deps,
                                std::size_t expectedSize)
{
  auto size = deps.size ();
  if (size != expectedSize)
  {
    failureDependencyNumberMismatch (contextNodeType, expectedSize, size);
  }
}

void checkDependencyVectorMinSize (const std::type_info& contextNodeType, const NodeRefVec& deps,
                                   std::size_t expectedMinSize)
{
  auto size = deps.size ();
  if (size < expectedMinSize)
  {
    failureDependencyNumberMismatch (contextNodeType, expectedMinSize, size);
  }
}

void checkDependenciesNotNull (const std::type_info& contextNodeType, const NodeRefVec& deps)
{
  for (std::size_t i = 0; i < deps.size (); ++i)
  {
    if (!deps[i])
    {
      failureEmptyDependency (contextNodeType, i);
    }
  }
}

void checkNthDependencyNotNull (const std::type_info& contextNodeType, const NodeRefVec& deps, std::size_t index)
{
  if (!deps[index])
    failureEmptyDependency (contextNodeType, index);
}

/*****************************************************************************
 * Node impl.
 */
Node_DF::Node_DF (const NodeRefVec& dependenciesArg) : dependencyNodes_ (dependenciesArg)
{
  // dependencyNodes_.erase(std::remove_if(dependencyNodes_.begin(),dependencyNodes_.end(),[](NodeRef nr){return nr==0;}), dependencyNodes_.end());
  for (auto& n : dependencyNodes_)
  {
    if (n)
      n->registerNode (this);
  }
}
Node_DF::Node_DF (NodeRefVec&& dependenciesArg) : dependencyNodes_ (std::move (dependenciesArg))
{
  // dependencyNodes_.erase(std::remove_if(dependencyNodes_.begin(),dependencyNodes_.end(),[](NodeRef nr){return nr==0;}), dependencyNodes_.end());
  for (auto& n : dependencyNodes_)
  {
    if (n)
      n->registerNode (this);
  }
}

Node_DF::~Node_DF ()
{
  for (auto& n : dependencyNodes_)
  {
    if (n)
      n->unregisterNode (this);
  }
}

std::string Node_DF::description() const
{
  std::string nodeType = prettyTypeName (typeid (*this));
  // Shorten displayed name by removing namespaces.
  nodeType = std::regex_replace(nodeType, std::regex("bpp::"), "");
  nodeType = std::regex_replace(nodeType, std::regex("Eigen::"), "");
  nodeType = std::regex_replace(nodeType, std::regex("std::"), "");
  nodeType = std::regex_replace(nodeType, std::regex("__cxx..::"), "");
  nodeType = std::regex_replace(nodeType, std::regex("(char_traits<[^>]*>)"), "char");
  nodeType = std::regex_replace(nodeType, std::regex("(allocator<[^>]*>)"), "");
  nodeType = std::regex_replace(nodeType, std::regex("(Matrix<)([^>]*>)"), "Matrix");
  nodeType = std::regex_replace(nodeType, std::regex("(basic_string<)([^>]*>)"), "string");
  return nodeType;
}

std::string Node_DF::debugInfo() const { return {}; }

bool Node_DF::hasNumericalProperty(NumericalProperty) const { return false; }

bool Node_DF::compareAdditionalArguments(const Node_DF&) const { return false; }
std::size_t Node_DF::hashAdditionalArguments() const { return 0; }

NodeRef Node_DF::derive(Context&, const Node_DF&)
{
  throw Exception ("Node does not support derivation: " + description ());
}

NodeRef Node_DF::recreate(Context&, NodeRefVec&&)
{
  throw Exception ("Node does not support recreate(deps): " + description ());
}

void Node_DF::computeRecursively ()
{
  // Compute the current node (and dependencies recursively) if needed
  if (isValid ())
    return;

  // Discover then recompute needed nodes
  std::stack<Node_DF*> nodesToVisit;
  std::stack<Node_DF*> nodesToRecompute;
  nodesToVisit.push (this);
  while (!nodesToVisit.empty ())
  {
    auto* n = nodesToVisit.top();
    nodesToVisit.pop();
    if (!n->isValid())
    {
      nodesToRecompute.push(n);
      for (auto& dep : n->dependencies())
      {
        if (dep)
          nodesToVisit.push(dep.get());
      }
    }
  }
  while (!nodesToRecompute.empty())
  {
    auto* n = nodesToRecompute.top();
    nodesToRecompute.pop ();
    if (!n->isValid())
    {
      n->compute();
      n->makeValid();
    }
  }
}

void Node_DF::invalidateRecursively () noexcept
{
  if (!isValid ())
    return;
  std::stack<Node_DF*> nodesToInvalidate;
  nodesToInvalidate.push (this);
  while (!nodesToInvalidate.empty ())
  {
    auto* n = nodesToInvalidate.top ();
    nodesToInvalidate.pop ();
    if (n->isValid ())
    {
      n->makeInvalid ();
      for (auto* dependent : n->dependentNodes_)
      {
        nodesToInvalidate.push (dependent);
      }
    }
  }
}

void Node_DF::registerNode (Node_DF* n)
{
  dependentNodes_.emplace_back (n);
}

void Node_DF::unregisterNode (const Node_DF* n)
{
  dependentNodes_.erase (std::remove (dependentNodes_.begin (), dependentNodes_.end (), n),
                         dependentNodes_.end ());
}

/*****************************************************************************
 * Free functions.
 */
bool isTransitivelyDependentOn (const Node_DF& searchedDependency, const Node_DF& node)
{
  std::stack<const Node_DF*> nodesToVisit;
  nodesToVisit.push (&node);
  while (!nodesToVisit.empty ())
  {
    auto* current = nodesToVisit.top ();
    nodesToVisit.pop ();
    if (current == &searchedDependency)
      return true;
    for (auto& dep : current->dependencies ())
    {
      nodesToVisit.push (dep.get ());
    }
  }
  return false;
}

NodeRef recreateWithSubstitution(Context& c, const NodeRef& node,
                                  const std::unordered_map<const Node_DF*, NodeRef>& substitutions)
{
  if (node == 0)
    return node;
  auto it = substitutions.find(node.get());
  if (it != substitutions.end())
  {
    // Substitute sub tree.
    return it->second;
  }
  else if (node->dependencies().empty())
  {
    // Leaves do no support rebuild: just copy them
    return node;
  }
  else
  {
    // Recursion : only rebuild if dependencies have changed
    NodeRefVec recreatedDeps(node->nbDependencies ());
    for (std::size_t i = 0; i < recreatedDeps.size (); ++i)
    {
      recreatedDeps[i] = recreateWithSubstitution(c, node->dependency(i), substitutions);
    }
    if (recreatedDeps == node->dependencies ())
    {
      return node;
    }
    else
    {
      return node->recreate(c, std::move(recreatedDeps));
    }
  }
}

/*****************************************************************************
 * Context.
 */

Context::Context() : nodeCache_(), zero_(ConstantZero<size_t>::create(*this, Dimension<size_t>()))
{
}

                           
NodeRef Context::cached(NodeRef&& newNode)
{
  assert (newNode != nullptr);
  // Try inserting it, which will fail if already present and return the old one
  auto r = nodeCache_.emplace(std::move (newNode));
  return r.first->ref;
}

NodeRef Context::cached (NodeRef& newNode)
{
  assert (newNode != nullptr);
  // First remove this object from the set if it is already
  for (auto it = nodeCache_.begin(); it != nodeCache_.end();)
  {
    if (it->ref == newNode)
      it = nodeCache_.erase(it);
    else
      ++it;
  }

  // Try inserting it, which will fail if already present and return the old one
  auto r = nodeCache_.emplace (newNode);
  return r.first->ref;
}

std::vector<const Node_DF*> Context::getAllNodes() const
{
  std::vector<const Node_DF*> ret;
  for (const auto& it : nodeCache_)
    ret.push_back(it.ref.get());
  return(ret);
}
    

bool Context::erase(const NodeRef& node)
{
  if (!node)
    return false;

  std::set<NodeRef> sNodes;
  sNodes.emplace(node);

  bool flag=true;
  bool ret=false;
  while(flag)
  {
    flag=false;
    for (auto n : sNodes)
    {
      if (n->nbDependentNodes()==0) {
        for (auto it = nodeCache_.begin(); it != nodeCache_.end();)
        {
          if (it->ref == n)
          {
            sNodes.erase(n);
            for (auto dep:n->dependencies())
            {
              if (!dep)
                continue;
              sNodes.emplace(dep);
              dep->unregisterNode (n.get());
            }
            it=nodeCache_.erase(it);
            flag=true;
            ret=true;
          }
          else
            ++it;
        }
      }
      if (flag)
        break;
    }
  }

  return ret;
}

/* Compare/hash the triplet (type, deps, additionalArgs).
 * type and deps are available directly from the Node*.
 * additionalArgs is handled through the two virtual methods.
 */
bool Context::CachedNodeRef::operator==(const CachedNodeRef& other) const
{
  const auto& lhs = *this->ref;
  const auto& rhs = *other.ref;
  return typeid (lhs) == typeid (rhs) && lhs.dependencies () == rhs.dependencies () &&
         lhs.compareAdditionalArguments (rhs);
}

std::size_t Context::CachedNodeRefHash::operator()(const CachedNodeRef& ref) const
{
  const auto& node = *ref.ref;
  std::size_t seed = typeid (node).hash_code ();
  for (const auto& dep : node.dependencies ())
  {
    combineHash (seed, dep);
  }
  combineHash (seed, node.hashAdditionalArguments ());
  return seed;
}

/*****************************************************************************
 * output dataflow graph to dot format.
 */

// Make enum behave as flags.
DotOptions operator|(DotOptions a, DotOptions b)
{
  using IntType = typename std::underlying_type<DotOptions>::type;
  return static_cast<DotOptions>(static_cast<IntType>(a) | static_cast<IntType>(b));
}
bool operator&(DotOptions a, DotOptions b)
{
  using IntType = typename std::underlying_type<DotOptions>::type;
  return static_cast<IntType>(a) & static_cast<IntType>(b);
}

// Escape text in a dot box label
static std::string dotLabelEscape (const std::string& s)
{
  std::string result;
  result.reserve (s.size ());
  static const char toEscape[] = "<>|{} ";
  for (const char c : s)
  {
    if (std::any_of (std::begin (toEscape), std::end (toEscape),
                     [c](const char c2) {
        return c == c2;
      }))
    {
      result.push_back ('\\');
    }
    result.push_back (c);
  }
  return result;
}

// Dot node id for dataflow Node: 'N' + hash as a string
static std::string dotIdentifier (const Node_DF& node)
{
  return 'N' + std::to_string (std::hash<const Node_DF*>{}(&node));
}

// Write line with node representation
static void writeDotNode (std::ostream& os, const Node_DF& node, DotOptions opt)
{
  os << '\t' << dotIdentifier (node);
  if (opt & DotOptions::DetailedNodeInfo)
  {
    os << " [style=filled, shape=Mrecord,label=\"{" << dotLabelEscape (node.description ())
       << "| valid=" << node.isValid () << ' ' << dotLabelEscape (node.debugInfo ()) << "}\"";
  }
  else
  {
    os << " [style=filled, shape=" << node.shape() << ",label=\"" << dotLabelEscape (node.description ()) << "\"";
  }
  os << ",fillcolor=\"" << node.color() << "\"]";
  os << ";\n";
}

// Write line with edge representation for n-th dependency of from
static void writeDotEdge (std::ostream& os, const Node_DF& from, std::size_t depIndex, DotOptions opt)
{
  const auto& to = *from.dependency (depIndex);
  os << '\t' << dotIdentifier (from) << " -> " << dotIdentifier (to);
  os << " [";
  if (opt & DotOptions::ShowDependencyIndex)
  {
    os << "label=\"" << depIndex << "\",";
  }

  auto fc = from.color() != "white" ? from.color() : "black";
  // auto tc=to.color()!="white"?to.color():"black";

  // std::vector<double> vpat={0.05,0.1,0.15,0.2};

  // std::string pat;
  // for (size_t i=0;i<vpat.size();i++)
  //   pat += fc + ";" + std::to_string(vpat[i]) + ":" + tc + ";" + std::to_string(vpat[vpat.size()-i-1])+":";

  os << "color=\"" << fc << "\"";
  os << " ]";
  os << ";\n";
}

// Write dot lines for graph structure, starting from the given entry points.
static void writeGraphStructure (std::ostream& os, const std::vector<const Node_DF*>& entryPoints,
                                 DotOptions opt)
{
  std::stack<const Node_DF*> nodesToVisit;
  std::unordered_set<const Node_DF*> discoveredNodes;

  const auto discover = [&nodesToVisit, &discoveredNodes](const Node_DF* n) {
                          if (!n)
                            return;
                          const bool discovered = discoveredNodes.find (n) != discoveredNodes.end ();
                          if (!discovered)
                          {
                            nodesToVisit.emplace (n);
                            discoveredNodes.emplace (n);
                          }
                        };

  for (const auto* n : entryPoints)
  {
    discover (n);
  }

  while (!nodesToVisit.empty ())
  {
    const auto* node = nodesToVisit.top ();

    nodesToVisit.pop ();
    writeDotNode (os, *node, opt);

    if (opt & DotOptions::FollowUpwardLinks)
    {
      for (const auto* dependent : node->dependentNodes ())
      {
        discover (dependent);
      }
    }
    for (std::size_t index = 0; index < node->nbDependencies (); ++index)
    {
      if (node->dependency(index))
      {
        writeDotEdge (os, *node, index, opt);
        discover (node->dependency (index).get ());
      }
    }
  }
}

void writeGraphToDot (std::ostream& os, const std::vector<const Node_DF*>& nodes, DotOptions opt)
{
  os << "digraph {\n";
  writeGraphStructure (os, nodes, opt);
  os << "}\n";
}

void writeGraphToDot (const std::string& filename, const std::vector<const Node_DF*>& nodes,
                      DotOptions opt)
{
  std::ofstream file{filename};
  writeGraphToDot (file, nodes, opt);
}
} // namespace bpp
