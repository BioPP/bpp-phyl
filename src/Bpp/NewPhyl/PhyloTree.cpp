//
// File: ConfiguredPhyloTree.cpp
// Created by: Laurent Guéguen
// Created on: Wed Jul 11 20:36 2012
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

#include <Bpp/Exceptions.h>
#include <Bpp/NewPhyl/PhyloTree.h>
#include <Bpp/NewPhyl/Parametrizable.h>
#include <Bpp/Phyl/NewLikelihood/ParametrizablePhyloTree.h>

using namespace bpp;
using namespace std;
using namespace dataflow;

using namespace std;

namespace bpp {
  namespace dataflow {

    // Tree node

    ConfiguredPhyloTree::ConfiguredPhyloTree (NodeRefVec && deps, std::unique_ptr<ParametrizablePhyloTree> && tree)
      : Value<ParametrizablePhyloTree*> (std::move (deps), tree.get ()), tree_(std::move(tree)) {}

    ConfiguredPhyloTree::~ConfiguredPhyloTree () = default;

    std::string ConfiguredPhyloTree::description () const { return "Tree";}

    std::string ConfiguredPhyloTree::debugInfo () const {
      return "nbLeaves=" + std::to_string (tree_->getNumberOfLeaves());
    }

    const std::string & ConfiguredPhyloTree::getParameterName (std::size_t index) {
      return tree_->getParameters ()[index].getName ();
    }
    
    std::size_t ConfiguredPhyloTree::getParameterIndex (const std::string & name) {
      return static_cast<std::size_t> (tree_->getParameters ().whichParameterHasName (name));
    }

    // Tree node additional arguments = (type of bpp::ParametrizablePhyloTree).
    // Everything else is determined by the node dependencies.
    bool ConfiguredPhyloTree::compareAdditionalArguments (const Node & other) const {
      const auto * derived = dynamic_cast<const Self *> (&other);
      if (derived == nullptr) {
        return false;
      } else {
        const auto & thisTree = *tree_;
        const auto & otherTree = *derived->tree_;
        return typeid (thisTree) == typeid (otherTree);
      }
    }
    
    std::size_t ConfiguredPhyloTree::hashAdditionalArguments () const {
      const auto & bppTree = *tree_;
      return typeid (bppTree).hash_code ();
    }

    NodeRef ConfiguredPhyloTree::recreate (Context & c, NodeRefVec && deps) {
      auto m = ConfiguredParametrizable::createConfigured<Target, Self> (c, std::move (deps), std::unique_ptr<Target>{tree_->clone ()});
      m->config = this->config; // Duplicate derivation config
      return m;
    }


    void ConfiguredPhyloTree::compute () {
      // Update each internal tree bpp::Parameter with the dependency
      auto & parameters = tree_->getParameters ();
      const auto nbParameters = this->nbDependencies ();
      for (std::size_t i = 0; i < nbParameters; ++i) {
        auto & v = accessValueConstCast<double> (*this->dependency (i));
        auto & p = parameters[i];
        if (p.getValue () != v) {
          // TODO improve bpp::Parametrizable interface to change values by index.
          tree_->setParameterValue (tree_->getParameterNameWithoutNamespace (p.getName ()), v);
        }
      }
    }
  }
}
