//
// File: ConfiguredPhyloTree.h
// Created by: Laurent Guéguen
// Created on: mardi 30 octobre 2018, à 15h 08
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

#ifndef _CONFIGURED_PHYLO_TREE_H_
#define _CONFIGURED_PHYLO_TREE_H_

#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/DataFlowNumeric.h>
#include <Bpp/Exceptions.h>
#include <functional>
#include <unordered_map>

namespace bpp
{
  class ParametrizablePhyloTree;

  namespace dataflow
  {
    
    /**
     * PhyloTree with DF Phylo Branches.
     *
     **/

    class ConfiguredPhyloTree:
      public Value<ParametrizablePhyloTree*>
    {
    public:
      using Self = ConfiguredPhyloTree;
      using Target = ParametrizablePhyloTree;

      ConfiguredPhyloTree (NodeRefVec && deps, std::unique_ptr<ParametrizablePhyloTree> && tree);
      ~ConfiguredPhyloTree ();

      std::string description () const final;
      std::string debugInfo () const final;

      /// Return the index of parameter with the given non namespaced name (or throw).
      std::size_t getParameterIndex (const std::string & name);
      /// Return the non namespaced name for parameter at the given index.
      const std::string & getParameterName (std::size_t index);

      bool compareAdditionalArguments (const Node & other) const;
      
      std::size_t hashAdditionalArguments () const;
      
      /// Configuration for numerical derivation of computation nodes using this Model.
      NumericalDerivativeConfiguration config;

      NodeRef recreate (Context & c, NodeRefVec && deps) final;
      
    private:
      void compute ();

      std::unique_ptr<ParametrizablePhyloTree> tree_;
      
    };
    
  } //end of namespace dataflow
  
} //end of namespace bpp

#endif //_CONFIGURED_PHYLO_TREE_H_

