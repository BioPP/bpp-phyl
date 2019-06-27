//
// File: ProcessTree.h
// Created by: Laurent Guéguen
// Created on: jeudi 15 septembre 2016, à 06h 40
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

#ifndef _PROCESS_TREE_H_
#define _PROCESS_TREE_H_

#include <Bpp/Phyl/NewLikelihood/ParametrizablePhyloTree.h>

#include <Bpp/Graph/AssociationTreeGraphImplObserver.h>

#include "Model.h"
#include "Parameter.h"

//From the stl:
#include <string>

namespace bpp
{

  namespace dataflow
  {

    /** Helper: create a map with mutable dataflow nodes for each
     *  branch of the tree.
     *  The map is indexed by branch ids.
     */

    /* Map (branch id, ConfiguredParameter of Branch Length) */
     
    using BrLenMap = std::map<uint, std::shared_ptr<ConfiguredParameter>>;

    /* Map (branch id, ConfiguredModel of Branch Length) */

    struct ModelAssign
    {
      /* ConfiguredModel on this branch */
      std::shared_ptr<ConfiguredModel> model_;

      /* vector of allowed submodels if model_ is mixed*/
      /* empty vector is model_ is not mixed, or considered as not
       * mixed (ie its transition matrix is a mixture of transition
       * matrices of submodels). */
   
      std::vector<size_t> modelNum_;

      ModelAssign() : model_(), modelNum_(){};

      ModelAssign(std::shared_ptr<ConfiguredModel> model) : model_(model), modelNum_(){};

      ModelAssign(std::shared_ptr<ConfiguredModel> model, std::vector<size_t> nMod) : model_(model), modelNum_(nMod){};
      
    };
      
    using ModelMap = std::map<uint, ModelAssign>;

    /** Helper: create a map with mutable dataflow nodes for each
     *  branch of the tree.
     *  The map is indexed by branch ids.
     */
    
    // Branch specific DataFlow objects
    struct ProcessEdge
    {
      std::shared_ptr<ConfiguredParameter> brlen_;
      //NodeRef brlen_;

      std::shared_ptr<ConfiguredModel> model_;

      /*
       *@brief Optional number of submodels, in case model_ is mixed
       * and a submodel is used.
       *
       */

      std::shared_ptr<NumericConstant<size_t>> nMod_;
      
    public:

      ProcessEdge(std::shared_ptr<ConfiguredParameter> brlen=0,
                  std::shared_ptr<ConfiguredModel> model=0,
                  std::shared_ptr<NumericConstant<size_t>> nMod=0) : brlen_(brlen), model_(model), nMod_(nMod){};

      std::shared_ptr<ConfiguredModel> getModel()
      {
        return model_;
      }

      void setModel(std::shared_ptr<ConfiguredModel> model)
      {
        model_=model;
      }

      std::shared_ptr<ConfiguredParameter> getBrLen()
      {
        return brlen_;
      }

      void setBrLen(std::shared_ptr<ConfiguredParameter> brlen)
      {
        brlen_=brlen;
      }

      void setNMod(std::shared_ptr<NumericConstant<size_t>> nMod)
      {
        nMod_=nMod;
      }

      std::shared_ptr<NumericConstant<size_t>> getNMod() 
      {
        return nMod_;
      }

    };

    // Node specific DataFlow objects
    class ProcessNode:
      public PhyloNode
    {
      /*
       * Mixing probabilities used at Node (if needed, such as for Mixture Nodes)
       *
       */
      
      std::shared_ptr<ProbabilitiesFromMixedModel> proba_;
      
    public:

      ProcessNode(const PhyloNode& node) :
        PhyloNode(node),
        proba_(0){}

      std::shared_ptr<ProbabilitiesFromMixedModel> getProba()
      {
        return proba_;
      }

      void setProba(std::shared_ptr<ProbabilitiesFromMixedModel> proba)
      {
        proba_=proba;
      }

    };

    
    class ProcessTree : public AssociationTreeGlobalGraphObserver<ProcessNode,ProcessEdge>
    {
      using Self = ProcessTree;

      Context context_;
      
    public:

      /*
       * @brief Build a ProcessTree with same topology as given
       * ParametrizablePhyloTree, and matching ConfiguredParameter
       * BrLen.
       *
       */
       
      ProcessTree(Context& context, const ParametrizablePhyloTree& tree, const BrLenMap& vrefmap);

      /*
       * @brief Build a ProcessTree given a basic topology and map of
       * models to branches. So the resulting topology may be
       * different from the given ParametrizablePhyloTree.
       *
       */

      ProcessTree(Context& context, const ParametrizablePhyloTree& tree, const BrLenMap& vrefmap,
                  ModelMap& modelmap) :
        AssociationTreeGlobalGraphObserver<ProcessNode,ProcessEdge>(tree.isRooted()), context_(context)
      {
        buildUnderNode_(tree, vrefmap, modelmap, tree.getRoot());
      }

    private:
      /*
       * Build the tree under a node of the ParametrizablePhyloTree
       * (which node is a priori a speciation node).
       *
       * In option, the new node will be linked from a father newNode,
       * through a newEdge.
       */
      
      void buildUnderNode_(const ParametrizablePhyloTree& tree, const BrLenMap& vrefmap, ModelMap& modelmap, std::shared_ptr<PhyloNode> node, std::shared_ptr<ProcessNode> newFather=0, std::shared_ptr<ProcessEdge> newEdge=0);

      /*
       * Build the tree under an edge of the ParametrizablePhyloTree.
       * A new branch will be built, which will support a model.
       *
       * If the model is a Mixture model met elsewhere, a mixture node
       * is built (at distance 0 from the newNode), and the
       * construction of ProcessEdge branches are set in the following
       * branches, for all submodels.
       *
       */
      
      void buildUnderEdgeFromNode_(const ParametrizablePhyloTree& tree, const BrLenMap& vrefmap, ModelMap& modelmap, std::shared_ptr<PhyloBranchParam> oldEdge, std::shared_ptr<ProcessNode> newNode);

    public:
      /*
       * Assigns new ConfiguredParameters for BrLen to relevant edges
       *
       * Names of the BrLen configured parameters are used to recover
       * the matching branches. So parameters in the new BrLenMap must
       * match the former parameters of the copied ProcessTree,
       * otherwise nothing is done and former parameters are shared.
       * 
       *
       */
      
      ProcessTree(const ProcessTree& tree, BrLenMap&& vrefmap) :
        AssociationTreeGlobalGraphObserver<ProcessNode,ProcessEdge>(tree)
      {
        std::vector<uint> vEdgesId=tree.getAllEdgesIndexes();

        for (auto& index:vEdgesId)
        {
          if (vrefmap.find(index)==vrefmap.end())
            throw Exception("ProcessTree::ProcessTree : unknown edge id in Map : " + std::to_string(index));
          getEdge(index)->setBrLen(vrefmap.at(index));
        }
      }

      ProcessTree* clone() const {
        throw Exception("ProcessTree::clone should not be called.");
      }

      ProcessTree(const ProcessTree& pTree): 
        AssociationTreeGlobalGraphObserver<ProcessNode,ProcessEdge>(pTree.getGraph())
      {
        throw Exception("ProcessTree::ProcessTree should not be called.");
      }
      
      ProcessTree& operator=(const ProcessTree& pTree)
      {
        throw Exception("ProcessTree::operator= should not be called.");
        //AssociationTreeGlobalGraphObserver<ProcessNode,Value<T>>::operator=(pTree);
        return *this;
      }

    };

    
    /**************************************/
    /* Methods to create/modify BrLenMap objects */
    
    inline BrLenMap createBrLenMap(Context & c, const ParametrizablePhyloTree& tree)
    {
      std::vector<std::shared_ptr<PhyloBranchParam> > vB=tree.getAllEdges();

      BrLenMap map;

      for (auto& branch:vB)
      {
        auto brl=NumericMutable<double>::create(c, branch->getLength());
        map.emplace (tree.getEdgeIndex(branch),
                     ConfiguredParameter::create (c, {std::move(brl)}, branch->getParameters()[0]));
      }

      return map;
    }

    inline BrLenMap createBrLenMap(Context & c, const ProcessTree& tree)
    {
      // Ids of the branches
      
      BrLenMap map;

      auto aEit=tree.allEdgesIterator();
        
      while (!aEit->end())
      {
        auto edge=**aEit;
        if (edge->getBrLen())
          map.emplace (tree.getEdgeIndex(edge),
                       std::dynamic_pointer_cast<ConfiguredParameter>(edge->getBrLen()->recreate(c,{edge->getBrLen()->dependency(0)})));
        else
          map.emplace (tree.getEdgeIndex(edge), std::shared_ptr<ConfiguredParameter>(0));
        
        aEit->next();
      }
      
      return map;
    }

    
    inline BrLenMap multiplyBrLenMap(Context & c, const ProcessTree& tree, ValueRef<double>&& rate)
    {
      // Ids of the branches

      BrLenMap map;

      auto aEit=tree.allEdgesIterator();
        
      while (!aEit->end())
      {
        auto edge=**aEit;
        
        if (edge->getBrLen())
        {
          auto mulref = CWiseMul<double, std::tuple<double, double>>::create (c, {edge->getBrLen()->dependency(0), rate}, Dimension<double>());

          map.emplace (tree.getEdgeIndex(edge),
                       std::dynamic_pointer_cast<ConfiguredParameter>(edge->getBrLen()->recreate(c,{std::move(mulref)})));
        }
        else
          map.emplace (tree.getEdgeIndex(edge), std::shared_ptr<ConfiguredParameter>(0));
//map.emplace (tree.getEdgeIndex(edge), NodeRef(0));
        
        aEit->next();
      }

      return map;
    }


    
    /***************************************************/
    /** Create a new Tree Node. Tree Node parameters are get from
     * ConfiguredParameters PREVIOUSLY built and stored in a
     * ParameterList*/
    
    inline std::shared_ptr<ProcessTree> makeProcessTree(Context& context, const ParametrizablePhyloTree& parTree, ModelMap& modelmap, ParameterList& parList, const std::string& suff = "")
    {
      std::vector<std::shared_ptr<PhyloBranchParam> > vB=parTree.getAllEdges();
      BrLenMap mapBr;
      
      for (auto& branch:vB)
      {
        const auto& bp=branch->getParameters()[0];
        
        std::string name=bp.getName()+suff;
        if (!parList.hasParameter(name) && suff=="")
        {
          if (!parList.hasParameter(bp.getName()+"_1"))
            throw Exception("makeTreeNode: unknown ConfiguredParameter " + name);
          else name=bp.getName()+"_1";
        }

        auto confPar=dynamic_cast<ConfiguredParameter*>(parList.getSharedParameter(name).get());
        if (!confPar)
          throw Exception("makeTreeNode: unknown ConfiguredParameter " + name);
        
        mapBr.emplace(parTree.getEdgeIndex(branch),
                      ConfiguredParameter::create(context, {confPar->dependency(0)}, bp));
      }
      
      return std::shared_ptr<ProcessTree>(new ProcessTree(context, parTree, mapBr, modelmap));
    }

    
  } //end of namespace dataflow

} //end of namespace bpp


#endif //_PHYLO_TREE_BRREF_H_

