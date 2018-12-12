//
// File: PhyloTree_BrRef.h
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

#ifndef _PHYLO_TREE_BRREF_H_
#define _PHYLO_TREE_BRREF_H_

#include <Bpp/Phyl/Tree/PhyloNode.h>
#include <Bpp/Phyl/Tree/PhyloTree.h>

#include <Bpp/Graph/AssociationTreeGraphImplObserver.h>

#include "DiscreteDistribution.h"
#include "Model.h"

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

    using NumericMutableMap = std::map<uint, std::shared_ptr<NumericMutable<double>>>;
    using ValueRefMap = std::map<uint, ValueRef<double>>;

    using ModelMap = std::map<uint, std::shared_ptr<ConfiguredModel>>;

    ValueRefMap createBrLenMap(Context & c, const PhyloTree& tree)
    {
      std::vector<std::shared_ptr<PhyloBranch> > vB=tree.getAllEdges();

      ValueRefMap map;

      for (auto& branch:vB)
      {
        map.emplace (tree.getEdgeIndex(branch),
                     NumericMutable<double>::create (c, branch->hasLength()?branch->getLength():0.1));
      }

      return map;
    }

    /** Helper: create a map with mutable dataflow nodes for each
     *  branch of the tree.
     *  The map is indexed by branch ids.
     */
    
    // Branch specific DataFlow objects
    struct BrRef
    {
      std::shared_ptr<ConfiguredModel> model_;
      ValueRef<double> brlen_;

    public:      
      std::shared_ptr<ConfiguredModel> getModel() const
      {
        return model_;
      }

      void setModel(std::shared_ptr<ConfiguredModel> model)
      {
        model_=model;
      }

      ValueRef<double> getBrLen()
      {
        return brlen_;
      }

      void setBrLen(ValueRef<double>&& brlen)
      {
        brlen_=brlen;
      }

    };
    
      
    class PhyloTree_BrRef : public AssociationTreeGlobalGraphObserver<PhyloNode,BrRef>
    {
    public:
      // Map to link branch Ids and ValueRef at attribution of NodeRef to branches
      
      using NumericMutableMap = std::map<uint, std::shared_ptr<NumericMutable<double>>>;
      
      // The Values are shared as they are in branches
      // PhyloTree_BrRef(const PhyloTree& tree, const ValueRefMap& vrefmap);
      
      // PhyloTree_BrRef(const PhyloTree_BrRef& pTree);
      
      // PhyloTree_BrRef& operator=(const PhyloTree_BrRef& pTree);
    
      // PhyloTree_BrRef* clone() const { return new PhyloTree_BrRef(*this); }
      
      PhyloTree_BrRef(const PhyloTree& tree, ValueRefMap&& vrefmap) :
        AssociationTreeGlobalGraphObserver<PhyloNode,BrRef>(tree.getGraph())
      {
        vector<uint> vNodesId=tree.getGraph()->getAllNodes();
        
        for (auto& index:vNodesId)
        {
          shared_ptr<PhyloNode> pn=tree.getNodeFromGraphid(index);
          
          associateNode(pn,index);
          setNodeIndex(pn,tree.getNodeIndex(pn));
        }
        
        // Ids of the branches in the graph may be different from the ids in the observer phyloTree
        vector<uint> vEdgesId=tree.getGraph()->getAllEdges();
        
        for (auto& index:vEdgesId)
        {
          // retrieve PhyloBranch id
          const std::shared_ptr<PhyloBranch> pb=tree.getEdgeFromGraphid(index);
          uint ids=tree.getEdgeIndex(pb);

          if (vrefmap.find(ids)!=vrefmap.end())
          {
            std::shared_ptr<BrRef> brref(new BrRef({0, std::move(vrefmap.at(ids))}));
            associateEdge(std::move(brref), index);
            setEdgeIndex(getEdgeFromGraphid(index),ids);
          }
          else
            throw Exception("PhyloTree_BrRef::PhyloTree_BrRef missing reference for branch " + TextTools::toString(ids));
        }
      }

      PhyloTree_BrRef(const PhyloTree& tree, const ValueRefMap& vrefmap) :
        AssociationTreeGlobalGraphObserver<PhyloNode,BrRef>(tree.getGraph())
      {
        vector<uint> vNodesId=tree.getGraph()->getAllNodes();
        
        for (auto& index:vNodesId)
        {
          shared_ptr<PhyloNode> pn=tree.getNodeFromGraphid(index);
          
          associateNode(pn,index);
          setNodeIndex(pn,tree.getNodeIndex(pn));
        }
        
        rootAt(tree.getRoot());

        // Ids of the branches in the graph may be different from the ids in the observer phyloTree
        vector<uint> vEdgesId=tree.getGraph()->getAllEdges();
        
        for (auto& index:vEdgesId)
        {
          // retrieve PhyloBranch id
          const std::shared_ptr<PhyloBranch> pb=tree.getEdgeFromGraphid(index);
          uint ids=tree.getEdgeIndex(pb);

          if (vrefmap.find(ids)!=vrefmap.end())
          {
            std::shared_ptr<BrRef> brref(new BrRef({0, vrefmap.at(ids)}));
            associateEdge(std::move(brref), index);
            setEdgeIndex(getEdgeFromGraphid(index),ids);
          }
          else
            throw Exception("PhyloTree_BrRef::PhyloTree_BrRef missing reference for branch " + TextTools::toString(ids));
        }
      }

      PhyloTree_BrRef(const PhyloTree_BrRef& tree, ValueRefMap&& vrefmap) :
        AssociationTreeGlobalGraphObserver<PhyloNode,BrRef>(tree)
      {
        auto aEit=allEdgesIterator();
        
        while (!aEit->end())
        {
          auto edge=**aEit;
          uint ids=getEdgeIndex(edge);

          if (vrefmap.find(ids)!=vrefmap.end())
          {
            edge->setBrLen(std::move(vrefmap.at(ids)));
          }
          else
            throw Exception("PhyloTree_BrRef::PhyloTree_BrRef missing length on branch " + TextTools::toString(ids));

          aEit->next();
        }
      }

      PhyloTree_BrRef* clone() const { return 0;}//new PhyloTree_BrRef(*this); }

      PhyloTree_BrRef(const PhyloTree_BrRef& pTree): 
        AssociationTreeGlobalGraphObserver<PhyloNode,BrRef>(pTree.getGraph())
      {
      }
      
      PhyloTree_BrRef& operator=(const PhyloTree_BrRef& pTree)
      {
        //AssociationTreeGlobalGraphObserver<PhyloNode,Value<T>>::operator=(pTree);
        return *this;
      }

      void setBranchModels(ModelMap& modelmap)
      {
        auto aEit=allEdgesIterator();
        
        while (!aEit->end())
        {
          auto edge=**aEit;
          // retrieve PhyloBranch id
          uint id=getEdgeIndex(edge);

          if (modelmap.find(id)!=modelmap.end())
            edge->setModel(modelmap[id]);
          else
            throw Exception("PhyloTree_BrRef::setBranchModels missing model on branch " + TextTools::toString(id));

          aEit->next();
        }
      }
      
    };

    ValueRefMap createBrLenMap(Context & c, const PhyloTree_BrRef& tree)
    {
      // Ids of the branches
      
      ValueRefMap map;

      auto aEit=tree.allEdgesIterator();
        
      while (!aEit->end())
      {
        auto edge=**aEit;
        map.emplace (tree.getEdgeIndex(edge),
                     NumericMutable<double>::create (c, edge->getBrLen()->getValue()));
        aEit->next();
      }
      
      return map;
    }

    
    ValueRefMap multiplyBrLenMap(Context & c, const PhyloTree_BrRef& tree, std::shared_ptr<CategoryFromDiscreteDistribution>&& rate)
    {
      // Ids of the branches

      ValueRefMap map;

      auto aEit=tree.allEdgesIterator();
        
      while (!aEit->end())
      {
        auto edge=**aEit;

        map.emplace (tree.getEdgeIndex(edge),
                     CWiseMul<double, std::tuple<double, double>>::create (c, {edge->getBrLen(), rate}, Dimension<double>()));

        aEit->next();
      }
      return map;
    }

  } //end of namespace dataflow

} //end of namespace bpp


#endif //_PHYLO_TREE_BRREF_H_

