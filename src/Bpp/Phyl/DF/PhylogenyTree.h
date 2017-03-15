// File: PhylogenyTree.h
// Authors:
//   Francois Gindraud (2017)
// Created: 15/03/2017

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

#ifndef _PHYLOGENY_TREE_H_
#define _PHYLOGENY_TREE_H_

#include <Bpp/Phyl/NewLikelihood/SubstitutionProcess.h>

#include <Bpp/Phyl/DF/DataFlowBaseClasses.h>
#include <Bpp/Phyl/DF/ForRange.h>

#include <memory>
#include <vector>

namespace bpp
{
  class PhylogenyTree
  {
    // Topology + Branch len
  public:
    using IndexType = std::size_t;

    class Node
    {
    private:
      IndexType id;

    public:
      virtual ~Node() = default;
    };
    class Branch
    {
    private:
    public:
      virtual ~Branch() = default;
    };

  private:
    std::vector<std::unique_ptr<Node>> treeNodes_;

  protected:
    virtual Node* createNode(void) { return new Node; }
    virtual Branch* createBranch(void) { return new Branch; }
  };

  class PhylogenyProcess : public virtual PhylogenyTree
  {
    // Add a model per branch
  public:
    class Node : public virtual PhylogenyTree::Node
    {
    };
    class Branch : public virtual PhylogenyTree::Branch
    {
    };

  protected:
    Node* createNode(void) override { return new Node; }
    Branch* createBranch(void) override { return new Branch; }
  };

} // end of namespace bpp.

#endif //_PHYLOGENY_TREE_H_
