#ifndef _PHYLOTREE_HPP_
#define _PHYLOTREE_HPP_

#include <Bpp/Graph/AssociationTreeGraphObserver.h>

#include "Tree.h"

namespace bpp
{

class PhyloNode;
class PhyloBranch;


  /**
   * @brief Defines a Phylogenetic Tree based on a TreeGraph & its associationObserver
   *
   * @author Thomas Bigot
   */
  
    
  class PhyloTree:
    public SimpleAssociationTreeGraphObserver<PhyloNode,PhyloBranch,SimpleTreeGraph<SimpleGraph> >
    {
    public:

      PhyloTree(const Tree&);

      
    };
    
}

#endif
