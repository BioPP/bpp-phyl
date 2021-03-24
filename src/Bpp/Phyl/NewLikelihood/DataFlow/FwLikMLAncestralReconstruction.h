
#ifndef _FWLIKMLANCESTRALRECONSTRUCTION_H_
#define _FWLIKMLANCESTRALRECONSTRUCTION_H_

#include "ForwardLikelihoodTree.h"
#include "Sequence_DF.h"

namespace bpp{
    using RootConditionalLikelihoodsRef = ValueRef<RowLik>;

    class FwLikMLAncestralReconstruction: public AssociationDAGlobalGraphObserver<ConditionalLikelihoodForward,ForwardLikelihoodBelow>{
        using DAClass = AssociationDAGlobalGraphObserver<ConditionalLikelihoodForward,ForwardLikelihoodBelow>;
        using MaxJointLik = MatrixMaxProduct<MatrixLik, MatrixLik, MatrixLik>;
        using LikelihoodRootConditional = MatrixMaxProduct<RowLik, RowLik, MatrixLik>;
        // This is a bit redundant... Myabe it is better to declare this class as inheriting from ForwardLikelihoodTree.
        // If so, these data members should be 'protected' not 'private' in ForwardLikelihoodTree. Should check if it is possible.
        protected:
            Context& context_;
            std::shared_ptr<ProcessTree> processTree_;
            MatrixDimension likelihoodMatrixDim_;
            const StateMap& statemap_;
            ValueRef<RowLik> rFreqs_;
            Eigen::Index nbState_;
            Eigen::Index nbSites_;
            std::map<Speciesindex, DAGindexes> mapNodesIndexes_; // For nodes that bring
            //information (ie not the empty ones)

            //std::map<Speciesindex, DAGindexes> mapEdgesIndexes_; // For edges that bring
            //information (ie not the empty ones)


        public:
            FwLikMLAncestralReconstruction(Context& c, 
                          std::shared_ptr<ProcessTree> tree,
                          const StateMap& statemap, ValueRef<RowLik> rootFreqs) : DAClass(), context_(c), processTree_(tree), likelihoodMatrixDim_(),
                                                         statemap_(statemap), rFreqs_(rootFreqs), nbState_(Eigen::Index(statemap.getNumberOfModelStates())),
                                                        nbSites_(0){
            }
            ConditionalLikelihoodForwardRef makeForwardComputationAtRoot(std::shared_ptr<ProcessNode> node, const AlignedValuesContainer & sites);


            ValueRef<double> initialize(const AlignedValuesContainer& sites);

        protected:

            ConditionalLikelihoodForwardRef makeForwardLikelihoodAtNode (std::shared_ptr<ProcessNode> node, const AlignedValuesContainer & sites);

            /*
            * @brief Compute ConditionalLikelihood for leaf.
            *
            */

            ConditionalLikelihoodForwardRef makeInitialConditionalLikelihood (const std::string & sequenceName, const AlignedValuesContainer & sites);


            void linkNodes(ConditionalLikelihoodForwardRef forwardNode, std::vector<std::shared_ptr<ProcessEdge>> childBranches, std::vector<ConditionalLikelihoodForwardRef> depE);
            void setNodeIndices(ConditionalLikelihoodForwardRef forwardNode, std::shared_ptr<ProcessNode> processNode, unsigned int spIndex);




    };


}

#endif // _FWLIKMLANCESTRALRECONSTRUCTION_H