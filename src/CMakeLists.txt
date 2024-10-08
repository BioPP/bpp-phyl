# SPDX-FileCopyrightText: The Bio++ Development Group
#
# SPDX-License-Identifier: CECILL-2.1

# File list
set (CPP_FILES
  Bpp/Phyl/App/PhylogeneticsApplicationTools.cpp
  Bpp/Phyl/Legacy/App/PhylogeneticsApplicationTools.cpp
  Bpp/Phyl/App/BppPhylogeneticsApplication.cpp
  Bpp/Phyl/Distance/AbstractAgglomerativeDistanceMethod.cpp
  Bpp/Phyl/Distance/BioNJ.cpp
  Bpp/Phyl/Distance/DistanceEstimation.cpp
  Bpp/Phyl/Distance/HierarchicalClustering.cpp
  Bpp/Phyl/Distance/NeighborJoining.cpp
  Bpp/Phyl/Distance/PGMA.cpp
  Bpp/Phyl/Graphics/AbstractDendrogramPlot.cpp
  Bpp/Phyl/Graphics/AbstractTreeDrawing.cpp
  Bpp/Phyl/Graphics/CladogramPlot.cpp
  Bpp/Phyl/Graphics/PhylogramPlot.cpp
  Bpp/Phyl/Graphics/TreeDrawingDisplayControler.cpp
  Bpp/Phyl/Graphics/TreeDrawingListener.cpp
  Bpp/Phyl/Io/BppOBranchModelFormat.cpp
  Bpp/Phyl/Io/BppOFrequencySetFormat.cpp
  Bpp/Phyl/Io/BppOMultiTreeReaderFormat.cpp
  Bpp/Phyl/Io/BppOMultiTreeWriterFormat.cpp
  Bpp/Phyl/Io/BppORateDistributionFormat.cpp
  Bpp/Phyl/Io/BppOSubstitutionModelFormat.cpp
  Bpp/Phyl/Io/BppOTransitionModelFormat.cpp
  Bpp/Phyl/Io/BppOTreeReaderFormat.cpp
  Bpp/Phyl/Io/BppOTreeWriterFormat.cpp
  Bpp/Phyl/Io/IoFrequencySetFactory.cpp
  Bpp/Phyl/Legacy/Io/IoPairedSiteLikelihoods.cpp
  Bpp/Phyl/Io/IoSubstitutionModelFactory.cpp
  Bpp/Phyl/Io/IoTreeFactory.cpp
  Bpp/Phyl/Io/ExtendedNewick.cpp
  Bpp/Phyl/Io/Newick.cpp
  Bpp/Phyl/Io/NexusIoTree.cpp
  Bpp/Phyl/Io/Nhx.cpp
  Bpp/Phyl/Legacy/Likelihood/AbstractDiscreteRatesAcrossSitesTreeLikelihood.cpp
  Bpp/Phyl/Legacy/Likelihood/AbstractHomogeneousTreeLikelihood.cpp
  Bpp/Phyl/Legacy/Likelihood/AbstractNonHomogeneousTreeLikelihood.cpp
  Bpp/Phyl/Legacy/Likelihood/AbstractTreeLikelihood.cpp
  Bpp/Phyl/Legacy/Likelihood/DRASDRTreeLikelihoodData.cpp
  Bpp/Phyl/Legacy/Likelihood/DRASRTreeLikelihoodData.cpp
  Bpp/Phyl/Legacy/Likelihood/DRHomogeneousMixedTreeLikelihood.cpp
  Bpp/Phyl/Legacy/Likelihood/DRHomogeneousTreeLikelihood.cpp
  Bpp/Phyl/Legacy/Likelihood/DRNonHomogeneousTreeLikelihood.cpp
  Bpp/Phyl/Legacy/Likelihood/DRTreeLikelihoodTools.cpp
  Bpp/Phyl/Legacy/Likelihood/GlobalClockTreeLikelihoodFunctionWrapper.cpp
  Bpp/Phyl/Legacy/Likelihood/MarginalAncestralStateReconstruction.cpp
  Bpp/Phyl/Legacy/Likelihood/NNIHomogeneousTreeLikelihood.cpp
  Bpp/Phyl/Legacy/Likelihood/PairedSiteLikelihoods.cpp
  Bpp/Phyl/Legacy/Likelihood/RASTools.cpp
  Bpp/Phyl/Legacy/Likelihood/RHomogeneousMixedTreeLikelihood.cpp
  Bpp/Phyl/Legacy/Likelihood/RHomogeneousTreeLikelihood.cpp
  Bpp/Phyl/Legacy/Likelihood/RNonHomogeneousMixedTreeLikelihood.cpp
  Bpp/Phyl/Legacy/Likelihood/RNonHomogeneousTreeLikelihood.cpp
  Bpp/Phyl/Legacy/Likelihood/TreeLikelihoodTools.cpp
  Bpp/Phyl/Legacy/Mapping/ProbabilisticRewardMapping.cpp 
  Bpp/Phyl/Legacy/Mapping/ProbabilisticSubstitutionMapping.cpp
  Bpp/Phyl/Legacy/Mapping/SubstitutionMappingTools.cpp
  Bpp/Phyl/Legacy/Mapping/RewardMappingTools.cpp 
  Bpp/Phyl/Legacy/Simulation/NonHomogeneousSequenceSimulator.cpp
  Bpp/Phyl/PseudoNewtonOptimizer.cpp
  Bpp/Phyl/Likelihood/DataFlow/BackwardLikelihoodTree.cpp
  Bpp/Phyl/Likelihood/DataFlow/CollectionNodes.cpp
  Bpp/Phyl/Likelihood/DataFlow/DataFlow.cpp
  Bpp/Phyl/Likelihood/DataFlow/DataFlowNumeric.cpp
  Bpp/Phyl/Likelihood/DataFlow/DataFlowCWise.cpp
  Bpp/Phyl/Likelihood/DataFlow/DataFlowCWiseComputing.cpp
  Bpp/Phyl/Likelihood/DataFlow/DiscreteDistribution.cpp
  Bpp/Phyl/Likelihood/DataFlow/ForwardLikelihoodTree.cpp
  Bpp/Phyl/Likelihood/DataFlow/ExtendedFloat.cpp
  Bpp/Phyl/Likelihood/DataFlow/FrequencySet.cpp
  Bpp/Phyl/Likelihood/DataFlow/LikelihoodCalculationOnABranch.cpp
  Bpp/Phyl/Likelihood/DataFlow/LikelihoodCalculationSingleProcess.cpp
  Bpp/Phyl/Likelihood/DataFlow/Model.cpp
  Bpp/Phyl/Likelihood/DataFlow/Parameter.cpp
  Bpp/Phyl/Likelihood/DataFlow/Parametrizable.cpp
  Bpp/Phyl/Likelihood/DataFlow/ProcessTree.cpp
  Bpp/Phyl/Likelihood/DataFlow/Simplex_DF.cpp
  Bpp/Phyl/Likelihood/DataFlow/TransitionMatrix.cpp
  Bpp/Phyl/Likelihood/ModelPath.cpp
  Bpp/Phyl/Likelihood/ModelScenario.cpp
  Bpp/Phyl/Likelihood/ProcessComputationTree.cpp
  Bpp/Phyl/Likelihood/AbstractSubstitutionProcess.cpp
  Bpp/Phyl/Likelihood/AbstractAutonomousSubstitutionProcess.cpp
  Bpp/Phyl/Likelihood/SimpleSubstitutionProcess.cpp
  Bpp/Phyl/Likelihood/NonHomogeneousSubstitutionProcess.cpp
  Bpp/Phyl/Likelihood/RateAcrossSitesSubstitutionProcess.cpp
  Bpp/Phyl/Likelihood/SubstitutionProcessCollection.cpp
  Bpp/Phyl/Likelihood/SubstitutionProcessCollectionMember.cpp
  Bpp/Phyl/Likelihood/OneProcessSequenceEvolution.cpp
  Bpp/Phyl/Likelihood/MultiProcessSequenceEvolution.cpp
  Bpp/Phyl/Likelihood/PartitionSequenceEvolution.cpp
  Bpp/Phyl/Likelihood/MixtureSequenceEvolution.cpp
  Bpp/Phyl/Likelihood/HmmSequenceEvolution.cpp
  Bpp/Phyl/Likelihood/AutoCorrelationSequenceEvolution.cpp
  Bpp/Phyl/Likelihood/ParametrizablePhyloTree.cpp
  Bpp/Phyl/Likelihood/PhyloLikelihoods/AlignedPhyloLikelihoodAutoCorrelation.cpp
  Bpp/Phyl/Likelihood/PhyloLikelihoods/AlignedPhyloLikelihoodHmm.cpp
  Bpp/Phyl/Likelihood/PhyloLikelihoods/AlignedPhyloLikelihoodMixture.cpp
  Bpp/Phyl/Likelihood/PhyloLikelihoods/AlignedPhyloLikelihoodProduct.cpp
  Bpp/Phyl/Likelihood/PhyloLikelihoods/AlignedPhyloLikelihoodSet.cpp
  Bpp/Phyl/Likelihood/PhyloLikelihoods/AutoCorrelationProcessPhyloLikelihood.cpp
  Bpp/Phyl/Likelihood/PhyloLikelihoods/HmmLikelihoodComputation.cpp
  Bpp/Phyl/Likelihood/PhyloLikelihoods/HmmLikelihood_DF.cpp
  Bpp/Phyl/Likelihood/PhyloLikelihoods/HmmPhyloEmissionProbabilities.cpp
  Bpp/Phyl/Likelihood/PhyloLikelihoods/HmmProcessPhyloLikelihood.cpp
  Bpp/Phyl/Likelihood/PhyloLikelihoods/MixtureProcessPhyloLikelihood.cpp
  Bpp/Phyl/Likelihood/PhyloLikelihoods/MultiProcessSequencePhyloLikelihood.cpp
  Bpp/Phyl/Likelihood/PhyloLikelihoods/OnABranchPhyloLikelihood.cpp
  Bpp/Phyl/Likelihood/PhyloLikelihoods/OneProcessSequencePhyloLikelihood.cpp
  Bpp/Phyl/Likelihood/PhyloLikelihoods/PartitionProcessPhyloLikelihood.cpp
  Bpp/Phyl/Likelihood/PhyloLikelihoods/PhyloLikelihoodFormula.cpp
  Bpp/Phyl/Likelihood/PhyloLikelihoods/PhyloLikelihoodSet.cpp
  Bpp/Phyl/Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.cpp
  Bpp/Phyl/Likelihood/MarginalAncestralReconstruction.cpp
  Bpp/Phyl/Mapping/DecompositionMethods.cpp
  Bpp/Phyl/Mapping/DecompositionReward.cpp
  Bpp/Phyl/Mapping/DecompositionSubstitutionCount.cpp
  Bpp/Phyl/Mapping/LaplaceSubstitutionCount.cpp
  Bpp/Phyl/Mapping/NaiveSubstitutionCount.cpp
  Bpp/Phyl/Mapping/OneJumpSubstitutionCount.cpp
  Bpp/Phyl/Mapping/PhyloMappings/AbstractSinglePhyloSubstitutionMapping.cpp
  Bpp/Phyl/Mapping/PhyloMappings/OneProcessSequenceSubstitutionMapping.cpp
  Bpp/Phyl/Mapping/PhyloMappings/SingleProcessSubstitutionMapping.cpp
  Bpp/Phyl/Mapping/ProbabilisticRewardMapping.cpp 
  Bpp/Phyl/Mapping/ProbabilisticSubstitutionMapping.cpp
  Bpp/Phyl/Mapping/RewardMappingTools.cpp
  Bpp/Phyl/Mapping/SubstitutionDistance.cpp
  Bpp/Phyl/Mapping/SubstitutionMappingTools.cpp
  Bpp/Phyl/Mapping/SubstitutionRegister.cpp
  Bpp/Phyl/Mapping/UniformizationSubstitutionCount.cpp
  Bpp/Phyl/Mapping/WeightedSubstitutionCount.cpp
  Bpp/Phyl/Mapping/StochasticMapping.cpp
  Bpp/Phyl/Model/AbstractBiblioMixedTransitionModel.cpp
  Bpp/Phyl/Model/AbstractBiblioSubstitutionModel.cpp
  Bpp/Phyl/Model/AbstractFromSubstitutionModelTransitionModel.cpp
  Bpp/Phyl/Model/AbstractKroneckerWordSubstitutionModel.cpp
  Bpp/Phyl/Model/AbstractMixedTransitionModel.cpp
  Bpp/Phyl/Model/AbstractSubstitutionModel.cpp
  Bpp/Phyl/Model/AbstractWordSubstitutionModel.cpp
  Bpp/Phyl/Model/BinarySubstitutionModel.cpp
  Bpp/Phyl/Model/D1WalkSubstitutionModel.cpp
  Bpp/Phyl/Model/TwoParameterBinarySubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/AbstractCodonAARateSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/AbstractCodonAAFitnessSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/AbstractCodonBGCSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/AbstractCodonClusterAASubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/AbstractCodonCpGSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/AbstractCodonDistanceSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/AbstractCodonFitnessSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/AbstractCodonFrequenciesSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/AbstractCodonPhaseFrequenciesSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/AbstractCodonSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/AbstractDFPSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/AbstractKroneckerCodonSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/CodonAdHocSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/CodonDistanceFrequenciesSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/CodonDistancePhaseFrequenciesSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/CodonDistanceSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/DFPDistanceFrequenciesSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/DFP07.cpp
  Bpp/Phyl/Model/Codon/GY94.cpp
  Bpp/Phyl/Model/Codon/KCM.cpp
  Bpp/Phyl/Model/Codon/KroneckerCodonDistanceFrequenciesSubstitutionModel.cpp 
  Bpp/Phyl/Model/Codon/KroneckerCodonDistanceSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/MG94.cpp
  Bpp/Phyl/Model/Codon/RELAX.cpp
  Bpp/Phyl/Model/Codon/CodonSameAARateSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/SENCA.cpp
  Bpp/Phyl/Model/Codon/TripletSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/YN98.cpp
  Bpp/Phyl/Model/Codon/YNGP_M1.cpp
  Bpp/Phyl/Model/Codon/YNGP_M10.cpp
  Bpp/Phyl/Model/Codon/YNGP_M2.cpp
  Bpp/Phyl/Model/Codon/YNGP_M3.cpp
  Bpp/Phyl/Model/Codon/YNGP_M7.cpp
  Bpp/Phyl/Model/Codon/YNGP_M8.cpp
  Bpp/Phyl/Model/Codon/YNGP_M9.cpp
  Bpp/Phyl/Model/EquiprobableSubstitutionModel.cpp
  Bpp/Phyl/Model/FrequencySet/CodonFrequencySet.cpp
  Bpp/Phyl/Model/FrequencySet/FrequencySet.cpp
  Bpp/Phyl/Model/FrequencySet/MvaFrequencySet.cpp
  Bpp/Phyl/Model/FrequencySet/NucleotideFrequencySet.cpp
  Bpp/Phyl/Model/FrequencySet/WordFrequencySet.cpp
  Bpp/Phyl/Model/FromMixtureSubstitutionModel.cpp
  Bpp/Phyl/Model/InMixedSubstitutionModel.cpp
  Bpp/Phyl/Model/KroneckerWordSubstitutionModel.cpp
  Bpp/Phyl/Model/MarkovModulatedSubstitutionModel.cpp
  Bpp/Phyl/Legacy/Model/MixedSubstitutionModelSet.cpp
  Bpp/Phyl/Model/MixtureOfATransitionModel.cpp
  Bpp/Phyl/Model/MixtureOfTransitionModels.cpp
  Bpp/Phyl/Model/MultinomialFromTransitionModel.cpp
  Bpp/Phyl/Model/Nucleotide/F81.cpp
  Bpp/Phyl/Model/Nucleotide/F84.cpp
  Bpp/Phyl/Model/Nucleotide/GTR.cpp
  Bpp/Phyl/Model/Nucleotide/HKY85.cpp
  Bpp/Phyl/Model/Nucleotide/JCnuc.cpp
  Bpp/Phyl/Model/Nucleotide/K80.cpp
  Bpp/Phyl/Model/Nucleotide/L95.cpp
  Bpp/Phyl/Model/Nucleotide/RN95.cpp
  Bpp/Phyl/Model/Nucleotide/RN95s.cpp
  Bpp/Phyl/Model/Nucleotide/SSR.cpp
  Bpp/Phyl/Model/Nucleotide/T92.cpp
  Bpp/Phyl/Model/Nucleotide/TN93.cpp
  Bpp/Phyl/Model/Nucleotide/YpR.cpp
  Bpp/Phyl/Model/Nucleotide/gBGC.cpp
  Bpp/Phyl/Model/OneChangeTransitionModel.cpp
  Bpp/Phyl/Model/OneChangeRegisterTransitionModel.cpp
  Bpp/Phyl/Model/POMO.cpp
  Bpp/Phyl/Model/Protein/Coala.cpp
  Bpp/Phyl/Model/Protein/CoalaCore.cpp
  Bpp/Phyl/Model/Protein/DSO78.cpp
  Bpp/Phyl/Model/Protein/JCprot.cpp
  Bpp/Phyl/Model/Protein/JTT92.cpp
  Bpp/Phyl/Model/Protein/LG08.cpp
  Bpp/Phyl/Model/Protein/LG10_EX_EHO.cpp
  Bpp/Phyl/Model/Protein/LGL08_CAT.cpp
  Bpp/Phyl/Model/Protein/LLG08_EHO.cpp
  Bpp/Phyl/Model/Protein/LLG08_EX2.cpp
  Bpp/Phyl/Model/Protein/LLG08_EX3.cpp
  Bpp/Phyl/Model/Protein/LLG08_UL2.cpp
  Bpp/Phyl/Model/Protein/LLG08_UL3.cpp
  Bpp/Phyl/Model/Protein/UserProteinSubstitutionModel.cpp
  Bpp/Phyl/Model/Protein/WAG01.cpp
  Bpp/Phyl/Model/RE08.cpp
  Bpp/Phyl/Model/RegisterRatesSubstitutionModel.cpp
  Bpp/Phyl/Model/StateMap.cpp
  Bpp/Phyl/Legacy/Model/SubstitutionModelSet.cpp
  Bpp/Phyl/Legacy/Model/SubstitutionModelSetTools.cpp
  Bpp/Phyl/Model/TransitionFromTransitionModel.cpp
  Bpp/Phyl/Model/WordSubstitutionModel.cpp
  Bpp/Phyl/OptimizationTools.cpp
  Bpp/Phyl/Legacy/OptimizationTools.cpp
  Bpp/Phyl/Parsimony/AbstractTreeParsimonyScore.cpp
  Bpp/Phyl/Parsimony/DRTreeParsimonyData.cpp
  Bpp/Phyl/Parsimony/DRTreeParsimonyScore.cpp
  Bpp/Phyl/PatternTools.cpp
  Bpp/Phyl/PhyloStatistics.cpp
  Bpp/Phyl/Simulation/MutationProcess.cpp
  Bpp/Phyl/Simulation/EvolutionSequenceSimulator.cpp
  Bpp/Phyl/Simulation/GivenDataSubstitutionProcessSequenceSimulator.cpp
  Bpp/Phyl/Simulation/GivenDataSubstitutionProcessSiteSimulator.cpp
  Bpp/Phyl/Simulation/SequenceSimulationTools.cpp
  Bpp/Phyl/Simulation/SimpleSubstitutionProcessSequenceSimulator.cpp
  Bpp/Phyl/Simulation/SimpleSubstitutionProcessSiteSimulator.cpp
  Bpp/Phyl/Simulation/SubstitutionProcessSequenceSimulator.cpp
  Bpp/Phyl/SitePatterns.cpp
  Bpp/Phyl/Tree/BipartitionList.cpp
  Bpp/Phyl/Tree/BipartitionTools.cpp
  Bpp/Phyl/Legacy/Tree/NNITopologySearch.cpp
  Bpp/Phyl/Tree/Node.cpp
  Bpp/Phyl/Legacy/Tree/AwareNode.cpp
  Bpp/Phyl/Tree/TreeExceptions.cpp
  Bpp/Phyl/Tree/TreeIterator.cpp
  Bpp/Phyl/Tree/TreeTemplateTools.cpp
  Bpp/Phyl/Tree/TreeTools.cpp
  Bpp/Phyl/Tree/PhyloDAG.cpp
  Bpp/Phyl/Tree/PhyloTree.cpp
  Bpp/Phyl/Tree/PhyloNode.cpp
  Bpp/Phyl/Tree/PhyloBranch.cpp
  Bpp/Phyl/Tree/PhyloBranchParam.cpp
  Bpp/Phyl/Tree/PhyloTreeTools.cpp
  Bpp/Phyl/Tree/PhyloTreeExceptions.cpp
  )

IF(BUILD_STATIC)
  # Build the static lib
  add_library (${PROJECT_NAME}-static STATIC ${CPP_FILES})
  target_include_directories (${PROJECT_NAME}-static PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )
  set_target_properties (${PROJECT_NAME}-static PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
  target_link_libraries (${PROJECT_NAME}-static ${BPP_LIBS_STATIC} Eigen3::Eigen)
ENDIF()

# Build the shared lib
add_library (${PROJECT_NAME}-shared SHARED ${CPP_FILES})
target_include_directories (${PROJECT_NAME}-shared PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )
set_target_properties (${PROJECT_NAME}-shared
  PROPERTIES OUTPUT_NAME ${PROJECT_NAME}
  MACOSX_RPATH 1
  VERSION ${${PROJECT_NAME}_VERSION}
  SOVERSION ${${PROJECT_NAME}_VERSION_MAJOR}
  )
target_link_libraries (${PROJECT_NAME}-shared ${BPP_LIBS_SHARED} Eigen3::Eigen)

# Install libs and headers
IF(BUILD_STATIC)
install (
  TARGETS ${PROJECT_NAME}-static
  EXPORT ${PROJECT_NAME}-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )
ENDIF()

install (
  TARGETS ${PROJECT_NAME}-shared
  EXPORT ${PROJECT_NAME}-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )
install (
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Bpp
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING PATTERN "*.h"
  )
# Install cmake file representing targets
install (EXPORT ${PROJECT_NAME}-targets DESTINATION ${cmake-package-location})
