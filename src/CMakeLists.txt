# CMake script for Bio++ PhylLib
# Authors:
#   Sylvain Gaillard
#   Julien Dutheil
#   Francois Gindraud (2017)
# Created: 20/08/2009

# File list
set (CPP_FILES
  Bpp/NewPhyl/DataFlow.cpp
  Bpp/NewPhyl/DataFlowNumeric.cpp
  Bpp/NewPhyl/Debug.cpp
  Bpp/NewPhyl/ImportMaster.cpp
  Bpp/NewPhyl/ImportNewlik.cpp
  Bpp/NewPhyl/Likelihood.cpp
  Bpp/NewPhyl/Model.cpp
  Bpp/NewPhyl/Phylogeny.cpp
  Bpp/NewPhyl/Topology.cpp
  Bpp/Phyl/App/PhylogeneticsApplicationTools.cpp
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
  Bpp/Phyl/Io/BppOFrequenciesSetFormat.cpp
  Bpp/Phyl/Io/BppOMultiTreeReaderFormat.cpp
  Bpp/Phyl/Io/BppOMultiTreeWriterFormat.cpp
  Bpp/Phyl/Io/BppORateDistributionFormat.cpp
  Bpp/Phyl/Io/BppOSubstitutionModelFormat.cpp
  Bpp/Phyl/Io/BppOTreeReaderFormat.cpp
  Bpp/Phyl/Io/BppOTreeWriterFormat.cpp
  Bpp/Phyl/Io/IoDistanceMatrixFactory.cpp
  Bpp/Phyl/Io/IoFrequenciesSetFactory.cpp
  Bpp/Phyl/Io/IoPairedSiteLikelihoods.cpp
  Bpp/Phyl/Io/IoSubstitutionModelFactory.cpp
  Bpp/Phyl/Io/IoTreeFactory.cpp
  Bpp/Phyl/Io/Newick.cpp
  Bpp/Phyl/Io/NexusIoTree.cpp
  Bpp/Phyl/Io/Nhx.cpp
  Bpp/Phyl/Io/PhylipDistanceMatrixFormat.cpp
  Bpp/Phyl/Likelihood/AbstractDiscreteRatesAcrossSitesTreeLikelihood.cpp
  Bpp/Phyl/Likelihood/AbstractHomogeneousTreeLikelihood.cpp
  Bpp/Phyl/Likelihood/AbstractNonHomogeneousTreeLikelihood.cpp
  Bpp/Phyl/Likelihood/AbstractTreeLikelihood.cpp
  Bpp/Phyl/Likelihood/DRASDRTreeLikelihoodData.cpp
  Bpp/Phyl/Likelihood/DRASRTreeLikelihoodData.cpp
  Bpp/Phyl/Likelihood/DRHomogeneousMixedTreeLikelihood.cpp
  Bpp/Phyl/Likelihood/DRHomogeneousTreeLikelihood.cpp
  Bpp/Phyl/Likelihood/DRNonHomogeneousTreeLikelihood.cpp
  Bpp/Phyl/Likelihood/DRTreeLikelihoodTools.cpp
  Bpp/Phyl/Likelihood/GlobalClockTreeLikelihoodFunctionWrapper.cpp
  Bpp/Phyl/Likelihood/MarginalAncestralStateReconstruction.cpp
  Bpp/Phyl/Likelihood/NNIHomogeneousTreeLikelihood.cpp
  Bpp/Phyl/Likelihood/PairedSiteLikelihoods.cpp
  Bpp/Phyl/Likelihood/PseudoNewtonOptimizer.cpp
  Bpp/Phyl/Likelihood/RASTools.cpp
  Bpp/Phyl/Likelihood/RHomogeneousClockTreeLikelihood.cpp
  Bpp/Phyl/Likelihood/RHomogeneousMixedTreeLikelihood.cpp
  Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.cpp
  Bpp/Phyl/Likelihood/RNonHomogeneousMixedTreeLikelihood.cpp
  Bpp/Phyl/Likelihood/RNonHomogeneousTreeLikelihood.cpp
  Bpp/Phyl/Likelihood/TreeLikelihoodTools.cpp
  Bpp/Phyl/Mapping/DecompositionReward.cpp
  Bpp/Phyl/Mapping/DecompositionSubstitutionCount.cpp
  Bpp/Phyl/Mapping/LaplaceSubstitutionCount.cpp
  Bpp/Phyl/Mapping/NaiveSubstitutionCount.cpp
  Bpp/Phyl/Mapping/OneJumpSubstitutionCount.cpp
  Bpp/Phyl/Mapping/ProbabilisticRewardMapping.cpp
  Bpp/Phyl/Mapping/ProbabilisticSubstitutionMapping.cpp
  Bpp/Phyl/Mapping/RewardMappingTools.cpp
  Bpp/Phyl/Mapping/SubstitutionMappingTools.cpp
  Bpp/Phyl/Mapping/SubstitutionRegister.cpp
  Bpp/Phyl/Mapping/UniformizationSubstitutionCount.cpp
  Bpp/Phyl/Mapping/WeightedSubstitutionCount.cpp
  Bpp/Phyl/Model/AbstractBiblioMixedSubstitutionModel.cpp
  Bpp/Phyl/Model/AbstractBiblioSubstitutionModel.cpp
  Bpp/Phyl/Model/AbstractKroneckerWordSubstitutionModel.cpp
  Bpp/Phyl/Model/AbstractMixedSubstitutionModel.cpp
  Bpp/Phyl/Model/AbstractSubstitutionModel.cpp
  Bpp/Phyl/Model/AbstractWordSubstitutionModel.cpp
  Bpp/Phyl/Model/BinarySubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/AbstractCodonCpGSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/AbstractCodonDistanceSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/AbstractCodonFitnessSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/AbstractCodonFrequenciesSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/AbstractCodonPhaseFrequenciesSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/AbstractCodonSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/AbstractKroneckerCodonSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/CodonDistanceCpGSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/CodonDistanceFrequenciesSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/CodonDistancePhaseFrequenciesSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/CodonDistanceSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/CodonRateFrequenciesSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/CodonRateSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/GY94.cpp
  Bpp/Phyl/Model/Codon/KCM.cpp
  Bpp/Phyl/Model/Codon/KroneckerCodonDistanceFrequenciesSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/KroneckerCodonDistanceSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/MG94.cpp
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
  Bpp/Phyl/Model/FrequenciesSet/CodonFrequenciesSet.cpp
  Bpp/Phyl/Model/FrequenciesSet/FrequenciesSet.cpp
  Bpp/Phyl/Model/FrequenciesSet/MvaFrequenciesSet.cpp
  Bpp/Phyl/Model/FrequenciesSet/NucleotideFrequenciesSet.cpp
  Bpp/Phyl/Model/FrequenciesSet/WordFrequenciesSet.cpp
  Bpp/Phyl/Model/FromMixtureSubstitutionModel.cpp
  Bpp/Phyl/Model/KroneckerWordSubstitutionModel.cpp
  Bpp/Phyl/Model/MarkovModulatedSubstitutionModel.cpp
  Bpp/Phyl/Model/MixedSubstitutionModelSet.cpp
  Bpp/Phyl/Model/MixtureOfASubstitutionModel.cpp
  Bpp/Phyl/Model/MixtureOfSubstitutionModels.cpp
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
  Bpp/Phyl/Model/StateMap.cpp
  Bpp/Phyl/Model/SubstitutionModelSet.cpp
  Bpp/Phyl/Model/SubstitutionModelSetTools.cpp
  Bpp/Phyl/Model/WordSubstitutionModel.cpp
  Bpp/Phyl/NewLikelihood/AbstractLikelihoodNode.cpp
  Bpp/Phyl/NewLikelihood/AbstractLikelihoodTree.cpp
  Bpp/Phyl/NewLikelihood/AbstractLikelihoodTreeCalculation.cpp
  Bpp/Phyl/NewLikelihood/AbstractSubstitutionProcess.cpp
  Bpp/Phyl/NewLikelihood/AutoCorrelationSequenceEvolution.cpp
  Bpp/Phyl/NewLikelihood/ComputingTree.cpp
  Bpp/Phyl/NewLikelihood/HmmProcessEmissionProbabilities.cpp
  Bpp/Phyl/NewLikelihood/HmmSequenceEvolution.cpp
  Bpp/Phyl/NewLikelihood/MarginalAncestralReconstruction.cpp
  Bpp/Phyl/NewLikelihood/MixtureSequenceEvolution.cpp
  Bpp/Phyl/NewLikelihood/MultiProcessSequenceEvolution.cpp
  Bpp/Phyl/NewLikelihood/NonHomogeneousSubstitutionProcess.cpp
  Bpp/Phyl/NewLikelihood/OneProcessSequenceEvolution.cpp
  Bpp/Phyl/NewLikelihood/ParametrizablePhyloTree.cpp
  Bpp/Phyl/NewLikelihood/ParametrizableTree.cpp
  Bpp/Phyl/NewLikelihood/PartitionSequenceEvolution.cpp
  Bpp/Phyl/NewLikelihood/PhyloLikelihoods/AbstractPhyloLikelihood.cpp
  Bpp/Phyl/NewLikelihood/PhyloLikelihoods/AutoCorrelationOfAlignedPhyloLikelihood.cpp
  Bpp/Phyl/NewLikelihood/PhyloLikelihoods/AutoCorrelationProcessPhyloLikelihood.cpp
  Bpp/Phyl/NewLikelihood/PhyloLikelihoods/FormulaOfPhyloLikelihood.cpp
  Bpp/Phyl/NewLikelihood/PhyloLikelihoods/HmmOfAlignedPhyloLikelihood.cpp
  Bpp/Phyl/NewLikelihood/PhyloLikelihoods/HmmPhyloEmissionProbabilities.cpp
  Bpp/Phyl/NewLikelihood/PhyloLikelihoods/HmmProcessPhyloLikelihood.cpp
  Bpp/Phyl/NewLikelihood/PhyloLikelihoods/MixtureOfAlignedPhyloLikelihood.cpp
  Bpp/Phyl/NewLikelihood/PhyloLikelihoods/MixtureProcessPhyloLikelihood.cpp
  Bpp/Phyl/NewLikelihood/PhyloLikelihoods/MultiProcessSequencePhyloLikelihood.cpp
  Bpp/Phyl/NewLikelihood/PhyloLikelihoods/OneProcessSequencePhyloLikelihood.cpp
  Bpp/Phyl/NewLikelihood/PhyloLikelihoods/PartitionProcessPhyloLikelihood.cpp
  Bpp/Phyl/NewLikelihood/PhyloLikelihoods/ProductOfAlignedPhyloLikelihood.cpp
  Bpp/Phyl/NewLikelihood/PhyloLikelihoods/SetOfAbstractPhyloLikelihood.cpp
  Bpp/Phyl/NewLikelihood/PhyloLikelihoods/SetOfAlignedPhyloLikelihood.cpp
  Bpp/Phyl/NewLikelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.cpp
  Bpp/Phyl/NewLikelihood/RateAcrossSitesSubstitutionProcess.cpp
  Bpp/Phyl/NewLikelihood/RecursiveLikelihoodTree.cpp
  Bpp/Phyl/NewLikelihood/RecursiveLikelihoodTreeCalculation.cpp
  Bpp/Phyl/NewLikelihood/SimpleSubstitutionProcess.cpp
  Bpp/Phyl/NewLikelihood/SpeciationComputingNode.cpp
  Bpp/Phyl/NewLikelihood/SubstitutionProcessCollection.cpp
  Bpp/Phyl/NewLikelihood/SubstitutionProcessCollectionMember.cpp
  Bpp/Phyl/OptimizationTools.cpp
  Bpp/Phyl/Parsimony/AbstractTreeParsimonyScore.cpp
  Bpp/Phyl/Parsimony/DRTreeParsimonyData.cpp
  Bpp/Phyl/Parsimony/DRTreeParsimonyScore.cpp
  Bpp/Phyl/PatternTools.cpp
  Bpp/Phyl/PhyloStatistics.cpp
  Bpp/Phyl/Simulation/EvolutionSequenceSimulator.cpp
  Bpp/Phyl/Simulation/MutationProcess.cpp
  Bpp/Phyl/Simulation/NonHomogeneousSequenceSimulator.cpp
  Bpp/Phyl/Simulation/SequenceSimulationTools.cpp
  Bpp/Phyl/Simulation/SubstitutionProcessSequenceSimulator.cpp
  Bpp/Phyl/SitePatterns.cpp
  Bpp/Phyl/Tree/AwareNode.cpp
  Bpp/Phyl/Tree/BipartitionList.cpp
  Bpp/Phyl/Tree/BipartitionTools.cpp
  Bpp/Phyl/Tree/NNITopologySearch.cpp
  Bpp/Phyl/Tree/Node.cpp
  Bpp/Phyl/Tree/PhyloBranch.cpp
  Bpp/Phyl/Tree/PhyloBranchParam.cpp
  Bpp/Phyl/Tree/PhyloNode.cpp
  Bpp/Phyl/Tree/PhyloTree.cpp
  Bpp/Phyl/Tree/PhyloTreeExceptions.cpp
  Bpp/Phyl/Tree/PhyloTreeTools.cpp
  Bpp/Phyl/Tree/TreeExceptions.cpp
  Bpp/Phyl/Tree/TreeTemplateTools.cpp
  Bpp/Phyl/Tree/TreeTools.cpp
  )

# Config file (currently not installed as it is only used in impl files).
# $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}> is due to this file.
configure_file (
  ${CMAKE_CURRENT_SOURCE_DIR}/Bpp/NewPhyl/Config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/Bpp/NewPhyl/Config.h
  )

# Build the static lib
add_library (${PROJECT_NAME}-static STATIC ${CPP_FILES})
target_include_directories (${PROJECT_NAME}-static PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
  )
set_target_properties (${PROJECT_NAME}-static PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
target_link_libraries (${PROJECT_NAME}-static ${BPP_LIBS_STATIC} Eigen3::Eigen)

# Build the shared lib
add_library (${PROJECT_NAME}-shared SHARED ${CPP_FILES})
target_include_directories (${PROJECT_NAME}-shared PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
  )
set_target_properties (${PROJECT_NAME}-shared
  PROPERTIES OUTPUT_NAME ${PROJECT_NAME}
  MACOSX_RPATH 1
  VERSION ${${PROJECT_NAME}_VERSION}
  SOVERSION ${${PROJECT_NAME}_VERSION_MAJOR}
  )
target_link_libraries (${PROJECT_NAME}-shared ${BPP_LIBS_SHARED} Eigen3::Eigen)

# Install libs and headers
install (
  TARGETS ${PROJECT_NAME}-static ${PROJECT_NAME}-shared
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
