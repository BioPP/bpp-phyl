# CMake script for Bio++ PhylLib
# Authors:
#   Sylvain Gaillard
#   Julien Dutheil
#   Francois Gindraud (2017)
# Created: 20/08/2009

# File list
set (CPP_FILES
  Bpp/Phyl/App/PhylogeneticsApplicationTools.cpp
  Bpp/Phyl/BipartitionList.cpp
  Bpp/Phyl/BipartitionTools.cpp
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
  Bpp/Phyl/Io/BppOFrequencySetFormat.cpp
  Bpp/Phyl/Io/BppOMultiTreeReaderFormat.cpp
  Bpp/Phyl/Io/BppOMultiTreeWriterFormat.cpp
  Bpp/Phyl/Io/BppORateDistributionFormat.cpp
  Bpp/Phyl/Io/BppOSubstitutionModelFormat.cpp
  Bpp/Phyl/Io/BppOTransitionModelFormat.cpp
  Bpp/Phyl/Io/BppOTreeReaderFormat.cpp
  Bpp/Phyl/Io/BppOTreeWriterFormat.cpp
  Bpp/Phyl/Io/IoDistanceMatrixFactory.cpp
  Bpp/Phyl/Io/IoFrequencySetFactory.cpp
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
  Bpp/Phyl/Likelihood/JointLikelihoodFunction.cpp
  Bpp/Phyl/Mapping/DecompositionMethods.cpp
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
  Bpp/Phyl/Mapping/StochasticMapping.cpp
  Bpp/Phyl/Model/AbstractBiblioMixedTransitionModel.cpp
  Bpp/Phyl/Model/AbstractBiblioSubstitutionModel.cpp
  Bpp/Phyl/Model/AbstractFromSubstitutionModelTransitionModel.cpp
  Bpp/Phyl/Model/AbstractKroneckerWordSubstitutionModel.cpp
  Bpp/Phyl/Model/AbstractMixedTransitionModel.cpp
  Bpp/Phyl/Model/AbstractSubstitutionModel.cpp
  Bpp/Phyl/Model/AbstractWordSubstitutionModel.cpp
  Bpp/Phyl/Model/BinarySubstitutionModel.cpp
  Bpp/Phyl/Model/TwoParameterBinarySubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/AbstractCodonAARateSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/AbstractCodonAAFitnessSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/AbstractCodonClusterAASubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/AbstractCodonCpGSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/AbstractCodonBGCSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/AbstractCodonDistanceSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/AbstractCodonFitnessSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/AbstractCodonFrequenciesSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/AbstractCodonPhaseFrequenciesSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/AbstractCodonSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/AbstractKroneckerCodonSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/CodonAdHocSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/CodonDistanceFrequenciesSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/CodonDistancePhaseFrequenciesSubstitutionModel.cpp
  Bpp/Phyl/Model/Codon/CodonDistanceSubstitutionModel.cpp
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
  Bpp/Phyl/Model/Codon/RELAX.cpp
  Bpp/Phyl/Model/FrequencySet/CodonFrequencySet.cpp
  Bpp/Phyl/Model/FrequencySet/FrequencySet.cpp
  Bpp/Phyl/Model/FrequencySet/MvaFrequencySet.cpp
  Bpp/Phyl/Model/FrequencySet/NucleotideFrequencySet.cpp
  Bpp/Phyl/Model/FrequencySet/WordFrequencySet.cpp
  Bpp/Phyl/Model/FromMixtureSubstitutionModel.cpp
  Bpp/Phyl/Model/InMixedSubstitutionModel.cpp
  Bpp/Phyl/Model/KroneckerWordSubstitutionModel.cpp
  Bpp/Phyl/Model/MarkovModulatedSubstitutionModel.cpp
  Bpp/Phyl/Model/MixedSubstitutionModelSet.cpp
  Bpp/Phyl/Model/MixtureOfATransitionModel.cpp
  Bpp/Phyl/Model/MixtureOfTransitionModels.cpp
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
  Bpp/Phyl/Model/SubstitutionModelSet.cpp
  Bpp/Phyl/Model/SubstitutionModelSetTools.cpp
  Bpp/Phyl/Model/WordSubstitutionModel.cpp
  Bpp/Phyl/NNITopologySearch.cpp
  Bpp/Phyl/Node.cpp
  Bpp/Phyl/OptimizationTools.cpp
  Bpp/Phyl/Parsimony/AbstractTreeParsimonyScore.cpp
  Bpp/Phyl/Parsimony/DRTreeParsimonyData.cpp
  Bpp/Phyl/Parsimony/DRTreeParsimonyScore.cpp
  Bpp/Phyl/PatternTools.cpp
  Bpp/Phyl/PhyloStatistics.cpp
  Bpp/Phyl/Simulation/MutationProcess.cpp
  Bpp/Phyl/Simulation/NonHomogeneousSequenceSimulator.cpp
  Bpp/Phyl/Simulation/SequenceSimulationTools.cpp
  Bpp/Phyl/SitePatterns.cpp
  Bpp/Phyl/TreeExceptions.cpp
  Bpp/Phyl/TreeTemplateTools.cpp
  Bpp/Phyl/TreeTools.cpp 
  Bpp/Phyl/TreeIterator.cpp
  )

IF(BUILD_STATIC)
  # Build the static lib
  add_library (${PROJECT_NAME}-static STATIC ${CPP_FILES})
  target_include_directories (${PROJECT_NAME}-static PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
    )
  set_target_properties (${PROJECT_NAME}-static PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
  target_link_libraries (${PROJECT_NAME}-static ${BPP_LIBS_STATIC})
ENDIF()

# Build the shared lib
add_library (${PROJECT_NAME}-shared SHARED ${CPP_FILES})
target_include_directories (${PROJECT_NAME}-shared PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
  )
set_target_properties (${PROJECT_NAME}-shared
  PROPERTIES OUTPUT_NAME ${PROJECT_NAME}
  MACOSX_RPATH 1
  VERSION ${${PROJECT_NAME}_VERSION}
  SOVERSION ${${PROJECT_NAME}_VERSION_MAJOR}
  )
target_link_libraries (${PROJECT_NAME}-shared ${BPP_LIBS_SHARED})

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
