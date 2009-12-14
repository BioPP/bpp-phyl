# CMake script for Bio++ NumCalc
# Author: Sylvain Gaillard
# Created: 20/08/2009

# File list
set(CPP_FILES
  AbstractSubstitutionModel.cpp
  AbstractTreeLikelihood.cpp
  AbstractDiscreteRatesAcrossSitesTreeLikelihood.cpp
  AbstractHomogeneousTreeLikelihood.cpp
  AbstractNonHomogeneousTreeLikelihood.cpp
  AbstractAgglomerativeDistanceMethod.cpp
  PhylogeneticsApplicationTools.cpp
  DistanceEstimation.cpp
  DSO78.cpp
  DSO78F.cpp
  PseudoNewtonOptimizer.cpp
  HKY85.cpp
  RHomogeneousTreeLikelihood.cpp
  RHomogeneousMixedTreeLikelihood.cpp
  RNonHomogeneousTreeLikelihood.cpp
  DRNonHomogeneousTreeLikelihood.cpp
  DRASRTreeLikelihoodData.cpp
  DRHomogeneousTreeLikelihood.cpp
  NNIHomogeneousTreeLikelihood.cpp
  RHomogeneousClockTreeLikelihood.cpp
  DRASDRTreeLikelihoodData.cpp
  TreeLikelihoodTools.cpp
  DRTreeLikelihoodTools.cpp
  JCnuc.cpp
  JCprot.cpp
  JCprotF.cpp
  JTT92.cpp
  JTT92F.cpp
  K80.cpp
  MarginalAncestralStateReconstruction.cpp
  MutationProcess.cpp
  NeighborJoining.cpp
  Newick.cpp
  Node.cpp
  OptimizationTools.cpp
  PatternTools.cpp
  PGMA.cpp
  PhylipDistanceMatrixFormat.cpp
  NonHomogeneousSequenceSimulator.cpp
  SubstitutionModel.cpp
  T92.cpp
  TN93.cpp
  TreeExceptions.cpp
  TreeTools.cpp
  TreeTemplateTools.cpp
  RASTools.cpp
  AbstractTreeParsimonyScore.cpp
  DRTreeParsimonyScore.cpp
  DRTreeParsimonyData.cpp
  SequenceSimulationTools.cpp
  UserProteinSubstitutionModel.cpp
  UserProteinSubstitutionModelF.cpp
  NNITopologySearch.cpp
  GTR.cpp
  SitePatterns.cpp
  ProbabilisticSubstitutionMapping.cpp
  AnalyticalSubstitutionCount.cpp
  SimpleSubstitutionCount.cpp
  OneJumpSubstitutionCount.cpp
  SubstitutionMappingTools.cpp
  BioNJ.cpp
  SubstitutionModelFactory.cpp
  RateDistributionFactory.cpp
  IOTreeFactory.cpp
  IODistanceMatrixFactory.cpp
  F84.cpp
  MarkovModulatedSubstitutionModel.cpp
  BipartitionList.cpp
  BipartitionTools.cpp
  SubstitutionModelSet.cpp
  SubstitutionModelSetTools.cpp
  FrequenciesSet.cpp
  AbstractWordReversibleSubstitutionModel.cpp
  WordReversibleSubstitutionModel.cpp
  TripletReversibleSubstitutionModel.cpp
  AbstractCodonReversibleSubstitutionModel.cpp
  CodonNeutralReversibleSubstitutionModel.cpp
  CodonAsynonymousReversibleSubstitutionModel.cpp
  CodonNeutralFrequenciesReversibleSubstitutionModel.cpp
  CodonAsynonymousFrequenciesReversibleSubstitutionModel.cpp
  AbstractCodonFrequenciesReversibleSubstitutionModel.cpp
  RE08.cpp
  L95.cpp
  NexusIOTree.cpp
  GY94.cpp
  MG94.cpp
  YN98.cpp
  AbstractTreeDrawing.cpp
  AbstractDendrogramPlot.cpp
  PhylogramPlot.cpp
  CladogramPlot.cpp
  PhyloStatistics.cpp
  MixedModel.cpp
  )
set(H_FILES
  AbstractSubstitutionModel.h
  AbstractTreeLikelihood.h
  AbstractDiscreteRatesAcrossSitesTreeLikelihood.h
  AbstractHomogeneousTreeLikelihood.h
  AbstractNonHomogeneousTreeLikelihood.h
  AbstractAgglomerativeDistanceMethod.h
  AncestralStateReconstruction.h
  MarginalAncestralStateReconstruction.h
  AgglomerativeDistanceMethod.h
  PhylogeneticsApplicationTools.h
  DistanceEstimation.h
  DSO78.h
  DSO78F.h
  PseudoNewtonOptimizer.h
  HKY85.h
  TreeLikelihood.h
  DRTreeLikelihood.h
  HomogeneousTreeLikelihood.h
  RHomogeneousTreeLikelihood.h
  RHomogeneousMixedTreeLikelihood.h
  DRHomogeneousTreeLikelihood.h
  NNIHomogeneousTreeLikelihood.h
  NonHomogeneousTreeLikelihood.h
  RNonHomogeneousTreeLikelihood.h
  DRNonHomogeneousTreeLikelihood.h
  treelikelihoods
  ClockTreeLikelihood.h
  RHomogeneousClockTreeLikelihood.h
  SitePartitionTreeLikelihood.h
  TreeLikelihoodTools.h
  DRTreeLikelihoodTools.h
  IODistanceMatrix.h
  iodistance
  IOTree.h
  iotree
  JCnuc.h
  JCprot.h
  JCprotF.h
  JTT92.h
  JTT92F.h
  K80.h
  models
  MutationProcess.h
  Newick.h
  NucleotideSubstitutionModel.h
  NeighborJoining.h
  NNISearchable.h
  OptimizationTools.h
  PatternTools.h
  PGMA.h
  PhylipDistanceMatrixFormat.h
  ProteinSubstitutionModel.h
  ProteinSubstitutionModelWithFrequencies.h
  SequenceSimulator.h
  SiteSimulator.h
  DetailedSiteSimulator.h
  NonHomogeneousSequenceSimulator.h
  HomogeneousSequenceSimulator.h
  SubstitutionModel.h
  FrequenciesSet.h
  T92.h
  TN93.h
  TreeExceptions.h
  Node.h
  NodeTemplate.h
  Tree.h
  TreeTemplate.h
  TreeLikelihoodData.h
  AbstractTreeLikelihoodData.h
  DRASRTreeLikelihoodData.h
  DRASDRTreeLikelihoodData.h
  TreeTools.h
  TreeTemplateTools.h
  DiscreteRatesAcrossSitesTreeLikelihood.h
  RASTools.h
  TreeParsimonyScore.h
  TreeParsimonyData.h
  AbstractTreeParsimonyScore.h
  AbstractTreeParsimonyData.h
  DRTreeParsimonyScore.h
  DRTreeParsimonyData.h
  SequenceSimulationTools.h
  UserProteinSubstitutionModel.h
  UserProteinSubstitutionModelF.h
  TopologySearch.h
  NNITopologySearch.h
  GTR.h
  trees
  SitePatterns.h
  SubstitutionMapping.h
  ProbabilisticSubstitutionMapping.h
  SubstitutionCount.h
  AnalyticalSubstitutionCount.h
  SimpleSubstitutionCount.h
  OneJumpSubstitutionCount.h
  IndexToCount.h
  SubstitutionMappingTools.h
  BioNJ.h
  SubstitutionModelFactory.h
  RateDistributionFactory.h
  IOTreeFactory.h
  IODistanceMatrixFactory.h
  F84.h
  MarkovModulatedSubstitutionModel.h
  G2001.h
  TS98.h
  BipartitionList.h
  BipartitionTools.h
  SubstitutionModelSet.h
  SubstitutionModelSetTools.h
  L95.h
  RE08.h
  AbstractWordReversibleSubstitutionModel.h
  WordReversibleSubstitutionModel.h
  TripletReversibleSubstitutionModel.h
  AbstractCodonReversibleSubstitutionModel.h
  CodonNeutralReversibleSubstitutionModel.h
  CodonAsynonymousReversibleSubstitutionModel.h
  AbstractCodonFrequenciesReversibleSubstitutionModel.h
  CodonNeutralFrequenciesReversibleSubstitutionModel.h
  CodonAsynonymousFrequenciesReversibleSubstitutionModel.h
  NexusIOTree.h
  GY94.h
  MG94.h
  YN98.h
  TreeDrawing.h
  AbstractTreeDrawing.h
  AbstractDendrogramPlot.h
  PhylogramPlot.h
  CladogramPlot.h
  PhyloStatistics.h
  MixedModel.h
  )

# Build the static lib
add_library(bppphyl-static STATIC ${CPP_FILES})
set_target_properties(bppphyl-static
  PROPERTIES OUTPUT_NAME bpp-phyl
  CLEAN_DIRECT_OUTPUT 1
  )
target_link_libraries(bppphyl-static ${LIBS})

# Build the shared lib
add_library(bppphyl-shared SHARED ${CPP_FILES})
set_target_properties(bppphyl-shared
  PROPERTIES OUTPUT_NAME bpp-phyl
  CLEAN_DIRECT_OUTPUT 1
  VERSION ${BPPPHYL_VERSION}
  SOVERSION ${BPPPHYL_VERSION_MAJOR}
  )
target_link_libraries(bppphyl-shared ${LIBS})

# Install libs
install(TARGETS bppphyl-static bppphyl-shared DESTINATION lib)

# Install headers
install(FILES ${H_FILES} DESTINATION include/Phyl)
