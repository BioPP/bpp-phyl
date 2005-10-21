H_FILES = \
AbstractSubstitutionModel.h\
AbstractTreeLikelihood.h\
AbstractDiscreteRatesAcrossSitesTreeLikelihood.h\
AbstractHomogeneousTreeLikelihood.h\
AbstractAgglomerativeDistanceMethod.h\
AncestralStateReconstruction.h\
MarginalAncestralStateReconstruction.h\
AgglomerativeDistanceMethod.h\
PhylogeneticsApplicationTools.h\
DistanceEstimation.h\
DistanceMatrix.h\
DSO78.h\
PseudoNewtonOptimizer.h\
HKY85.h\
HomogeneousTreeLikelihood.h\
DRHomogeneousTreeLikelihood.h\
DRTreeLikelihoodTools.h\
IODistanceMatrix.h\
IOTree.h\
iotree\
JCnuc.h\
JCprot.h\
JTT92.h\
K80.h\
models\
MutationProcess.h\
Newick.h\
NewtonBrentMetaOptimizer.h\
NonHomogeneousTreeLikelihood.h\
NucleotideSubstitutionModel.h\
NeighborJoining.h\
NNISearchable.h\
OptimizationTools.h\
PatternTools.h\
PGMA.h\
PhylipDistanceMatrixFormat.h\
ProteinSubstitutionModel.h\
SequenceSimulator.h\
SiteSimulator.h\
PreciseSequenceSimulator.h\
HomogeneousSequenceSimulator.h\
SubstitutionModel.h\
T92.h\
TN93.h\
TreeExceptions.h\
Node.h\
NodeTemplate.h\
Tree.h\
TreeTemplate.h\
TreeLikelihood.h\
TreeTools.h\
DiscreteRatesAcrossSitesTreeLikelihood.h\
RASTools.h\
TreeParsimonyScore.h\
AbstractTreeParsimonyScore.h\
DRTreeParsimonyScore.h\
SequenceSimulationTools.h\
UserProteinSubstitutionModel.h\
TopologySearch.h\
NNITopologySearch.h


CPP_FILES = \
AbstractSubstitutionModel.cpp\
AbstractTreeLikelihood.cpp\
AbstractDiscreteRatesAcrossSitesTreeLikelihood.cpp\
AbstractHomogeneousTreeLikelihood.cpp\
AbstractAgglomerativeDistanceMethod.cpp\
PhylogeneticsApplicationTools.cpp\
DistanceEstimation.cpp\
DSO78.cpp\
PseudoNewtonOptimizer.cpp\
HKY85.cpp\
HomogeneousTreeLikelihood.cpp\
DRHomogeneousTreeLikelihood.cpp\
DRTreeLikelihoodTools.cpp\
JCnuc.cpp\
JCprot.cpp\
JTT92.cpp\
K80.cpp\
MarginalAncestralStateReconstruction.cpp\
MutationProcess.cpp\
NeighborJoining.cpp\
Newick.cpp\
NewtonBrentMetaOptimizer.cpp\
Node.cpp\
NonHomogeneousTreeLikelihood.cpp\
NucleotideSubstitutionModel.cpp\
OptimizationTools.cpp\
PatternTools.cpp\
PGMA.cpp\
PhylipDistanceMatrixFormat.cpp\
ProteinSubstitutionModel.cpp\
HomogeneousSequenceSimulator.cpp\
SubstitutionModel.cpp\
T92.cpp\
TN93.cpp\
TreeExceptions.cpp\
TreeLikelihood.cpp\
TreeTools.cpp\
RASTools.cpp\
AbstractTreeParsimonyScore.cpp\
DRTreeParsimonyScore.cpp\
SequenceSimulationTools.cpp\
UserProteinSubstitutionModel.cpp\
NNITopologySearch.cpp

CPP_EXTRA = \
__DSO78FrequenciesCode\
__DSO78ExchangeabilityCode\
__JTT92FrequenciesCode\
__JTT92ExchangeabilityCode

lib_LTLIBRARIES = libphyl.la
libphyl_la_SOURCES = $(H_FILES) $(CPP_FILES)
EXTRA_libphyl_la_SOURCES = $(CPP_EXTRA)
libphyl_la_LDFLAGS = -version-info 1:0:0

library_includedir=$(includedir)/Phyl
library_include_HEADERS = $(H_FILES)
