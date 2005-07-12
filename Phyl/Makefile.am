H_FILES = \
AbstractSubstitutionModel.h\
AbstractTreeLikelihood.h\
AbstractDiscreteRatesAcrossSitesTreeLikelihood.h\
AbstractHomogeneousTreeLikelihood.h\
AbstractAgglomerativeDistanceMethod.h\
AncestralStateReconstruction.h\
MarginalAncestralStateReconstruction.h\
AgglomerativeDistanceMethod.h\
ApplicationTools.h\
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
OptimizationTools.h\
PatternTools.h\
PGMA.h\
PhylipDistanceMatrixFormat.h\
ProteinSubstitutionModel.h\
SequenceSimulator.h\
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
RASTools.h


CPP_FILES = \
AbstractSubstitutionModel.cpp\
AbstractTreeLikelihood.cpp\
AbstractDiscreteRatesAcrossSitesTreeLikelihood.cpp\
AbstractHomogeneousTreeLikelihood.cpp\
AbstractAgglomerativeDistanceMethod.cpp\
ApplicationTools.cpp\
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
NonHomogeneousTreeLikelihood.cpp\
NucleotideSubstitutionModel.cpp\
OptimizationTools.cpp\
PatternTools.cpp\
PGMA.cpp\
PhylipDistanceMatrixFormat.cpp\
ProteinSubstitutionModel.cpp\
SequenceSimulator.cpp\
SubstitutionModel.cpp\
T92.cpp\
TN93.cpp\
Tree.cpp\
TreeExceptions.cpp\
TreeLikelihood.cpp\
TreeTools.cpp\
RASTools.cpp

CPP_EXTRA = \
DSO78FrequenciesCode.cpp\
DSO78ExchangeabilityCode.cpp\
JTT92FrequenciesCode.cpp\
JTT92ExchangeabilityCode.cpp

lib_LTLIBRARIES = libphyl.la
libphyl_la_SOURCES = $(H_FILES) $(CPP_FILES)
EXTRA_libphyl_la_SOURCES = $(CPP_EXTRA)
libphyl_la_LDFLAGS = -version-info 1:0:0

library_includedir=$(includedir)/Phyl
library_include_HEADERS = $(H_FILES)
