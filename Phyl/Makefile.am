H_FILES = \
AbstractSubstitutionModel.h\
AbstractTreeLikelihood.h\
ApplicationTools.h\
DSO78.h\
GaltierNewtonOptimizer.h\
HKY85.h\
HomogeneousTreeLikelihood.h\
IOTree.h\
JCnuc.h\
JTT92.h\
K80.h\
MutationProcess.h\
Newick.h\
NewtonBrentMetaOptimizer.h\
NonHomogeneousTreeLikelihood.h\
NucleotideSubstitutionModel.h\
OptimizationTools.h\
PatternTools.h\
ProteinSubstitutionModel.h\
SequenceSimulator.h\
SubstitutionModel.h\
T92.h\
TN93.h\
TreeExceptions.h\
Tree.h\
TreeLikelihood.h\
TreeTools.h\
DiscreteRatesAcrossSites.h\
RASTools.h


CPP_FILES = \
AbstractSubstitutionModel.cpp\
AbstractTreeLikelihood.cpp\
ApplicationTools.cpp\
DSO78.cpp\
GaltierNewtonOptimizer.cpp\
HKY85.cpp\
HomogeneousTreeLikelihood.cpp\
IOTree.cpp\
JCnuc.cpp\
JTT92.cpp\
K80.cpp\
MutationProcess.cpp\
Newick.cpp\
NewtonBrentMetaOptimizer.cpp\
NonHomogeneousTreeLikelihood.cpp\
NucleotideSubstitutionModel.cpp\
OptimizationTools.cpp\
PatternTools.cpp\
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
