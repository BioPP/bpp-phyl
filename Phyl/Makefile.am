H_FILES = \
AbstractSubstitutionModel.h\
AbstractTreeLikelihood.h\
ApplicationTools.h\
HKY85.h\
HomogeneousTreeLikelihood.h\
IOTree.h\
JCnuc.h\
JTT.h\
K80.h\
MutationProcess.h\
Newick.h\
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
TreeTools.h


CPP_FILES = \
AbstractSubstitutionModel.cpp\
AbstractTreeLikelihood.cpp\
ApplicationTools.cpp\
HKY85.cpp\
HomogeneousTreeLikelihood.cpp\
IOTree.cpp\
JCnuc.cpp\
JTT.cpp\
K80.cpp\
MutationProcess.cpp\
Newick.cpp\
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
TreeTools.cpp


lib_LTLIBRARIES = libphyl.la
libphyl_la_SOURCES = $(H_FILES) $(CPP_FILES)
libphyl_la_LDFLAGS = -version-info 1:0:0

library_includedir=$(includedir)/Phyl
library_include_HEADERS = $(H_FILES)
