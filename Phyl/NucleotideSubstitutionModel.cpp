//
// File: NucleicSubstitutionModel.h
// Created by:  <@bogdanof>
// Created on: Tue May 27 11:03:53 2003
//

#include "NucleotideSubstitutionModel.h"

NucleotideSubstitutionModel::NucleotideSubstitutionModel(const NucleicAlphabet * alpha):
	AbstractSubstitutionModel(alpha) {}
	
NucleotideSubstitutionModel::~NucleotideSubstitutionModel() {}
