//
// File: NucleicSubstitutionModel.h
// Created by:  <@bogdanof>
// Created on: Tue May 27 11:03:53 2003
//

#include "NucleotideSubstitutionModel.h"

NucleotideSubstitutionModel::NucleotideSubstitutionModel(const Alphabet * alpha):
	AbstractSubstitutionModel(alpha) {}
	
NucleotideSubstitutionModel::~NucleotideSubstitutionModel() {}
