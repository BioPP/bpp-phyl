//
// File: SubstitutionMatrix.h
// Created by:  <Julien.Dutheil@univ-montp2.fr>
// Created on: Mon May 26 14:52:34 2003
//

#include "SubstitutionModel.h"

/******************************************************************************
 *                        SubstitutionModel exceptions:                       *
 ******************************************************************************/

SubstitutionModelException::SubstitutionModelException(const char *   text, const SubstitutionModel * sm) :
	Exception("SubstitutionModelException: " + string(text) + (sm != NULL ? "(" + sm -> getName() + ")" : "")),
	sm(sm) {};
SubstitutionModelException::SubstitutionModelException(const string & text, const SubstitutionModel * sm) :
	Exception("SubstitutionModelException: " + text + (sm != NULL ? "(" + sm -> getName() + ")" : "")),
	sm(sm) {};
SubstitutionModelException::~SubstitutionModelException() throw() {};
const SubstitutionModel * SubstitutionModelException::getSubstitutionModel() const { return sm; }
