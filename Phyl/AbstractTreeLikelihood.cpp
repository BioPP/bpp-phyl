//
// File: AbstractTreeLikelihood.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Fri Oct 17 17:57:21 2003
//

#include "AbstractTreeLikelihood.h"

/******************************************************************************/

AbstractTreeLikelihood::~AbstractTreeLikelihood() {}

/******************************************************************************/	
	
void AbstractTreeLikelihood::computeTreeLikelihood() {
	computeSubtreeLikelihood(_tree -> getRootNode());
}

/******************************************************************************/

const SiteContainer * AbstractTreeLikelihood::getData() const { return _data; }


/******************************************************************************/

ParameterList AbstractTreeLikelihood::getParameters() const 
throw (Exception)
{
	return _parameters;
}

/******************************************************************************/

double AbstractTreeLikelihood::getParameter(const string & name) const
throw (ParameterNotFoundException)
{
	Parameter * p = _parameters.getParameter(name);
	if(p == NULL) throw ParameterNotFoundException("AbstractTreeLikelihood::getParameter", name);
	return p -> getValue();
}

/******************************************************************************/

void AbstractTreeLikelihood::setAllParametersValues(const ParameterList & params)
throw (ParameterNotFoundException, ConstraintException)
{
	_parameters.setAllParametersValues(params);
}

/******************************************************************************/

void AbstractTreeLikelihood::setParameterValue(const string & name, double value)
throw (ParameterNotFoundException, ConstraintException)
{
	_parameters.setParameterValue(name, value);
}

/******************************************************************************/

void AbstractTreeLikelihood::setParametersValues(const ParameterList & params)
throw (ParameterNotFoundException, ConstraintException)
{
	_parameters.setParametersValues(params);
}

/******************************************************************************/

void AbstractTreeLikelihood::matchParametersValues(const ParameterList & params)
throw (ConstraintException)
{
	_parameters.matchParametersValues(params);
}

/******************************************************************************/

Tree * AbstractTreeLikelihood::getTree() const { return _tree; }

/******************************************************************************/

Vdouble AbstractTreeLikelihood::getLikelihoodForEachSite() const {
	Vdouble l(_data -> getNumberOfSites());
	for(unsigned int i = 0; i < l.size(); i++) l[i] = getLikelihoodForASite(i);
	return l;
}

/******************************************************************************/

Vdouble AbstractTreeLikelihood::getLogLikelihoodForEachSite() const {
	Vdouble l(_data -> getNumberOfSites());
	for(unsigned int i = 0; i < l.size(); i++) l[i] = getLogLikelihoodForASite(i);
	return l;
}

/******************************************************************************/
