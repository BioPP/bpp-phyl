//
// File: AbstractTreeLikelihood.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Fri Oct 17 17:57:21 2003
//

#include "AbstractTreeLikelihood.h"

/******************************************************************************/

AbstractTreeLikelihood::~AbstractTreeLikelihood() {}

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
	fireParameterChanged(_parameters);
}

/******************************************************************************/

void AbstractTreeLikelihood::setParameterValue(const string & name, double value)
throw (ParameterNotFoundException, ConstraintException)
{
	_parameters.setParameterValue(name, value);
	fireParameterChanged(_parameters.subList(name));
}

/******************************************************************************/

void AbstractTreeLikelihood::setParametersValues(const ParameterList & params)
throw (ParameterNotFoundException, ConstraintException)
{
	_parameters.setParametersValues(params);
	fireParameterChanged(params);
}

/******************************************************************************/

void AbstractTreeLikelihood::matchParametersValues(const ParameterList & params)
throw (ConstraintException)
{
	_parameters.matchParametersValues(params);
	fireParameterChanged(_parameters); //Conservative, should be the sublist...
}

/******************************************************************************/

Tree * AbstractTreeLikelihood::getTree() const { return _tree; }

/******************************************************************************/
Vdouble AbstractTreeLikelihood::getLikelihoodForEachSite() const {
	Vdouble l(getNumberOfSites());
	for(unsigned int i = 0; i < l.size(); i++) l[i] = getLikelihoodForASite(i);
	return l;
}

/******************************************************************************/

Vdouble AbstractTreeLikelihood::getLogLikelihoodForEachSite() const {
	Vdouble l(getNumberOfSites());
	for(unsigned int i = 0; i < l.size(); i++) l[i] = getLogLikelihoodForASite(i);
	return l;
}

/******************************************************************************/

VVdouble AbstractTreeLikelihood::getLikelihoodForEachSiteForEachState() const {
	VVdouble l(getNumberOfSites());
	for(unsigned int i = 0; i < l.size(); i++) {
		Vdouble * l_i = & l[i];
		l_i -> resize(getNumberOfStates());
		for(unsigned int x = 0; x < l_i -> size(); x++) {
			(* l_i)[x] = getLikelihoodForASiteForAState(i, x);
		}
	}
	return l;
}

/******************************************************************************/

VVdouble AbstractTreeLikelihood::getLogLikelihoodForEachSiteForEachState() const {
	VVdouble l(getNumberOfSites());
	for(unsigned int i = 0; i < l.size(); i++) {
		Vdouble * l_i = & l[i];
		l_i -> resize(getNumberOfStates());
		for(unsigned int x = 0; x < l_i -> size(); x++) {
			(* l_i)[x] = getLogLikelihoodForASiteForAState(i, x);
		}
	}
	return l;
}

/******************************************************************************/

