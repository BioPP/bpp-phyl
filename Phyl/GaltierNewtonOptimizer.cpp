//
// File: GaltierNewtonOptimizer.cpp
// Created by; jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: ue Nov 16 12:33 2004
//

/*
Copyright ou © ou Copr. Julien Dutheil, (16 Novembre 2004) 

Julien.Dutheil@univ-montp2.fr

Ce logiciel est un programme informatique servant à fournir des classes
pour l'analyse de données phylogénétiques.

Ce logiciel est régi par la licence CeCILL soumise au droit français et
respectant les principes de diffusion des logiciels libres. Vous pouvez
utiliser, modifier et/ou redistribuer ce programme sous les conditions
de la licence CeCILL telle que diffusée par le CEA, le CNRS et l'INRIA 
sur le site "http://www.cecill.info".

En contrepartie de l'accessibilité au code source et des droits de copie,
de modification et de redistribution accordés par cette licence, il n'est
offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
seule une responsabilité restreinte pèse sur l'auteur du programme,  le
titulaire des droits patrimoniaux et les concédants successifs.

A cet égard  l'attention de l'utilisateur est attirée sur les risques
associés au chargement,  à l'utilisation,  à la modification et/ou au
développement et à la reproduction du logiciel par l'utilisateur étant 
donné sa spécificité de logiciel libre, qui peut le rendre complexe à 
manipuler et qui le réserve donc à des développeurs et des professionnels
avertis possédant  des  connaissances  informatiques approfondies.  Les
utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
logiciel à leurs besoins dans des conditions permettant d'assurer la
sécurité de leurs systèmes et ou de leurs données et, plus généralement, 
à l'utiliser et l'exploiter dans les mêmes conditions de sécurité. 

Le fait que vous puissiez accéder à cet en-tête signifie que vous avez 
pris connaissance de la licence CeCILL, et que vous en avez accepté les
termes.
*/

/*
Copyright or © or Copr. Julien Dutheil, (November 16, 2004)

Julien.Dutheil@univ-montp2.fr

This software is a computer program whose purpose is to provide classes
for phylogenetic data analysis.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

/**************************************************************************/

#include "GaltierNewtonOptimizer.h"

// From NumCalc:
#include <NumCalc/VectorTools.h>

/**************************************************************************/
     
GaltierNewtonOptimizer::GNStopCondition::GNStopCondition(GaltierNewtonOptimizer * gno):
	AbstractOptimizationStopCondition(gno) {}

GaltierNewtonOptimizer::GNStopCondition::~GNStopCondition() {}
          
bool GaltierNewtonOptimizer::GNStopCondition::isToleranceReached() const {
	return abs<double>(
			dynamic_cast<const GaltierNewtonOptimizer *>(_optimizer) -> _currentValue -
			dynamic_cast<const GaltierNewtonOptimizer *>(_optimizer) -> _previousValue)
		< _tolerance; 
}
   
/**************************************************************************/
	
GaltierNewtonOptimizer::GaltierNewtonOptimizer(const DerivableSecondOrder * function) :
  AbstractOptimizer(function)
{
	_defaultStopCondition = new GNStopCondition(this);
	_stopCondition = _defaultStopCondition;
}

/**************************************************************************/

GaltierNewtonOptimizer::~GaltierNewtonOptimizer() {
	delete _defaultStopCondition;
}

/**************************************************************************/

void GaltierNewtonOptimizer::init(const ParameterList & params) throw (Exception)
{
	AbstractOptimizer::init(params);
	_function -> setParameters(_parameters);
	_currentValue = _function -> getValue();
	_n = _parameters.size();
	_params = _parameters.getParameterNames();
	for (int j = 0; j < _n; j++) {
		profile(_parameters[j] -> getName() + "\t"); 
	}
	profileln("Function");
	printPoint(_parameters, _currentValue);
}

/**************************************************************************/

inline double GaltierNewtonOptimizer::step() throw (Exception)
{
	cout << "*"; cout.flush();
	// Compute derivative at current point:
	Vdouble movements(_n);
	ParameterList newPoint = _parameters;
	double newValue;
	for(unsigned int i = 0; i < _n; i++) {
	  double  firstOrderDerivative = dynamic_cast<const DerivableSecondOrder *>(_function) -> getFirstOrderDerivative(_params[i]);
		double secondOrderDerivative = dynamic_cast<const DerivableSecondOrder *>(_function) -> getSecondOrderDerivative(_params[i]);
		if(secondOrderDerivative < 0) movements[i] = 0;  // We want to reach a minimum, not a maximum!
		else movements[i] = firstOrderDerivative / secondOrderDerivative;
		newPoint.setParameterValue(_params[i], _parameters.getParameter(_params[i]) -> getValue() - movements[i]);
	}
	newValue = _function -> f(newPoint);

	// Check newValue:
	unsigned int count = 0;
	while(newValue > _currentValue) {
		printMessage("!!! Function at new point is greater than at current point. Applying Felsenstein-Churchill correction.");
		for(unsigned int i = 0; i < _n; i++) {
			movements[i] = movements[i] / 2;
			newPoint.setParameterValue(_params[i], _parameters.getParameter(_params[i]) -> getValue() - movements[i]);
			count++;
			if(count > 1000) throw Exception("GaltierNewtonOptimizer::step(). Felsenstein-Churchill correction applied more than 1000 times.");
		}
		newValue = _function -> f(newPoint);
	}

	_previousPoint = _parameters;
	_previousValue = _currentValue;
	_parameters   = newPoint; // Function as been set to newPoint by the call of f(newPoint).
	_currentValue = newValue;

	_tolIsReached = _stopCondition -> isToleranceReached();

	printPoint(_parameters, _currentValue);
	return _currentValue;
}

/**************************************************************************/

double GaltierNewtonOptimizer::optimize() throw (Exception)
{
	_tolIsReached = false;
	for (_nbEval = 0; _nbEval < _nbEvalMax && ! _tolIsReached; _nbEval++) {
		step();
	}
}

/**************************************************************************/

double GaltierNewtonOptimizer::getFunctionValue() const { return _currentValue; }

/**************************************************************************/

