//
// File: PseudoNewtonOptimizer.cpp
// Created by: Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Tue Nov 16 12:33 2004
//

/*
Copyright ou � ou Copr. CNRS, (16 Novembre 2004) 

Julien.Dutheil@univ-montp2.fr

Ce logiciel est un programme informatique servant � fournir des classes
pour l'analyse de donn�es phylog�n�tiques.

Ce logiciel est r�gi par la licence CeCILL soumise au droit fran�ais et
respectant les principes de diffusion des logiciels libres. Vous pouvez
utiliser, modifier et/ou redistribuer ce programme sous les conditions
de la licence CeCILL telle que diffus�e par le CEA, le CNRS et l'INRIA 
sur le site "http://www.cecill.info".

En contrepartie de l'accessibilit� au code source et des droits de copie,
de modification et de redistribution accord�s par cette licence, il n'est
offert aux utilisateurs qu'une garantie limit�e.  Pour les m�mes raisons,
seule une responsabilit� restreinte p�se sur l'auteur du programme,  le
titulaire des droits patrimoniaux et les conc�dants successifs.

A cet �gard  l'attention de l'utilisateur est attir�e sur les risques
associ�s au chargement,  � l'utilisation,  � la modification et/ou au
d�veloppement et � la reproduction du logiciel par l'utilisateur �tant 
donn� sa sp�cificit� de logiciel libre, qui peut le rendre complexe � 
manipuler et qui le r�serve donc � des d�veloppeurs et des professionnels
avertis poss�dant  des  connaissances  informatiques approfondies.  Les
utilisateurs sont donc invit�s � charger  et  tester  l'ad�quation  du
logiciel � leurs besoins dans des conditions permettant d'assurer la
s�curit� de leurs syst�mes et ou de leurs donn�es et, plus g�n�ralement, 
� l'utiliser et l'exploiter dans les m�mes conditions de s�curit�. 

Le fait que vous puissiez acc�der � cet en-t�te signifie que vous avez 
pris connaissance de la licence CeCILL, et que vous en avez accept� les
termes.
*/

/*
Copyright or � or Copr. CNRS, (November 16, 2004)

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

#include "PseudoNewtonOptimizer.h"

// From NumCalc:
#include <NumCalc/VectorTools.h>

// From Utils:
#include <Utils/TextTools.h>

/**************************************************************************/
     
PseudoNewtonOptimizer::PNStopCondition::PNStopCondition(PseudoNewtonOptimizer * pno):
	AbstractOptimizationStopCondition(pno) {}

PseudoNewtonOptimizer::PNStopCondition::~PNStopCondition() {}
          
bool PseudoNewtonOptimizer::PNStopCondition::isToleranceReached() const
{
	return abs<double>(
			dynamic_cast<const PseudoNewtonOptimizer *>(_optimizer) -> _currentValue -
			dynamic_cast<const PseudoNewtonOptimizer *>(_optimizer) -> _previousValue)
		< _tolerance; 
}
   
/**************************************************************************/
	
PseudoNewtonOptimizer::PseudoNewtonOptimizer(DerivableSecondOrder * function) :
  AbstractOptimizer(function)
{
	_defaultStopCondition = new PNStopCondition(this);
	_stopCondition = _defaultStopCondition;
}

/**************************************************************************/

PseudoNewtonOptimizer::~PseudoNewtonOptimizer()
{
	delete _defaultStopCondition;
}

/**************************************************************************/

void PseudoNewtonOptimizer::init(const ParameterList & params) throw (Exception)
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

	// Initialize stop condition:
  //_stopCondition -> isToleranceReached();
}

/**************************************************************************/

inline double PseudoNewtonOptimizer::step() throw (Exception)
{
	if(_verbose > 0) { cout << "*"; cout.flush(); }
	// Compute derivative at current point:
	Vdouble movements(_n);
	ParameterList newPoint = _parameters;
	double newValue;
	_function -> setParameters(_parameters);
	for(unsigned int i = 0; i < _n; i++) {
	  double  firstOrderDerivative = dynamic_cast<const DerivableSecondOrder *>(_function) -> getFirstOrderDerivative(_params[i]);
		double secondOrderDerivative = dynamic_cast<const DerivableSecondOrder *>(_function) -> getSecondOrderDerivative(_params[i]);
		if(secondOrderDerivative < 0) {
			printMessage("!!! Second order derivative is negative for parameter " + _params[i] + ". No move performed.");
			movements[i] = 0;  // We want to reach a minimum, not a maximum!
		}
		else movements[i] = firstOrderDerivative / secondOrderDerivative;
		//DEBUG: cout << "PN[" << i << "]=" << movements[i] << "\t " << firstOrderDerivative << "\t" << secondOrderDerivative << endl;
		newPoint.setParameterValue(_params[i], _parameters.getParameter(_params[i]) -> getValue() - movements[i]);
	}
	newValue = _function -> f(newPoint);

	// Check newValue:
	unsigned int count = 0;
	while(newValue > _currentValue) {
		printMessage("!!! Function at new point is greater than at current point: " + TextTools::toString(newValue) + ">" + TextTools::toString(_currentValue) + ". Applying Felsenstein-Churchill correction.");
		for(unsigned int i = 0; i < _n; i++) {
			movements[i] = movements[i] / 2;
			newPoint.setParameterValue(_params[i], _parameters.getParameter(_params[i]) -> getValue() - movements[i]);
		}
		count++;
		if(count > 10000) throw Exception("PseudoNewtonOptimizer::step(). Felsenstein-Churchill correction applied more than 10000 times.");
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

double PseudoNewtonOptimizer::optimize() throw (Exception)
{
	_tolIsReached = false;
	for (_nbEval = 0; _nbEval < _nbEvalMax && ! _tolIsReached; _nbEval++) {
		step();
	}
	return _currentValue;
}

/**************************************************************************/

double PseudoNewtonOptimizer::getFunctionValue() const throw (Exception) { return _currentValue; }

/**************************************************************************/

