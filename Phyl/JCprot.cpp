//
// File: JCprot.cpp
// Created by:  <Julien.Dutheil@univ-montp2.fr>
// Created on: Tue May 27 16:04:36 2003
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

#include "JCprot.h"

#include <cmath>

/******************************************************************************/

JCprot::JCprot(const ProteicAlphabet * alpha): ProteinSubstitutionModel(alpha), AbstractSubstitutionModel(alpha)
{
	_parameters = ParameterList(); //no parameters for this model.	
	
	// Frequencies:
	for(unsigned int i = 0; i < 20; i++) _freq[i] = 1. / 20.;	
	updateMatrices();
}

JCprot::~JCprot() {}

/******************************************************************************/
	
void JCprot::updateMatrices()
{
	// Generator:
	for(unsigned int i = 0; i < 20; i++) {
		for(unsigned int j = 0; j < 20; j++) {
			_generator(i, j) = (i == j) ? -1. : 1./19.;
		}
	}
	
	// Eigen values:
	_eigenValues[0] = 0;
	for(unsigned int i = 1; i < 20; i++) _eigenValues[i] = -20. / 19.;
	
	// Eigen vectors:
	// todo!
}
	
/******************************************************************************/

double JCprot::Pij_t(int i, int j, double d) const {
	if(i == j) return 1./20. + 19./20. * exp(- 20./19. * d);
	else       return 1./20. -  1./20. * exp(- 20./19. * d);
}

/******************************************************************************/

double JCprot::dPij_dt(int i, int j, double d) const {
	if(i == j) return -        exp(- 20./19. * d);
	else       return 1./19. * exp(- 20./19. * d);
}

/******************************************************************************/

double JCprot::d2Pij_dt2(int i, int j, double d) const {
	if(i == j) return   20./19.  * exp(- 20./19. * d);
	else       return - 20./361. * exp(- 20./19. * d);
}

/******************************************************************************/

Mat JCprot::getPij_t(double d) const {
	Mat p(_size, _size);
	for(unsigned int i = 0; i < _size; i++) {
		for(unsigned int j = 0; j < _size; j++) {
			p(i,j) = (i==j) ? 1./20. + 19./20. * exp(- 20./19. * d) : 1./20. - 1./20. * exp(- 20./19. * d);
		}
	}
	return p;
}

Mat JCprot::getdPij_dt(double d) const {
	Mat p(_size, _size);
	for(unsigned int i = 0; i < _size; i++) {
		for(unsigned int j = 0; j < _size; j++) {
			p(i,j) = (i==j) ? - exp(- 20./19. * d) : 1./19. * exp(- 20./19. * d);
		}
	}
	return p;
}

Mat JCprot::getd2Pij_dt2(double d) const {
	Mat p(_size, _size);
	for(unsigned int i = 0; i < _size; i++) {
		for(unsigned int j = 0; j < _size; j++) {
			p(i,j) = (i==j) ? 20./19. * exp(- 20./19. * d) : - 20./361. * exp(- 420./19. * d);
		}
	}
	return p;
}

/******************************************************************************/

string JCprot::getName() const {
	return string("Jukes and Cantor (1969) for proteins");
}

/******************************************************************************/
