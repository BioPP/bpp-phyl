//CodonBiasUse

//
// File: CodonFitnessSubstitutionModel.cpp
// Created by: Fanny Pouyet 
// Created on: February 2012
//
 
/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

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
  with loading,  using,  modifying and/or developi_ng or reproducing the
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


#include "CodonFitnessSubstitutionModel.h"
using namespace bpp;
using namespace std;

/*
  Ok, alors soit le modèle O avec F (fitness identique pour les codons synonymes) ddl=19
  Soit le modèle avec 61 F différents. 
  Puis soit otimisé soit observé. Commencons par coder l'observé.
  Le plus simple,
  1/ Creer une liste des Fitness de chaque codon selon l'un des 2 modèles proposé - 0 ou non. On doit avoir des règles précises ---> FrequenciesSet.h (cf AbstractCodonPhaseFrequenciesSubstitutionModel)

  2/ Récupérer le a et les pi_* qui devraient etres respectivement disponible, dans GTR et dans ??(mutation bias nt ??)
  Pour les a, structure a voir dans GTR mais c'est réversible !!
  Pour les pi_*, est-ce qu'ils changent en fonction de leur position dans le codon - 1 à 3 ? la somme fait 1.
  ---> voir HKY / YN98 / K80
  3/ La formule c'est qij= a*(Fj-Fi)*pi_k* exp(Fj)/(exp(Fj)-exp(Fi))
  En moyenne, on doit avoir pi_k1*pi_k2*pi_k3*exp(Fj)=1 a une constante près


  ok mais qu'est ce qu'on retourne ??????????


*/
//Dans l'équation 3, les a(ik,jk) correspondent au modèle nucléotidique de base. Donc ce n'est pas à implémenter, mais appelé lorsque le modèle complet est construit. Par exemple, si tu regardes la ligne 52 de YN98.cpp, tu vois un exemple de telle construction, avec le même K80 modèle en chaque position. Et d'ailleurs le constructeur de CodonDistance.... a besoin d'un modèle (ou de 3) nucléotidique pour être construit (cf les constructeurs dans les .h). 
// ils considèrent des fitness relatives entre codons synonymes, donc une est fixée à 0 et les autres sont paramétrisés. Si tu regardes le constructeur de AbstractCodonDistanceSubstitutionModel tu vois comment on ajoute des paramètres à un modèle. En fait tu peux complètement te baser sur les fichiers de cette classe, ce ne sera pas beaucoup plus compliqué. Il faut juste y mettre plus de paramètres, et la fonction getMulRate tient compte de l'usage du code qui est paramétré.

//pour normaliser les F_I est de dire que leur somme vaut 1, comme des fréquences. A ce moment, tu peux dire que le modèle dépend d'un FrequenciesSet , et donc tu n'as pas à gérer les 60 
/******************************************************************************/
CodonFitnessSubstitutionModel::CodonFitnessSubstitutionModel(const GeneticCode* palph,
                                                             NucleotideSubstitutionModel* pmod,
                                                             FrequenciesSet* pfit,
                                                             FrequenciesSet* pfreq,
                                                             const AlphabetIndex2<double>* pdist) :
  AbstractParameterAliasable("CodonDistFitFreq."),
  AbstractSubstitutionModel(dynamic_cast<const CodonAlphabet*>(palph->getSourceAlphabet()), "CodonDistFitFreq."),
  AbstractWordSubstitutionModel(dynamic_cast<const CodonAlphabet*>(palph->getSourceAlphabet()), "CodonDistFitFreq."),
  AbstractCodonSubstitutionModel(dynamic_cast<const CodonAlphabet*>(palph->getSourceAlphabet()), pmod, "CodonDistFitFreq."),
  AbstractCodonDistanceSubstitutionModel(palph, pdist, "CodonDistFitFreq."),
  AbstractCodonPhaseFrequenciesSubstitutionModel(pfreq, "CodonDistFitFreq."),
  AbstractCodonFitnessSubstitutionModel(pfit, "CodonDistFitFreq.")
{
  updateMatrices();
}

CodonFitnessSubstitutionModel::CodonFitnessSubstitutionModel(const GeneticCode* palph,
                                                             NucleotideSubstitutionModel* pmod1,
                                                             NucleotideSubstitutionModel* pmod2,
                                                             NucleotideSubstitutionModel* pmod3,
                                                             FrequenciesSet* pfit,
                                                             FrequenciesSet* pfreq,
                                                             const AlphabetIndex2<double>* pdist) :
  AbstractParameterAliasable("CodonDistFitFreq."),
  AbstractSubstitutionModel(dynamic_cast<const CodonAlphabet*>(palph->getSourceAlphabet()), "CodonDistFitFreq."),
  AbstractWordSubstitutionModel(dynamic_cast<const CodonAlphabet*>(palph->getSourceAlphabet()), "CodonDistFitFreq."),
  AbstractCodonSubstitutionModel(dynamic_cast<const CodonAlphabet*>(palph->getSourceAlphabet()), pmod1, pmod2, pmod3, "CodonDistFitFreq."),
  AbstractCodonDistanceSubstitutionModel(palph, pdist, "CodonDistFitFreq."),
  AbstractCodonPhaseFrequenciesSubstitutionModel(pfreq, "CodonDistFitFreq."),
  AbstractCodonFitnessSubstitutionModel(pfit,"CodonDistFitFreq.")
{
  updateMatrices();
}

string CodonFitnessSubstitutionModel::getName() const
{
  return ("CodonDistFitFreq.");
}

void CodonFitnessSubstitutionModel::fireParameterChanged(const ParameterList& parameters)
{
  AbstractCodonDistanceSubstitutionModel::fireParameterChanged(parameters);
  AbstractCodonPhaseFrequenciesSubstitutionModel::fireParameterChanged(parameters);
  AbstractCodonFitnessSubstitutionModel::fireParameterChanged(parameters);

  AbstractCodonSubstitutionModel::fireParameterChanged(parameters);
}

double CodonFitnessSubstitutionModel::getCodonsMulRate(unsigned int i, unsigned int j) const
{
  return AbstractCodonDistanceSubstitutionModel::getCodonsMulRate(i,j)
    * AbstractCodonSubstitutionModel::getCodonsMulRate(i,j)
    * AbstractCodonPhaseFrequenciesSubstitutionModel::getCodonsMulRate(i,j)
    * AbstractCodonFitnessSubstitutionModel::getCodonsMulRate(i,j);
}

void CodonFitnessSubstitutionModel::setNamespace(const std::string& st)
{
  AbstractParameterAliasable::setNamespace(st);
  AbstractCodonSubstitutionModel::setNamespace(st);
  AbstractCodonDistanceSubstitutionModel::setNamespace(st);
  AbstractCodonPhaseFrequenciesSubstitutionModel::setNamespace(st); 
  AbstractCodonFitnessSubstitutionModel::setNamespace(st);
}

void CodonFitnessSubstitutionModel::setFreq(map<int,double>& frequencies)
{
  AbstractCodonPhaseFrequenciesSubstitutionModel::setFreq(frequencies);
}


