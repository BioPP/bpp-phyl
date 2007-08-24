//
// File: SubstitutionModelFactory.cpp
// Created by: Julien Dutheil
//             Vincent Ranwez
// Created on: Fri apr 14 11:11 2006
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004, 2005, 2006)

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

#include "SubstitutionModelFactory.h"

// From SeqLib:
#include <Seq/AlphabetTools.h>

// From the STL:
#include <algorithm>

const string SubstitutionModelFactory::JUKES_CANTOR            = "JC69";
const string SubstitutionModelFactory::KIMURA_2P               = "K80";
const string SubstitutionModelFactory::HASEGAWA_KISHINO_YANO   = "HKY85";
const string SubstitutionModelFactory::TAMURA_NEI              = "TN93";
const string SubstitutionModelFactory::GENERAL_TIME_REVERSIBLE = "HKY85";
const string SubstitutionModelFactory::TAMURA                  = "T92";
const string SubstitutionModelFactory::FELSENSTEIN84           = "F84";
const string SubstitutionModelFactory::JOHN_TAYLOR_THORNTON    = "JTT92";
const string SubstitutionModelFactory::DAYHOFF_SCHWARTZ_ORCUTT = "DSO78";

SubstitutionModel * SubstitutionModelFactory::createModel(const string& modelName) const throw (AlphabetException, Exception)
{
  if(modelName == JUKES_CANTOR)
  {
    if(AlphabetTools::isNucleicAlphabet(_alphabet))
      return new JCnuc(dynamic_cast<const NucleicAlphabet *>(_alphabet));
    else
      return new JCprot(dynamic_cast<const ProteicAlphabet *>(_alphabet));
  }
  else if(modelName == KIMURA_2P)
  {
    try { 
      return new K80(dynamic_cast<const NucleicAlphabet *>(_alphabet));
    } catch(Exception & e) {
      throw AlphabetException("SubstitutionModelFactory::createInstanceOf(). K80 model requires a nucleotide alphabet.", _alphabet);
    }
  }
  else if(modelName == TAMURA)
  {
    try { 
      return new T92(dynamic_cast<const NucleicAlphabet *>(_alphabet));
    } catch(Exception & e) {
      throw AlphabetException("SubstitutionModelFactory::createInstanceOf(). T92 model requires a nucleotide alphabet.", _alphabet);
    }
  }
  else if(modelName == FELSENSTEIN84)
  {
    try { 
      return new F84(dynamic_cast<const NucleicAlphabet *>(_alphabet));
    } catch(Exception & e) {
      throw AlphabetException("SubstitutionModelFactory::createInstanceOf(). T92 model requires a nucleotide alphabet.", _alphabet);
    }
  }
  else if(modelName == TAMURA_NEI)
  {
    try { 
      return new TN93(dynamic_cast<const NucleicAlphabet *>(_alphabet));
    } catch(Exception & e) {
      throw AlphabetException("SubstitutionModelFactory::createInstanceOf(). TN93 model requires a nucleotide alphabet.", _alphabet);
    }
  }
  else if(modelName == HASEGAWA_KISHINO_YANO)
  {
    try { 
      return new HKY85(dynamic_cast<const NucleicAlphabet *>(_alphabet));
    } catch(Exception & e) {
      throw AlphabetException("SubstitutionModelFactory::createInstanceOf(). HKY85 model requires a nucleotide alphabet.", _alphabet);
    }
  }
  else if(modelName == GENERAL_TIME_REVERSIBLE)
  {
    try { 
      return new GTR(dynamic_cast<const NucleicAlphabet *>(_alphabet));
    } catch(Exception & e) {
      throw AlphabetException("SubstitutionModelFactory::createInstanceOf(). GTR model requires a nucleotide alphabet.", _alphabet);
    }
  }
  else if(modelName == JOHN_TAYLOR_THORNTON)
  {
    try { 
      return new JTT92(dynamic_cast<const ProteicAlphabet *>(_alphabet));
    } catch(Exception & e) {
      throw AlphabetException("SubstitutionModelFactory::createInstanceOf(). JTT92 model requires a protein alphabet.", _alphabet);
    }
  }
  else if(modelName == DAYHOFF_SCHWARTZ_ORCUTT)
  {
    try { 
      return new DSO78(dynamic_cast<const ProteicAlphabet *>(_alphabet));
    } catch(Exception & e) {
      throw AlphabetException("SubstitutionModelFactory::createInstanceOf(). DSO78 model requires a protein alphabet.", _alphabet);
    }
  }
  else throw Exception("SubstitutionModelFactory::createModel(). Unknown model: " + modelName);
}

SubstitutionModelSet * SubstitutionModelFactory::createHomogeneousModelSet(const string& modelName, const Tree * tree) const throw (AlphabetException, Exception)
{
  SubstitutionModel * model = createModel(modelName);
  SubstitutionModelSet * modelSet = new SubstitutionModelSet(_alphabet);
  //We assign this model to all nodes in the tree (excepted root node), and link all parameters with it.
  vector<int> ids = tree->getNodesId();
  int rootId = tree->getRootId();
  remove(ids.begin(), ids.end(), rootId);
  modelSet->addModel(model, ids, model->getParameters().getParameterNames());
  return modelSet;
}

SubstitutionModelSet * SubstitutionModelFactory::createNonHomogeneousModelSet(const string& modelName, const Tree * tree, const vector<string> & globalParameterNames) const throw (AlphabetException, Exception)
{
  SubstitutionModel * model = createModel(modelName); //template model.
  ParameterList globalParameters, branchParameters;
  globalParameters = model->getParameters();
  for(unsigned int i = globalParameters.size(); i > 0; i--)
  {
    if(find(globalParameterNames.begin(), globalParameterNames.end(), globalParameters[i]->getName()) == globalParameterNames.end())
    {
      //not a global parameter:
      branchParameters.addParameter(*globalParameters[i]);
      globalParameters.erase(globalParameters.begin() + i);
    }
  }
  SubstitutionModelSet * modelSet = new SubstitutionModelSet(_alphabet);
  //We assign a copy of this model to all nodes in the tree (excepted root node), and link all parameters with it.
  vector<int> ids = tree->getNodesId();
  int rootId = tree->getRootId();
  remove(ids.begin(), ids.end(), rootId);
  for(unsigned int i = 0; i < ids.size(); i++)
  {
    modelSet->addModel(dynamic_cast<SubstitutionModel *>(model->clone()), vector<int>(1, ids[i]), branchParameters.getParameterNames());
  }
  //Now add global parameters to all nodes:
  modelSet->addParameters(globalParameters, ids);
  delete model; //delete template model.
  return modelSet;
}

