//
// File: MixedModel.cpp
// Created by: David Fournier, Laurent Gueguen
//

/*
  Copyright or © or Copr. CNRS, (November 16, 2004)

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

#include "MixedModel.h"
#include <string>

using namespace bpp;
using namespace std;


MixedModel::MixedModel(const Alphabet * alpha,
                       SubstitutionModel * _model,  
                       map<string,DiscreteDistribution*>  _parametersdistributionslist)
  
  :AbstractSubstitutionModel(alpha,"MixedModel.")
{ 

  unsigned int c, i;
  double d;
  string s1, s2, t;
  map<string, DiscreteDistribution*>::iterator it;
  
  // Initialization of distributionmap_.

  vector<string> parnames=_model->getParameters().getParameterNames();
  
  for(i=0; i<_model->getNumberOfParameters(); i++){

    s1=parnames[i];
    s2=_model->getParameterNameWithoutNamespace(s1);
    
    if (_parametersdistributionslist.find(s2)!=_parametersdistributionslist.end())
      distributionmap_["MixedModel."+s1]=dynamic_cast<DiscreteDistribution*>(_parametersdistributionslist.find(s2)->second->clone());
    else
      distributionmap_["MixedModel."+s1]= new ConstantDistribution(_model->getParameterValue(s2));
    
    
    if (dynamic_cast<ConstantDistribution*>(distributionmap_["MixedModel."+s1])==NULL)
      distributionmap_["MixedModel."+s1]->setNamespace("MixedModel."+s1+"."+distributionmap_["MixedModel."+s1]->getNamespace());
    else
      distributionmap_["MixedModel."+s1]->setNamespace("MixedModel."+s1+".");
  }

  
  // Initialization of modelscontainer_.

  c=1;

  for(it=distributionmap_.begin(); it!= distributionmap_.end(); it++)
    c*=it->second->getNumberOfCategories();

  for(i=0; i<c ; i++){
    modelscontainer_.push_back(_model->clone());
    modelscontainer_[i]->setNamespace("MixedModel."+_model->getNamespace());
  }

  // Initialization of _parameters.

  Constraint* pc;
  
  for(it=distributionmap_.begin(); it!=distributionmap_.end(); it++){

    if (dynamic_cast<ConstantDistribution*>(it->second)==NULL){
      for(i=0; i!=it->second->getNumberOfParameters(); i++){ 
        
	t=it->second->getParameters().getParameterNames()[i];
      	d=it->second->getParameters().getParameter(t).getValue();
        
        addParameter_(Parameter(t,d));
      }
    }
    else{
      t=it->first;
      pc=_model->getParameter(_model->getParameterNameWithoutNamespace(getParameterNameWithoutNamespace(t))).getConstraint()->clone();
      d=it->second->getCategory(0);
      addParameter_(Parameter(t,d,pc,true));

    }

  }

  updateMatrices();

}


SubstitutionModel* MixedModel::getNModel(int i)
{
  return modelscontainer_[i];
}

int MixedModel::getNumberOfModels() const
{
  return modelscontainer_.size();
}

MixedModel::~MixedModel()
{
  unsigned int i;
  map<string, DiscreteDistribution*>::iterator it;
  
  for (it= distributionmap_.begin(); it!= distributionmap_.end(); it++)
    delete it->second;
  
  for(i = 0; i < modelscontainer_.size(); i++){
    delete modelscontainer_[i];
  }
  
}

void MixedModel::updateMatrices()
{
  string s,t;
  unsigned int i, j, l;
  double d;
  ParameterList pl;
  map<string, DiscreteDistribution*>::iterator it;

  // Update of distribution parameters from the parameters_ member
  // data. (reverse operation compared to what has been done in the
  // constructor).

  for( it=distributionmap_.begin(); it!=distributionmap_.end(); it++){

    if (dynamic_cast<ConstantDistribution*>(it->second)==NULL){
      for(i=0; i<it->second->getNumberOfParameters(); i++){ 
        
        t=it->second->getParameters().getParameterNames()[i];
        d=getParameter(getParameterNameWithoutNamespace(t)).getValue();
        it->second->setParameterValue(it->second->getParameterNameWithoutNamespace(t),d); 
      }    
    }
    else{
      t=it->second->getNamespace();
      d=getParameter(getParameterNameWithoutNamespace(t.substr(0,t.length()-1))).getValue();
      it->second->setParameterValue("value",d); 

    }
  }

  
  for(i=0; i<modelscontainer_.size(); i++){
    
    j=i;
    for( it=distributionmap_.begin(); it!=distributionmap_.end(); it++){
      s=it->first;
      l=j%it->second->getNumberOfCategories();
      
      d=distributionmap_.find(s)->second->getCategory(l);
      
      if (pl.hasParameter(s))
	pl.setParameterValue(s,d);
      else
	pl.addParameter(Parameter(s,d));

      j=j/it->second->getNumberOfCategories();
      
    }

    modelscontainer_[i]->matchParametersValues(pl);
  }

  
  for (i=0;i<getNumberOfStates();i++){
    freq_[i]=0;
    for( j=0; j<modelscontainer_.size(); j++)
      freq_[i]+=modelscontainer_[j]->freq(i);
    freq_[i]/=modelscontainer_.size();
  }

}


void MixedModel::setFreq(std::map<int,double>& m){
  throw Exception("setFreq method is not available for MixedModel.");
}

unsigned int MixedModel::getNumberOfStates() const
{
  return modelscontainer_[0]->getNumberOfStates();
}

