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
                       AbstractSubstitutionModel * _model,  
                       map<string,AbstractDiscreteDistribution*>  _parametersdistributionslist)
  
  :AbstractSubstitutionModel(alpha,"MixedModel.")
{ 

  unsigned int c, i;
  double d;
  string s1, s2, t;
  map<string, AbstractDiscreteDistribution*>::iterator it;
  
  // Initialization of _distributionmap.

  vector<string> parnames=_model->getParameters().getParameterNames();
  
  for(i=0; i<_model->getNumberOfParameters(); i++){

    s1=parnames[i];
    s2=_model->getParameterNameWithoutNamespace(s1);
    
    if (_parametersdistributionslist.find(s2)!=_parametersdistributionslist.end())
      _distributionmap["MixedModel."+s1]=dynamic_cast<AbstractDiscreteDistribution*>(_parametersdistributionslist.find(s2)->second->clone());
    else
      _distributionmap["MixedModel."+s1]= new ConstantDistribution(_model->getParameterValue(s2));
    
    
    if (dynamic_cast<ConstantDistribution*>(_distributionmap["MixedModel."+s1])==NULL)
      _distributionmap["MixedModel."+s1]->setNamespace("MixedModel."+s1+"."+_distributionmap["MixedModel."+s1]->getNamespace());
    else
      _distributionmap["MixedModel."+s1]->setNamespace("MixedModel."+s1+".");
  }

  
  // Initialization of _modelscontainer.

  c=1;

  for(it=_distributionmap.begin(); it!= _distributionmap.end(); it++)
    c*=it->second->getNumberOfCategories();

  for(i=0; i<c ; i++){
    _modelscontainer.push_back(_model->clone());
    _modelscontainer[i]->setNamespace("MixedModel."+_model->getNamespace());
  }

  // Initialization of _parameters.

  Constraint* pc;
  
  for(it=_distributionmap.begin(); it!=_distributionmap.end(); it++){

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


AbstractSubstitutionModel* MixedModel::getNModel(int i)
{
  return _modelscontainer[i];
}

int MixedModel::getNumberOfModels() const
{
  return _modelscontainer.size();
}

MixedModel::~MixedModel()
{
  unsigned int i;
  map<string, AbstractDiscreteDistribution*>::iterator it;
  
  for (it= _distributionmap.begin(); it!= _distributionmap.end(); it++)
    delete it->second;
  
  for(i = 0; i < _modelscontainer.size(); i++){
    delete _modelscontainer[i];
  }
  
}

void MixedModel::updateMatrices()
{
  string s,t;
  unsigned int i, j, l;
  double d;
  ParameterList pl;
  map<string, AbstractDiscreteDistribution*>::iterator it;

  // Update of distribution parameters from the parameters_ member
  // data. (reverse operation compared to what has been done in the
  // constructor).

  for( it=_distributionmap.begin(); it!=_distributionmap.end(); it++){

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

  
  for(i=0; i<_modelscontainer.size(); i++){
    
    j=i;
    for( it=_distributionmap.begin(); it!=_distributionmap.end(); it++){
      s=it->first;
      l=j%it->second->getNumberOfCategories();
      
      d=_distributionmap.find(s)->second->getCategory(l);
      
      if (pl.hasParameter(s))
	pl.setParameterValue(s,d);
      else
	pl.addParameter(Parameter(s,d));

      j=j/it->second->getNumberOfCategories();
      
    }

    _modelscontainer[i]->matchParametersValues(pl);
  }

  
  for (i=0;i<getNumberOfStates();i++){
    freq_[i]=0;
    for( j=0; j<_modelscontainer.size(); j++)
      freq_[i]+=_modelscontainer[j]->freq(i);
    freq_[i]/=_modelscontainer.size();
  }

}


void MixedModel::setFreq(std::map<int,double>& m){
  throw Exception("setFreq method is not available for MixedModel.");
}

unsigned int MixedModel::getNumberOfStates() const
{
  return _modelscontainer[0]->getNumberOfStates();
}

