//
// File: MappingTools.cpp
// Created by: Julien Dutheil, Laurent Guéguen
// Created on: lundi 15 avril 2013, à 15h 49
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 16, 2004, 2005, 2006)

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


#include "MappingTools.h"

#include "UniformizationSubstitutionCount.h"
#include "SubstitutionMappingTools.h"
#include "DecompositionReward.h"
#include "ProbabilisticRewardMapping.h"
#include "RewardMappingTools.h"

#include "../Likelihood/DRTreeLikelihoodTools.h"

#include <Bpp/Seq/AlphabetIndex/UserAlphabetIndex1.h>

using namespace std;
using namespace bpp;

/**************************************************************************************************/

vector< vector<double> > MappingTools::getCountsPerBranch(
                                                          DRTreeLikelihood& drtl,
                                                          const vector<int>& ids,
                                                          SubstitutionModel* model,
                                                          const SubstitutionRegister& reg,
                                                          double threshold)
{
  auto_ptr<SubstitutionCount> count(new UniformizationSubstitutionCount(model, reg.clone()));

  auto_ptr<ProbabilisticSubstitutionMapping> mapping(SubstitutionMappingTools::computeSubstitutionVectors(drtl, *count, false));

  vector< vector<double> > counts(ids.size());
  size_t nbSites = mapping->getNumberOfSites();
  size_t nbTypes = mapping->getNumberOfSubstitutionTypes();
  for (size_t k = 0; k < ids.size(); ++k) {
    vector<double> countsf(nbTypes, 0);
    vector<double> tmp(nbTypes, 0);
    size_t nbIgnored = 0;
    bool error = false;
    for (size_t i = 0; !error && i < nbSites; ++i) {
      double s = 0;
      for (size_t t = 0; t < nbTypes; ++t) {
        tmp[t] = (*mapping)(k, i, t);
        error = isnan(tmp[t]);
        if (error)
          goto ERROR;
        s += tmp[t];
      }
      if (threshold >= 0) {
        if (s <= threshold)
          countsf += tmp;
        else {
          nbIgnored++;
        }
      } else {
        countsf += tmp;
      }
    }

  ERROR:
    if (error) {
      //We do nothing. This happens for small branches.
      ApplicationTools::displayWarning("On branch " + TextTools::toString(ids[k]) + ", counts could not be computed.");
      for (size_t t = 0; t < nbTypes; ++t)
        countsf[t] = 0;
    } else {
      if (nbIgnored > 0) {
        ApplicationTools::displayWarning("On branch " + TextTools::toString(ids[k]) + ", " + TextTools::toString(nbIgnored) + " sites (" + TextTools::toString(ceil(static_cast<double>(nbIgnored * 100) / static_cast<double>(nbSites))) + "%) have been ignored because they are presumably saturated.");
      }
    }
    
    counts[k].resize(countsf.size());
    for (size_t j = 0; j < countsf.size(); ++j) {
      counts[k][j] = countsf[j]; 
    }
  }
  return counts;
}

/**************************************************************************************************/

vector< vector<double> > MappingTools::getNormalizationsPerBranch(
                                                                  DRTreeLikelihood& drtl,
                                                                  const vector<int>& ids,
                                                                  const SubstitutionModel* nullModel,
                                                                  const SubstitutionRegister& reg)
{
  size_t nbTypes = reg.getNumberOfSubstitutionTypes();
  size_t nbStates = nullModel->getAlphabet()->getSize();
  size_t nbSites = drtl.getNumberOfSites();
  
  // compute the AlphabetIndex for each substitutionType
  vector<UserAlphabetIndex1 > usai(nbTypes, UserAlphabetIndex1(nullModel->getAlphabet()));

  for (size_t i=0; i<nbStates; i++)
    for (size_t j=0; j<nbStates; j++)
      if (i!=j)
        {
          size_t nbt=reg.getType((int)i,(int)j);
          if (nbt!=0)
            usai[nbt-1].setIndex((int)i,usai[nbt-1].getIndex((int)i)+nullModel->Qij((int)i,(int)j)); 
        }
  
  // compute the normalization for each substitutionType
  vector< vector<double> > rewards(ids.size());

  for (size_t k = 0; k < ids.size(); ++k)
    rewards[k].resize(nbTypes);

  for (size_t nbt=0; nbt< nbTypes; nbt++)
    {
      auto_ptr<Reward> reward(new DecompositionReward(nullModel, &usai[nbt]));

      auto_ptr<ProbabilisticRewardMapping> mapping(RewardMappingTools::computeRewardVectors(drtl, *reward, false));
      
      for (size_t k = 0; k < ids.size(); ++k) {
        double s = 0;
        for (size_t i = 0; i < nbSites; ++i) {
          double tmp = (*mapping)(k, i);
          if (isnan(tmp))
            {
              ApplicationTools::displayWarning("On branch " + TextTools::toString(ids[k]) + ", reward for type " + reg.getTypeName(nbt) + " could not be computed.");
              s=0;
              break;
            }
          s += tmp;
        }
        rewards[k][nbt]=s;
      }
      reward.release();
      mapping.release();
    }
  return rewards;
}

/**************************************************************************************************/

vector< vector<double> > MappingTools::getNormalizationsPerBranch(
                                                                  DRTreeLikelihood& drtl,
                                                                  const vector<int>& ids,
                                                                  const SubstitutionModelSet* nullModelSet,
                                                                  const SubstitutionRegister& reg)
{
  size_t nbTypes = reg.getNumberOfSubstitutionTypes();
  size_t nbStates = nullModelSet->getAlphabet()->getSize();
  size_t nbSites = drtl.getNumberOfSites();
  size_t nbModels = nullModelSet->getNumberOfModels();

  // compute the AlphabetIndex for each substitutionType
  vector<vector<UserAlphabetIndex1 > > usai(nbModels, vector<UserAlphabetIndex1>(nbTypes, UserAlphabetIndex1(nullModelSet->getAlphabet())));

  for (size_t nbm=0; nbm<nbModels; nbm++)
    {
      const SubstitutionModel* modn=nullModelSet->getModel(nbm);
      for (size_t i=0; i<nbStates; i++)
        for (size_t j=0; j<nbStates; j++)
          if (i!=j)
            {
              size_t nbt=reg.getType((int)i,(int)j);
              if (nbt!=0)
                usai[nbm][nbt-1].setIndex((int)i,usai[nbm][nbt-1].getIndex((int)i)+modn->Qij((int)i,(int)j)); 
            }
    }

  
  // compute the normalization for each substitutionType
  vector< vector<double> > rewards(ids.size());

  for (size_t k = 0; k < ids.size(); ++k)
    rewards[k].resize(nbTypes);

  for (size_t nbt=0; nbt< nbTypes; nbt++)
    {
      
      for (size_t k = 0; k < ids.size(); ++k) {
        auto_ptr<Reward> reward(new DecompositionReward(nullModelSet->getModelForNode(ids[k]), &usai[nullModelSet->getModelIndexForNode(ids[k])][nbt]));
        auto_ptr<ProbabilisticRewardMapping> mapping(RewardMappingTools::computeRewardVectors(drtl, *reward, false));
        double s = 0;
        for (size_t i = 0; i < nbSites; ++i) {
          double tmp = (*mapping)(k, i);
          if (isnan(tmp))
            {
              ApplicationTools::displayWarning("On branch " + TextTools::toString(ids[k]) + ", reward for type " + reg.getTypeName(nbt) + " could not be computed.");
              s=0;
              break;
            }
          s += tmp;
        }
        rewards[k][nbt]=s;
        reward.release();
        mapping.release();
      }
    }
  return rewards;
}

/**************************************************************************************************/

vector< vector<double> >  MappingTools::getNormalizedCountsPerBranch(
                                                                     DRTreeLikelihood& drtl,
                                                                     const vector<int>& ids,
                                                                     SubstitutionModel* model,
                                                                     SubstitutionModel* nullModel,
                                                                     const SubstitutionRegister& reg,
                                                                     double threshold)
{
  CompleteSubstitutionRegister compreg=CompleteSubstitutionRegister(reg);

  vector< vector<double> > counts=getCountsPerBranch(drtl, ids, model, compreg, threshold);

  vector< vector<double> > factors=getNormalizationsPerBranch(drtl, ids, nullModel, compreg);

  size_t nbTypes = counts[0].size();
    
  for (size_t k = 0; k < ids.size(); ++k) {
    for (size_t t = 0; t < nbTypes; ++t) {
      if (factors[k][t]!=0)
        counts[k][t] /= factors[k][t];
    }
  }
    
  // Sets the sum of complete counts to 1
  for (size_t ibr=0; ibr<counts.size(); ibr++){
    double x=VectorTools::sum(counts[ibr]);
    if (x!=0)
      counts[ibr]/=x;
  }

  // Removes the completion class if needed

  return counts;

}

/**************************************************************************************************/

vector< vector<double> >  MappingTools::getNormalizedCountsPerBranch(
                                                                     DRTreeLikelihood& drtl,
                                                                     const vector<int>& ids,
                                                                     SubstitutionModelSet* modelSet,
                                                                     SubstitutionModelSet* nullModelSet,
                                                                     const SubstitutionRegister& reg,
                                                                     double threshold)
{
  CompleteSubstitutionRegister compreg=CompleteSubstitutionRegister(reg);

  vector< vector<double> > counts=getCountsPerBranch(drtl, ids, modelSet->getModel(0), compreg, threshold);

  size_t nbTypes = counts[0].size();
  vector< vector<double> > factors=getNormalizationsPerBranch(drtl, ids, nullModelSet, compreg);

  for (size_t k = 0; k < ids.size(); ++k) {
    for (size_t t = 0; t < nbTypes; ++t) {
      if (factors[k][t]!=0)
        counts[k][t] /= factors[k][t];
    }
  }
    
  // Sets the sum of complete counts to 1
  for (size_t ibr=0; ibr<counts.size(); ibr++){
    double x=VectorTools::sum(counts[ibr]);
    if (x!=0)
      counts[ibr]/=x;
  }

  // Removes the completion class if needed

  if (reg.getNumberOfSubstitutionTypes()!=nbTypes)
    for (size_t ibr=0; ibr<counts.size(); ibr++)
      counts[ibr].pop_back();

  return counts;
}

/**************************************************************************************************/

vector< vector<double> >  MappingTools::getRelativeCountsPerBranch(
                                                                   DRTreeLikelihood& drtl,
                                                                   const vector<int>& ids,
                                                                   SubstitutionModel* model,
                                                                   const SubstitutionRegister& reg,
                                                                   bool stationarity,
                                                                   double threshold)
{
  vector< vector<double> > counts=getCountsPerBranch(drtl, ids, model, reg, threshold);

  const CategorySubstitutionRegister* creg;  
  if (!stationarity) {
    try {
      creg = &dynamic_cast<const CategorySubstitutionRegister&>(reg);
    } catch (Exception& ex) {
      throw Exception("The stationarity option can only be used with a category substitution register.");
    }
    
    size_t nbTypes = counts[0].size();

    for (size_t k = 0; k < ids.size(); ++k) {
      vector<double> freqs = DRTreeLikelihoodTools::getPosteriorStateFrequencies(drtl, ids[k]);
      //Compute frequencies for types:
      vector<double> freqsTypes(creg->getNumberOfCategories());
      for (size_t i = 0; i < freqs.size(); ++i) {
        size_t c = creg->getCategory(static_cast<int>(i));
        freqsTypes[c - 1] += freqs[i];
      }
      
      //We devide the counts by the frequencies and rescale:
      double s = VectorTools::sum(counts[k]);
      for (size_t t = 0; t < nbTypes; ++t) {
        counts[k][t] /= freqsTypes[creg->getCategoryFrom(t + 1) - 1];
      }
       
      double s2 = VectorTools::sum(counts[k]);
      //Scale:
      counts[k] = (counts[k] / s2) * s;
    }
  }
  
  return counts;
}

/**************************************************************************************************/

void MappingTools::outputTotalCountsPerBranchPerSite(
                                                     string& filename,
                                                     DRTreeLikelihood& drtl,
                                                     const vector<int>& ids,
                                                     SubstitutionModel* model,
                                                     const SubstitutionRegister& reg)
{
  auto_ptr<SubstitutionCount> count(new UniformizationSubstitutionCount(model, reg.clone()));
  auto_ptr<ProbabilisticSubstitutionMapping> smap(SubstitutionMappingTools::computeSubstitutionVectors(drtl, *count, false));

  ofstream file;
  file.open(filename.c_str());

  size_t nbSites = smap->getNumberOfSites();
  size_t nbBr = ids.size();

  vector<int> sdi(nbBr);  //reverse of ids
  for (size_t i = 0; i < nbBr; ++i) {
    for (size_t j = 0; j < nbBr; ++j) {
      if (ids[j] == static_cast<int>(i)) {
        sdi[i] = static_cast<int>(j);
        break;
      }
    }
  }

  file << "sites";
  for (size_t i = 0; i < nbBr; ++i)
    file << "\t" << i ;
  file << endl;
  
  for (size_t k = 0; k < nbSites; ++k) {
    vector<double> countsf = SubstitutionMappingTools::computeTotalSubstitutionVectorForSite(*smap, k);
    file << k;
    for (size_t i = 0; i < nbBr; ++i)
      file << "\t" << countsf[sdi[i]];
    file << endl;
  }
  file.close();
}
