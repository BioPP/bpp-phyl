//
// File: SubstitutionProcessCollection.h
// Created by: Laurent Guéguen
// Created on: mercredi 12 juin 2013, à 14h 07
//

/*
   Copyright or (c) or Copr. Bio++ Development Team, (November 16, 2004)

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

#ifndef _SUBSTITUTIONPROCESSCOLLECTION_H_
#define _SUBSTITUTIONPROCESSCOLLECTION_H_


#include "../Model/SubstitutionModel.h"
#include "../Model/FrequenciesSet/FrequenciesSet.h"

#include "SubstitutionProcessCollectionMember.h"
#include "SubstitutionProcess.h"


#include "ParametrizableTree.h"

#include <Bpp/Numeric/ParametrizableCollection.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/VectorTools.h>

// From Seqlib:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>


namespace bpp
{
/**
 * @brief Collection of Substitution Process, which owns all the
 * necessary objects: Substitution models, frequencies sets, rate
 * distributions and trees.
 *
 * This class contains a set of substitution models sets, linking
 * models, rates and trees.
 *
 * This object has the same parameters names as the owned objects. The
 * updating of the parameters is done through the
 * ParametrizableCollection objects.
 *
 * This class also deals with the parameters associated to the models,
 * distribution, frequencies.
 *
 * It will deal with the parameters associated to the trees when they
 * are parametrized.
 *
 */

class SubstitutionProcessCollection :
  public AbstractParameterAliasable
{
private:
  /**
   * A collection of Substitution Models
   */
  ParametrizableCollection<SubstitutionModel> modelColl_;

  /**
   * A map from each SubstitutionModel number to the SubProcess
   * members that are linked to it.
   *
   */
  
  std::map<size_t, std::vector<size_t> > mModelToSubPro_;
  
  /*
   * A collection of Frequencies Sets
   */
  
  ParametrizableCollection<FrequenciesSet> freqColl_;

  /**
   * A map from each FrequenciesSet number to the SubProcess members
   * that are linked to it.
   *
   */
  
  std::map<size_t, std::vector<size_t> > mFreqToSubPro_;
  
  /*
   * A collection of DiscreteDistributions
   */
  
  ParametrizableCollection<DiscreteDistribution> distColl_;

  /**
   * A map from the DiscreteDistribution numbers to the numbers of
   * categories that are used independently as ConstantDistributions.
   *
   * These ConstantDistributions are stored in distColl_ with numbers
   * 10000*(numberOfDiscreteDistribution+1) + numberOfTheCategory.
   */

  std::map<size_t, std::vector<size_t> > mVConstDist_;
  
   
  /**
   * A map from each DiscreteDistribution number to the SubProcess
   * members that are linked to it.
   *
   */
  
  std::map<size_t, std::vector<size_t> > mDistToSubPro_;
  
  /*
   * A collection of trees
   *
   */

  ParametrizableCollection<ParametrizableTree> treeColl_;

  /**
   * A map from each Tree number to the SubProcess members that are
   * linked to it.
   *
   */
  
  std::map<size_t, std::vector<size_t> > mTreeToSubPro_;
  
  /*
   * A map of SubstitutionProcessCollectionMember
   */

  std::map<size_t, SubstitutionProcessCollectionMember*> mSubProcess_;

public:
  /**
   * @brief Create empty collections.
   *
   */
  
  SubstitutionProcessCollection():
    AbstractParameterAliasable(""),
    modelColl_(),
    mModelToSubPro_(),
    freqColl_(),
    mFreqToSubPro_(),
    distColl_(),
    mVConstDist_(),
    mDistToSubPro_(),
    treeColl_(),
    mTreeToSubPro_(),
    mSubProcess_()
  {
  }


  SubstitutionProcessCollection(const SubstitutionProcessCollection& set);

  SubstitutionProcessCollection& operator=(const SubstitutionProcessCollection& set);

  virtual ~SubstitutionProcessCollection()
  {
    clear();
  }

#ifndef NO_VIRTUAL_COV
  SubstitutionProcessCollection*
#else
  Clonable*
#endif
  clone() const { return new SubstitutionProcessCollection(*this); }


  /**
   * @brief Resets all the information contained in this object.
   *
   */
   
  void clear();

  /**
   * @brief Add a new parametrizable to the matching collection with a
   * given number.
   *
   * @throw Exception if the number is already used. See replace
   * function instead.
   *
   * @param parametrizable A pointer toward a parametrizable, that will added to
   * the collection.
   * 
   * WARNING! The collection will now be the owner of the pointer, and will destroy it if needed!
   * Copy the parametrizable first if you don't want it to be lost!
   
   * @param parametrizableIndex The number of the parametrizable in
   * the Collection
   * @param withParameters boolean if the parameters of the object
   *         should be added (default: true)
   * 
   */

  void addParametrizable(Parametrizable* parametrizable, size_t parametrizableIndex, bool withParameters = true);

  /*
   * @brief specific methods to add specific objects.
   *
   */


  /**
   * ? operator[] ?
   */
  
  void addModel(SubstitutionModel* model, size_t modelIndex)
  {
    addParametrizable(model, modelIndex);
  }

  void addFrequencies(FrequenciesSet* frequencies, size_t frequenciesIndex)
  {
    addParametrizable(frequencies, frequenciesIndex);
  }

  void addDistribution(DiscreteDistribution* distribution, size_t distributionIndex)
  {
    addParametrizable(distribution, distributionIndex, (distributionIndex<10000));

    if (distributionIndex>=10000)
      mVConstDist_[distributionIndex/10000-1].push_back(distributionIndex%10000);
  }
  
  void addTree(ParametrizableTree* tree, size_t treeIndex)
  {
    addParametrizable(tree , treeIndex);
  }

  /**
   * @brief Get a SubstitutionModel from the collection.
   *
   * @param modelIndex The index of the model in the collection.
   * @return the got SubstitutionModel*. 
   */
  
  SubstitutionModel& getModel(size_t modelIndex)
  {
    return *(dynamic_cast<SubstitutionModel*>(modelColl_[modelIndex]));
  }

  const SubstitutionModel& getModel(size_t modelIndex) const
  {
    return *(dynamic_cast<const SubstitutionModel*>(modelColl_[modelIndex]));
  }

  /**
   * @brief Get a FrequenciesSet from the collection.
   *
   * @param frequenciesIndex The index of the frequencies set in the collection.
   * @return the got FrequenciesSet*. 
   */
  
  FrequenciesSet& getFrequencies(size_t frequenciesIndex)
  {
    return *(dynamic_cast<FrequenciesSet*>(freqColl_[frequenciesIndex]));
  }

  const FrequenciesSet& getFrequencies(size_t frequenciesIndex) const
  {
    return *(dynamic_cast<const FrequenciesSet*>(freqColl_[frequenciesIndex]));
  }

  /**
   * @brief Get a DiscreteDistribution from the collection.
   *
   * @param distributionIndex The index of the distribution in the collection.
   * @return the got DiscreteDistribution*. 
   */
  
  DiscreteDistribution& getRateDistribution(size_t distributionIndex)
  {
    return *(dynamic_cast<DiscreteDistribution*>(distColl_[distributionIndex]));
  }

  const DiscreteDistribution& getRateDistribution(size_t distributionIndex) const 
  {
    return *(dynamic_cast<const DiscreteDistribution*>(distColl_[distributionIndex]));
  }

  /**
   * @brief Get a tree from the set.
   *
   * @param treeIndex The index of the model in the set.
   * @return the got ParametrizableTree*. 
   */
  
  ParametrizableTree& getTree(size_t treeIndex)
  {
    return *(dynamic_cast<ParametrizableTree*>(treeColl_[treeIndex]));
  }
  
  const ParametrizableTree& getTree(size_t treeIndex) const 
  {
    return *(dynamic_cast<const ParametrizableTree*>(treeColl_[treeIndex]));
  }
  

  /**
   * @brief Get the numbers of the specified objects from the collections.
   *
   */
  
  std::vector<size_t> getModelNumbers() const
  {
    return modelColl_.keys();
  }

  std::vector<size_t> getFrequenciesNumbers() const
  {
    return freqColl_.keys();
  }

  std::vector<size_t> getRateDistributionNumbers() const 
  {
    return distColl_.keys();
  }

  std::vector<size_t> getTreeNumbers() const
  {
    return treeColl_.keys();
  }

  std::vector<size_t> getSubstitutionProcessNumbers() const
  {
    std::vector<size_t> vn;
    
    std::map<size_t, SubstitutionProcessCollectionMember*>::const_iterator it;
    for (it = mSubProcess_.begin(); it != mSubProcess_.end(); it++)
      vn.push_back(it->first);

    return vn;
  }

  bool hasModelNumber(size_t n) const
  {
    return modelColl_.hasObject(n);
  }

  bool hasFrequenciesNumber(size_t n) const
  {
    return freqColl_.hasObject(n);
  }

  bool hasDistributionNumber(size_t n) const 
  {
    return distColl_.hasObject(n);
  }

  bool hasTreeNumber(size_t n) const
  {
    return treeColl_.hasObject(n);
  }

  bool hasSubstitutionProcessNumber(size_t n) const 
  {
    return (mSubProcess_.find(n)!=mSubProcess_.end());
  }

  /**
   * @brief AbstractParameterAliasable functions, redirected towards
   * the process members.
   *
   * @{
   **/
  

  /**
   * @brief To be called when a parameter has changed. This will call
   * fireParameterChanged on the collections and the process members.
   *
   * @param parameters The modified parameters.
   */
  
  void fireParameterChanged(const ParameterList& parameters);

  void setNamespace(const std::string& prefix);

  void aliasParameters(const std::string& p1, const std::string& p2) throw (ParameterNotFoundException, Exception);

  void unaliasParameters(const std::string& p1, const std::string& p2) throw (ParameterNotFoundException, Exception);

  void aliasParameters(std::map<std::string, std::string>& unparsedParams, bool verbose) throw (ParameterNotFoundException, Exception);
  

  /**
   * @}
   **/
  
  /*
   * @brief Method to add a SubstitutionProcess.
   *
   * @param nProc the number of the Substitution process
   * @param mModBr a map associating numbers of models (from the collection) and numbers of branches
   * @param nTree the number of a Tree (from the collection)
   * @param nRate the number of a Distribution Rate (from the collection)
   * @param nFreq the number of a FrequenciesSet for the root (from the collection)
   *
   * @throw an Exception if the built SubstitutionModelSet is not complete or well built.
   *
   */

  void addSubstitutionProcess(size_t nProc, std::map<size_t, std::vector<int> > mModBr, size_t nTree, size_t nRate, size_t nFreq);

  /*
   * @brief Method to add a SubstitutionProcess.
   *
   * @param nProc the number of the Substitution process
   * @param mModBr a map associating numbers of models (from the collection) and numbers of branches
   * @param nTree the number of a Tree (from the collection)
   * @param nRate the number of a Distribution Rate (from the collection)
   *
   * @throw an Exception if the built SubstitutionModelSet is not complete or well built.
   *
   */

  void addSubstitutionProcess(size_t nProc, std::map<size_t, std::vector<int> > mModBr, size_t nTree, size_t nRate);

  /*
   * @brief Methods to retrieve Substitution Process
   *
   */

  size_t getNumberOfSubstitutionProcess() const { return mSubProcess_.size(); }

  SubstitutionProcess& getSubstitutionProcess(size_t  i)
  {
    return *dynamic_cast<SubstitutionProcess*>(mSubProcess_[i]);
  }

  const SubstitutionProcess& getSubstitutionProcess(size_t i) const
  {
    std::map<size_t, SubstitutionProcessCollectionMember*>::const_iterator it(mSubProcess_.find(i));
    
    return *dynamic_cast<const SubstitutionProcess*>(it->second);
  }

  
  /*
   * @brief Methods to retrieve parameters
   *
   */

  /*
   * @brief Get the branch lengths parameters.
   *
   * @return A ParameterList with all branch lengths.
   */
  
  ParameterList getBranchLengthParameters(bool independent) const
  {
    ParameterList pl=treeColl_.getParameters();

    if (independent)
      return pl.getCommonParametersWith(getIndependentParameters());
    else
      return pl;
  }

  ParameterList getBranchLengthParameters(size_t nTree, bool independent) const
  {
    ParameterList pl=treeColl_.getParametersForObject(nTree);
    
    if (independent)
      return pl.getCommonParametersWith(getIndependentParameters());
    else
      return pl;
  }
  
  bool hasBranchLengthParameter(const std::string& name) const;
  
  /**
   * @brief Get the parameters associated to substitution process(es).
   *
   * @return A ParameterList.
   */
  
  ParameterList getSubstitutionProcessParameters(size_t nProc, bool independent) const
  {
    ParameterList pl;
    
    if (mSubProcess_.find(nProc)!=mSubProcess_.end())
      pl=mSubProcess_.find(nProc)->second->getParameters();
    else
      return ParameterList();
    
    if (independent)
      return pl.getCommonParametersWith(getIndependentParameters());
    else
      return pl;
  }

  ParameterList getSubstitutionProcessParameters() const;

  /**
   * @brief Get the parameters associated to substitution model(s).
   *
   * @return A ParameterList.
   */
  
  ParameterList getSubstitutionModelParameters(size_t nMod, bool independent) const
  {
    ParameterList pl=modelColl_.getParametersForObject(nMod);

    if (independent)
      return pl.getCommonParametersWith(getIndependentParameters());
    else
      return pl;
  }

  ParameterList getSubstitutionModelParameters(bool independent) const
  {
    ParameterList pl=modelColl_.getParameters();

    if (independent)
      return pl.getCommonParametersWith(getIndependentParameters());
    else
      return pl;
  }
  

  /**
   * @brief Get the Non-derivable parameters
   *
   * @return A ParameterList.
   */
  
  ParameterList getNonDerivableParameters() const;
  
  /**
   * @brief Get the parameters associated to the rate distribution(s).
   *
   * @return A ParameterList.
   */
  
  ParameterList getRateDistributionParameters(size_t nRate, bool independent) const
  {
    ParameterList pl=distColl_.getParametersForObject(nRate);

    if (independent)
      return pl.getCommonParametersWith(getIndependentParameters());
    else
      return pl;
  }

  ParameterList getRateDistributionParameters(bool independent) const
  {
    ParameterList pl=distColl_.getParameters();

    if (independent)
      return pl.getCommonParametersWith(getIndependentParameters());
    else
      return pl;
  }
  
  /**
   * @brief Get the parameters associated to the root frequencies(s).
   *
   * @return A ParameterList.
   */
  
  ParameterList getRootFrequenciesParameters(size_t nFreq, bool independent) const
  {
    ParameterList pl=freqColl_.getParametersForObject(nFreq);

    if (independent)
      return pl.getCommonParametersWith(getIndependentParameters());
    else
      return pl;
  }

  ParameterList getRootFrequenciesParameters(bool independent) const
  {
    ParameterList pl=freqColl_.getParameters();

    if (independent)
      return pl.getCommonParametersWith(getIndependentParameters());
    else
      return pl;
  }
  

};
} // end of namespace bpp.

#endif // _SUBSTITUTIONPROCESSCOLLECTION_H_

