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

  /*
   * A collection of Frequencies Sets
   */
  
  ParametrizableCollection<FrequenciesSet> freqColl_;

  /*
   * A collection of DiscreteDistributions
   */
  
  ParametrizableCollection<DiscreteDistribution> distColl_;

  /*
   * A collection of trees
   *
   */

  ParametrizableCollection<ParametrizableTree> treeColl_;

  /*
   * A vector of SubstitutionProcessCollectionMember
   */

  std::vector<std::auto_ptr<SubstitutionProcessCollectionMember > > vSubProcess_;
  
public:
  /**
   * @brief Create empty collections.
   *
   * @param alpha The alphabet to use for this set.
   */
  
  SubstitutionProcessCollection():
    AbstractParameterAliasable(""),
    modelColl_(),
    freqColl_(),
    distColl_(),
    treeColl_(),
    vSubProcess_()
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
   
   * @param parametrizableIndex The number of the parametrizable in the Collection
   * 
   */

  void addParametrizable(Parametrizable* parametrizable, unsigned int parametrizableIndex);

  /*
   * @brief specific methods to add specific objects.
   *
   */
  
  void addModel(SubstitutionModel* model, unsigned int modelIndex)
  {
    addParametrizable(model, modelIndex);
  }

  void addFrequencies(FrequenciesSet* frequencies, unsigned int frequenciesIndex)
  {
    addParametrizable(frequencies, frequenciesIndex);
  }

  void addDistribution(DiscreteDistribution* distribution, unsigned int distributionIndex)
  {
    addParametrizable(distribution, distributionIndex);
  }
  
  void addTree(ParametrizableTree* tree, unsigned int treeIndex)
  {
    addParametrizable(tree , treeIndex);
  }

  /**
   * @brief Get a SubstitutionModel from the collection.
   *
   * @param modelIndex The index of the model in the collection.
   * @return the got SubstitutionModel*. 
   */
  
  SubstitutionModel* getModel(unsigned int modelIndex)
  {
    return (dynamic_cast<SubstitutionModel*>(modelColl_.getObject(modelIndex)));
  }

  const SubstitutionModel* getModel(unsigned int modelIndex) const
  {
    return (dynamic_cast<const SubstitutionModel*>(modelColl_.getObject(modelIndex)));
  }

  /**
   * @brief Get a FrequenciesSet from the collection.
   *
   * @param frequenciesIndex The index of the frequencies set in the collection.
   * @return the got FrequenciesSet*. 
   */
  
  FrequenciesSet* getFrequencies(unsigned int frequenciesIndex)
  {
    return (dynamic_cast<FrequenciesSet*>(freqColl_.getObject(frequenciesIndex)));
  }

  const FrequenciesSet* getFrequencies(unsigned int frequenciesIndex) const
  {
    return (dynamic_cast<const FrequenciesSet*>(freqColl_.getObject(frequenciesIndex)));
  }

  /**
   * @brief Get a DiscreteDistribution from the collection.
   *
   * @param distributionIndex The index of the distribution in the collection.
   * @return the got DiscreteDistribution*. 
   */
  
  DiscreteDistribution* getDistribution(unsigned int distributionIndex)
  {
    return (dynamic_cast<DiscreteDistribution*>(distColl_.getObject(distributionIndex)));
  }

  const DiscreteDistribution* getDistribution(unsigned int distributionIndex) const 
  {
    return (dynamic_cast<const DiscreteDistribution*>(distColl_.getObject(distributionIndex)));
  }

  /**
   * @brief Get a tree from the set.
   *
   * @param treeIndex The index of the model in the set.
   * @return the got ParametrizableTree*. 
   */
  
  ParametrizableTree* getTree(unsigned int treeIndex)
  {
    return (dynamic_cast<ParametrizableTree*>(treeColl_.getObject(treeIndex)));
  }
  
  const ParametrizableTree* getTree(unsigned int treeIndex) const 
  {
    return (dynamic_cast<ParametrizableTree*>(treeColl_.getObject(treeIndex)));
  }
  

  /**
   * @brief Get the numbers of the specified objects from the collections.
   *
   */
  
  std::vector<unsigned int> getModelNumbers() const
  {
    return modelColl_.keys();
  }

  std::vector<unsigned int> getFrequenciesNumbers() const
  {
    return freqColl_.keys();
  }

  std::vector<unsigned int> getDistributionNumbers() const 
  {
    return distColl_.keys();
  }

  std::vector<unsigned int> getTreeNumbers() const
  {
    return treeColl_.keys();
  }


  /**
   * @brief Remove a SubstitutionModel from the collection.
   *
   * @param modelIndex The index of the model in the collection.
   * @return the removed SubstitutionModel*. 
   */
  
  SubstitutionModel* removeModel(unsigned int modelIndex)
  {
    return (dynamic_cast<SubstitutionModel*>(modelColl_.removeObject(modelIndex)));
  }

  /**
   * @brief Remove a FrequenciesSet from the collection.
   *
   * @param frequenciesIndex The index of the frequencies set in the collection.
   * @return the removed FrequenciesSet*. 
   */
  
  FrequenciesSet* removeFrequencies(unsigned int frequenciesIndex)
  {
    return (dynamic_cast<FrequenciesSet*>(freqColl_.removeObject(frequenciesIndex)));
  }

  /**
   * @brief Remove a DiscreteDistribution from the collection.
   *
   * @param distributionIndex The index of the distribution in the collection.
   * @return the removed DiscreteDistribution*. 
   */
  
  DiscreteDistribution* removeDistribution(unsigned int distributionIndex)
  {
    return (dynamic_cast<DiscreteDistribution*>(distColl_.removeObject(distributionIndex)));
  }

  /**
   * @brief Remove a tree from the set.
   *
   * @param treeIndex The index of the model in the set.
   * @return the removed Tree*. 
   */
  
  ParametrizableTree* removeTree(unsigned int treeIndex)
  {
    return (dynamic_cast<ParametrizableTree*>(treeColl_.removeObject(treeIndex)));
  }
  
  /**
   * @brief Replace a parametrizable in the set, and returns the replaced one.
   *
   * @param parametrizableIndex The index of the model to be replaced in the set.
   * @param parametrizable the replacing Parametrizable
   * @return the replaced Parametrizable*. 
   */
  
  Parametrizable* replaceParametrizable(Parametrizable* parametrizable, unsigned int parametrizableIndex);

  /*
   * @brief specific methods to replace specific objects.
   *
   */
  
  SubstitutionModel* replaceModel(SubstitutionModel* model, unsigned int modelIndex)
  {
    return dynamic_cast<SubstitutionModel*>(replaceParametrizable(model, modelIndex));
  }

  FrequenciesSet* replaceFrequencies(FrequenciesSet* frequencies, unsigned int frequenciesIndex)
  {
    return dynamic_cast<FrequenciesSet*>(replaceParametrizable(frequencies, frequenciesIndex));
  }

  DiscreteDistribution* replaceDistribution(DiscreteDistribution* distribution, unsigned int distributionIndex)
  {
    return dynamic_cast<DiscreteDistribution*>(replaceParametrizable(distribution, distributionIndex));
  }
  
  ParametrizableTree* replaceTree(ParametrizableTree* tree, unsigned int treeIndex)
  {
    return dynamic_cast<ParametrizableTree*>(replaceParametrizable(tree, treeIndex));
  }

  /**
   * @brief To be called when a parameter has changed. This will call
   * fireParameterChanged on the collections.
   *
   * @param parameters The modified parameters.
   */
  
  void fireParameterChanged(const ParameterList& parameters);


  /*
   * Method to build a SubstitutionModelSet.
   *
   * @param mModBr a map associating numbers of models (from the collection) and numbers of branches
   * @param nTree the number of a Tree (from the collection)
   * @param nRate the number of a Distribution Rate (from the collection)
   * @param nFreq the number of a FrequenciesSet for the root (from the collection)
   *
   * @throw an Exception if the built SubstitutionModelSet is not complete or well built.
   *
   */

  void addSubstitutionProcess(std::map<std::vector<unsigned int> > mModBr, unsigned int nTree, unsigned int nRate, unsigned int nFreq);

  /*
   * Method to build a stationary SubstitutionModelSet.
   *
   * @param mModBr a map associating numbers of models (from the collection) and numbers of branches
   * @param nTree the number of a Tree (from the collection)
   * @param nRate the number of a Distribution Rate (from the collection)
   *
   * @throw an Exception if the built SubstitutionModelSet is not complete or well built.
   *
   */

  void addSubstitutionProcess(std::map<std::vector<unsigned int> > mModBr, unsigned int nTree, unsigned int nRate);

  
};
} // end of namespace bpp.

#endif // _SUBSTITUTIONPROCESSCOLLECTION_H_

