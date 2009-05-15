//
// File: ProteinSubstitutionModelWithFrequencies.h
// Created by: Julien Dutheil
// Created on: Mon Sept 01 13:30 2008
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

#ifndef _PROTEINSUBSTITUTIONMODELWITHFREQUENCIES_H_
#define _PROTEINSUBSTITUTIONMODELWITHFREQUENCIES_H_

#include "ProteinSubstitutionModel.h"
#include "FrequenciesSet.h"

//From SeqLib:
#include <Seq/SequenceContainerTools.h>

namespace bpp
{

/**
 * @brief Basal class for protein substitution model with estimated equilibrium frequencies (+F models).
 */
class ProteinSubstitutionModelWithFrequencies:
  public virtual ProteinSubstitutionModel
{
  protected:
    ProteinFrequenciesSet* freqSet_;

	public:
    /**
     * @brief Create a model with a given alphabet and frequencies set.
     *
     * @param alpha The alphabet to use.
     * @param freqSet The frequencies set object to use.
     */
		ProteinSubstitutionModelWithFrequencies(const ProteicAlphabet * alpha, const ProteinFrequenciesSet & freqSet, const string& prefix) :
      ProteinSubstitutionModel(alpha, prefix), freqSet_(dynamic_cast<ProteinFrequenciesSet *>(freqSet.clone()))
    {
      freq_ = freqSet_->getFrequencies();
      addParameters_(freqSet_->getParameters());
    }
	
    /**
     * @brief Create a model with a given alphabet.
     * A FullProteinFrequenciesSet object will be used in order to parametrize the equilibrium frequencies.
     *
     * @param alpha The alphabet to use.
     */
    ProteinSubstitutionModelWithFrequencies(const ProteicAlphabet * alpha, const string& prefix):
      ProteinSubstitutionModel(alpha, prefix)
    {
      freqSet_ = new FullProteinFrequenciesSet(alpha);
      freq_ = freqSet_->getFrequencies();
      addParameters_(freqSet_->getParameters());
    }
		
    ProteinSubstitutionModelWithFrequencies(const ProteinSubstitutionModelWithFrequencies & model):
      ProteinSubstitutionModel(model)
    {
      freqSet_ = dynamic_cast<ProteinFrequenciesSet *>(model.freqSet_->clone());
      freq_ = freqSet_->getFrequencies();
    }

    ProteinSubstitutionModelWithFrequencies & operator=(const ProteinSubstitutionModelWithFrequencies & model)
    {
      ProteinSubstitutionModel::operator=(model);
      freqSet_ = dynamic_cast<ProteinFrequenciesSet *>(model.freqSet_->clone());
      return *this;
    }

		virtual ~ProteinSubstitutionModelWithFrequencies()
    {
      delete freqSet_;
    }

#ifndef NO_VIRTUAL_COV
    ProteinSubstitutionModelWithFrequencies*
#else
    Clonable*
#endif
    clone() const = 0;
	
  public:
    void fireParameterChanged(const ParameterList & parameters)
    {
      freqSet_->matchParametersValues(parameters);
      freq_ = freqSet_->getFrequencies();
      ProteinSubstitutionModel::fireParameterChanged(parameters);
    }

    void setFrequenciesSet(const ProteinFrequenciesSet & freqSet)
    {
      delete freqSet_;
      freqSet_ = dynamic_cast<ProteinFrequenciesSet *>(freqSet.clone());
      resetParameters_();
      addParameters_(freqSet_->getParameters());
    }

    const ProteinFrequenciesSet & getFrequenciesSet() const { return *freqSet_; }

    void setFreqFromData(const SequenceContainer & data)
    {
      map<int, double> freqs = SequenceContainerTools::getFrequencies(data);
      double t = 0;
      for(unsigned int i = 0; i < size_; i++) t += freqs[i];
      for(unsigned int i = 0; i < size_; i++) freq_[i] = freqs[i] / t;
      freqSet_->setFrequencies(freq_);
      //Update parametrers and re-compute generator and eigen values:
      matchParametersValues(freqSet_->getParameters());
    }

};

} //end of namespace bpp.

#endif //_PROTEINSUBSTITUTIONMODELWITHFREQUENCIES_H_

