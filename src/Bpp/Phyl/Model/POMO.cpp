//
// File: POMO.cpp
// Authors:
//   Laurent Gueguen
// Created: jeudi 23 décembre 2021, à 15h 36
//

/*
  Copyright or ÃÂ© or Copr. CNRS, (November 16, 2004)
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/


#include "POMO.h"

using namespace bpp;

using namespace std;

/******************************************************************************/

POMO::POMO(const AllelicAlphabet* allAlph,
           std::shared_ptr<SubstitutionModel> pmodel,
           std::shared_ptr<FrequencySet> pfitness):
  AbstractParameterAliasable("POMO."),
  AbstractSubstitutionModel(allAlph, std::shared_ptr<const StateMap>(new CanonicalStateMap(allAlph, false)), "POMO."),
  nbAlleles_(allAlph->getNbAlleles()),
  pmodel_(pmodel),
  pfitness_(pfitness)
{

  const auto& alph=allAlph->getStateAlphabet();
  
  if (alph.getAlphabetType() != pmodel_->getAlphabet()->getAlphabetType())
    throw AlphabetMismatchException("POMO mismatch state alphabet for model.", &alph, pmodel_->getAlphabet());

  if (pfitness_ && alph.getAlphabetType() != pfitness_->getAlphabet()->getAlphabetType())
    throw AlphabetMismatchException("POMO mismatch state alphabet for fitness.", &alph, pfitness_->getAlphabet());

  if (pfitness_)
    pfitness_->setNamespace("POMO.fit_" + pfitness_->getNamespace());
  pmodel_->setNamespace("POMO." + pmodel_->getNamespace());

  pmodel_->enableEigenDecomposition(pmodel_->computeFrequencies()); // enable eigen if needed for pmodel->freq_
  if (pfitness_)
    addParameters_(pfitness_->getParameters());
  addParameters_(pmodel_->getParameters());

  computeFrequencies(false); // freq_ analytically defined
  updateMatrices();
}

void POMO::updateMatrices()
{
  auto nbStates=pmodel_->getNumberOfStates();
  auto nbAlleles=dynamic_cast<const AllelicAlphabet*>(getAlphabet())->getNbAlleles();

  const auto& Q=pmodel_->getGenerator();

  const Vdouble* fit=pfitness_?&pfitness_->getFrequencies():0;
  
  // Per couple of alleles

  // position of the bloc of alleles 
  size_t nbloc = nbStates;

  // for all couples starting with i
  for (size_t i=0;i<nbStates;i++)
  {
    double phi_i = fit?(*fit)[i]:0.25;
    // then for all couples ending with j
    for (size_t j=i+1;j<nbStates;j++)
    {
      double phi_j = fit?(*fit)[j]:0.25;
      
      // mutations
      generator_(i, nbloc) = nbAlleles * Q(i,j);
      generator_(j, nbloc + nbAlleles-2)= nbAlleles * Q(j,i);

      // drift + selection
      for (size_t a=1;a<nbAlleles;a++)
      {
        double rap = (double)(a * (nbAlleles - a ))/((double)a * phi_i  + (double)(nbAlleles - a) * phi_j);
        
        // drift towards i:  a -> a+1
        generator_(nbloc+a-1 , (a>1?nbloc+a-2:i)) = phi_i * rap;

        // drift towards j: nbAlleles - a -> nbAlleles -a +1
        generator_(nbloc+a-1 , (a<nbAlleles-1?nbloc+a:j)) = phi_j  * rap;
      }
      // setting for the next couple of alleles
      nbloc += nbAlleles-1;
    }
  }

  setDiagonal();

  // Stationnary distribution
  // and variables used for rate of substitutions
  
  const auto& pFreq = pmodel_->getFrequencies();

  Vdouble pN(nbStates), pNm(nbStates), pNm2(nbStates), p(nbStates);
  for (size_t i=0;i<nbStates;i++)
  {
    pNm2[i]=fit?std::pow((*fit)[i],nbAlleles-2):1;
    pNm[i]=fit?pNm2[i] * (*fit)[i]:1;
    pN[i]=fit?pNm[i] * (*fit)[i]:1;
    freq_[i]=pFreq[i]*pNm[i];
  }

  size_t k=nbStates;
  double sden=0;
  double snum=0;
  
  for (size_t i = 0; i < nbStates-1; ++i)
  {
    double phi_i = fit?(*fit)[i]:0.25;
    // then for all couples ending with j
    for (size_t j=i+1;j<nbStates;j++)
    {
      double phi_j = fit?(*fit)[j]:0.25;

      auto mu=Q(i,j)*pFreq[i];   // alleles i & j

      double rfreq=pNm2[i];
      
      double rat=fit?(*fit)[j]/(*fit)[i]:1;
                            
      for (size_t n=1 ; n<nbAlleles; n++)   // i_{N-n}j_{n}
      {
        freq_[k++]=mu*rfreq*((int)(nbAlleles-n) * phi_i  + (int)n * phi_j)*nbAlleles/((int)n*(int)(nbAlleles-n));
        rfreq*=rat;
      }

      snum+=2*mu*pNm2[i] * pN[j] * (pN[i]!=pN[j]?(phi_i-phi_j)/(pN[i]-pN[j]):(1./nbAlleles));
      sden+=mu * 2 * pNm[i] + ((phi_i+phi_j)*
                               ((phi_i!=phi_j)?((pNm[i]-pNm[j])/(phi_i-phi_j)):(nbAlleles-1)));
    }
  }
  
  // stationary freq
  double x = VectorTools::sum(freq_);
  freq_ /= x;


  // Specific normalization in numbers of substitutions (from appendix D of Genetics. 2019 Aug; 212(4): 1321–1336.)
  // s is the probability of substitution on a duration of 1 generation (ie the actual scale time of the model). 
  double s=snum/sden;
  
  // And everything for exponential
  AbstractSubstitutionModel::updateMatrices();

  setScale(1/s);
}
  
  
void POMO::fireParameterChanged(const ParameterList& parameters)
{
  if (pfitness_)
    pfitness_->matchParametersValues(parameters);
  pmodel_->matchParametersValues(parameters);

  updateMatrices();
}


void POMO::setFreq(map<int, double>& frequencies)
{
  // pfreqset_->setFrequenciesFromAlphabetStatesFrequencies(frequencies);
  // matchParametersValues(pfreqset_->getParameters());
}

