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
#include <Bpp/Numeric/Matrix/MatrixTools.h>

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
  pfitness_(pfitness),
  fixation_times_()
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
    double phi_i = fit?(*fit)[i]:1./(double)nbStates;
    // then for all couples ending with j
    for (size_t j=i+1;j<nbStates;j++)
    {
      double phi_j = fit?(*fit)[j]:1./(double)nbStates;
      
      // mutations
      generator_(i, nbloc) = nbAlleles * Q(i,j);
      generator_(j, nbloc + nbAlleles-2)= nbAlleles * Q(j,i);

      // drift + selection
      if (std::abs(generator_(i, nbloc)) < NumConstants::TINY() &&
          std::abs(generator_(j, nbloc + nbAlleles-2) < NumConstants::TINY())) // No mutation between states
      {
        for (size_t a=1;a<nbAlleles;a++)
        {
          generator_(nbloc+a-1 , (a>1?nbloc+a-2:i)) = 0;
          generator_(nbloc+a-1 , (a<nbAlleles-1?nbloc+a:j)) = 0;
        }
      }
      else
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

  double invnbStates=1./(double)nbStates;
  Vdouble pN(nbStates), pNm(nbStates), pNm2(nbStates), p(nbStates);
  for (size_t i=0;i<nbStates;i++)
  {
    pNm2[i]=fit?std::pow((*fit)[i],nbAlleles-2):invnbStates;
    pNm[i]=pNm2[i] * (fit? (*fit)[i]:invnbStates);
    pN[i]=pNm[i] * (fit? (*fit)[i]:invnbStates);
    freq_[i]=pFreq[i]*pNm[i];
  }

  size_t k=nbStates;
  for (size_t i = 0; i < nbStates-1; ++i)
  {
    double phi_i = fit?(*fit)[i]:invnbStates;
    // then for all couples ending with j
    for (size_t j=i+1;j<nbStates;j++)
    {
      double phi_j = fit?(*fit)[j]:invnbStates;

      auto mu=Q(i,j)*pFreq[i];   // alleles i & j

      double rfreq=pNm2[i];
      
      double rat=fit?(*fit)[j]/(*fit)[i]:1;
                            
      for (size_t n=1 ; n<nbAlleles; n++)   // i_{N-n}j_{n}
      {
        freq_[k++]=mu*rfreq*((int)(nbAlleles-n) * phi_i  + (int)n * phi_j)*nbAlleles/((int)n*(int)(nbAlleles-n));
        rfreq*=rat;
      }

    }
  }


  // normalizing stationary freq
  double x = VectorTools::sum(freq_);
  freq_ /= x;

  // Compute scaling factor
  double sden=0;
  
  double s = 0;
  for (size_t i = 0; i < nbStates-1; ++i)
  {
    double phi_i = fit?(*fit)[i]:invnbStates;
    // then for all couples ending with j
    for (size_t j=i+1;j<nbStates;j++)
    {
      double phi_j = fit?(*fit)[j]:invnbStates;

      double rat=phi_j/phi_i;
                            
      double fixtij = fixation_time(rat);
      double fixtji = fixation_time(1/rat);

      auto mij = 2 * Q(i,j) * pFreq[i] * pNm[i];
      auto mji = 2 * Q(j,i) * pFreq[j] * pNm[j];

      auto hij = rat==1.?1./(int)nbAlleles:(pN[j]/(pN[i]-pN[j])*(1-rat)/rat);
      auto hji = rat==1.?1./(int)nbAlleles:(pN[i]/(pN[j]-pN[i])*(rat-1));

      s += mij*hij/fixtij + mji*hji/fixtji;

      sden += Q(i,j) * pFreq[i] * (2 * pNm[i]
                                   + (phi_i+phi_j)* ((phi_i!=phi_j)?((pNm[i]-pNm[j])/(phi_i-phi_j)):(pNm2[i]*(nbAlleles-1))));
      sden += Q(j,i) * pFreq[j] * (2 * pNm[j]
                                   + (phi_i+phi_j)* ((phi_i!=phi_j)?((pNm[i]-pNm[j])/(phi_i-phi_j)):(pNm2[j]*(nbAlleles-1))));    

    }
  }

  // s is the probability of substitution on a duration of 1 generation (ie the actual scale time of the model).
  s/=sden;

  // And everything for exponential
  auto r=getScale();
  
  AbstractSubstitutionModel::updateMatrices();

  // Finally, set correct rate
  setScale(r/s);
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

double POMO::fixation_time(double fitness) const
{
  if (fixation_times_.find(fitness)==fixation_times_.end())
  {
    Vdouble pN(nbAlleles_+1,1); // pN[i] = fitness^i

    if (fitness!=1)
    {
      pN[0]=1.;
      for (size_t i=1;i<(size_t)nbAlleles_;i++)
      {
        pN[i]=fitness * pN[i-1];
      }
    }

    auto f = std::make_shared<vector<double>>(nbAlleles_);
    auto g = std::make_shared<vector<double>>(nbAlleles_);

    for (size_t i=0;i<nbAlleles_;i++)
      (*f)[i]=double(2*nbAlleles_-i-1);

    size_t cpt=0;

    while(cpt++<1000)
    {
      (*g)[nbAlleles_-1]=0;
      (*g)[0]=(1+(nbAlleles_-1)*fitness)/((nbAlleles_-1)*(1+fitness)) + (*f)[1];
      for (size_t i=1;i<nbAlleles_-1;i++)
      {
        auto n = double(i+1);
        if (fitness!=1)
        {
          double den = (1-pN[i+1])*(1+fitness);
          (*g)[i]=(n+((double)nbAlleles_-n)*fitness)/(n*((double)nbAlleles_-n)*(1+fitness))
            + (1-pN[i+2])/den * (*f)[i+1]
            + (fitness-pN[i+1])/den * (*f)[i-1];
        }
        else
        {
          (*g)[i]=(double)nbAlleles_/(2*n*((double)nbAlleles_-n))
            + (n+1)/(2*n) * (*f)[i+1]
            + (n-1)/(2*n) * (*f)[i-1];
        }
      }
      auto h=f;
      f=g;
      g=h;      

      bool flag=true;
      for (size_t i=0;i<f->size();i++)
        if (abs((*f)[i]-(*g)[i])> NumConstants::NANO())
        {
          flag=false;          
          break;
        }
      if (flag)
        break;
    }
    
    fixation_times_[fitness]=(*f)[0];
  }
  
  return (fixation_times_[fitness]);
}
