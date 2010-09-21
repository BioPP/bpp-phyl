//
// File: AbstractWordReversibleSubstitutionModel.cpp
// Created by:  Laurent Gueguen
// Created on: Jan 2009
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

#include "AbstractWordReversibleSubstitutionModel.h"

// From SeqLib:
#include <Seq/SequenceContainerTools.h>
#include <Seq/WordAlphabet.h>
#include <Seq/AlphabetTools.h>

// From NumCalc:
#include <NumCalc/MatrixTools.h>
#include <NumCalc/VectorTools.h>
#include <NumCalc/EigenValue.h>

using namespace bpp;

// From the STL:
#include <cmath>

using namespace std;

/******************************************************************************/

AbstractWordReversibleSubstitutionModel::AbstractWordReversibleSubstitutionModel(
  const std::vector<SubstitutionModel*>& modelVector,
  const std::string& st) :
  AbstractReversibleSubstitutionModel(AbstractWordReversibleSubstitutionModel::extractAlph(modelVector), st),
  new_alphabet_ (true),
  VSubMod_      (),
  VnestedPrefix_(),
  rate_         (modelVector.size()),
  p_            (getNumberOfStates(), getNumberOfStates())
{
  enableEigenDecomposition(false);
  unsigned int i, j;
  unsigned int n = modelVector.size();

  // test whether two models are identical

  bool flag = 0;
  i = 0;
  j = 1;
  while (!flag && i < (n - 1))
  {
    if (modelVector[i] == modelVector[j])
      flag = 1;
    else
    {
      j++;
      if (j == n)
      {
        i++;
        j = i + 1;
      }
    }
  }

  if (!flag)
  {
    for (i = 0; i < n; i++)
    {
      VSubMod_.push_back(modelVector[i]);
      VnestedPrefix_.push_back(modelVector[i]->getNamespace());
      VSubMod_[i]->setNamespace(st + TextTools::toString(i+1) + "_" + VnestedPrefix_[i]);
      addParameters_(VSubMod_[i]->getParameters());
    }
  }
  else
  {
   string t = "";
    for (i = 0; i < n; i++)
    {
      VSubMod_.push_back(modelVector[0]);
      VnestedPrefix_.push_back(modelVector[0]->getNamespace());
      t += TextTools::toString(i+1);
    }
    VSubMod_[0]->setNamespace(st + t + "_" + VnestedPrefix_[0]);
    addParameters_(VSubMod_[0]->getParameters());
  }

  for (i = 0; i < n; i++)
  {
    rate_[i] = 1.0 / n;
  }
}

AbstractWordReversibleSubstitutionModel::AbstractWordReversibleSubstitutionModel(
  const Alphabet* alph,
  const std::string& st) :
  AbstractReversibleSubstitutionModel(alph, st),
  new_alphabet_ (false),
  VSubMod_      (),
  VnestedPrefix_(),
  rate_         (0),
  p_            (getNumberOfStates(), getNumberOfStates())
{
  enableEigenDecomposition(false);
}

AbstractWordReversibleSubstitutionModel::AbstractWordReversibleSubstitutionModel(
  SubstitutionModel* pmodel,
  unsigned int num,
  const std::string& st) :
  AbstractReversibleSubstitutionModel(new WordAlphabet(pmodel->getAlphabet(), num),st),
  new_alphabet_ (true),
  VSubMod_      (),
  VnestedPrefix_(),
  rate_         (num),
  p_            (getNumberOfStates(), getNumberOfStates())
{
  enableEigenDecomposition(false);
  unsigned int i;

  string t = "";
  for (i = 0; i < num; i++)
  {
   VSubMod_.push_back(pmodel);
   VnestedPrefix_.push_back(pmodel->getNamespace());
    rate_[i] = 1.0 / num;
    t += TextTools::toString(i+1);
  }

  pmodel->setNamespace(st + t + "_" + VnestedPrefix_[0]);
  addParameters_(pmodel->getParameters());
}

AbstractWordReversibleSubstitutionModel::AbstractWordReversibleSubstitutionModel(
  const AbstractWordReversibleSubstitutionModel& wrsm) :
  AbstractReversibleSubstitutionModel(wrsm),
  new_alphabet_ (wrsm.new_alphabet_),
  VSubMod_      (),
  VnestedPrefix_(wrsm.VnestedPrefix_),
  rate_         (wrsm.rate_),
  p_            (wrsm.p_)
{
   unsigned int i;
   unsigned int num = wrsm.VSubMod_.size();

  if (wrsm.new_alphabet_)
    alphabet_ = new WordAlphabet(*(dynamic_cast<const WordAlphabet*>(wrsm.getAlphabet())));

  SubstitutionModel* pSM = 0;
  if ((num > 1) & (wrsm.VSubMod_[0] == wrsm.VSubMod_[1]))
    pSM = wrsm.VSubMod_[0]->clone();

  for (i = 0; i < num; i++)
  {
   VSubMod_.push_back(pSM ? pSM : wrsm.VSubMod_[i]->clone());
  }
}

AbstractWordReversibleSubstitutionModel& AbstractWordReversibleSubstitutionModel::operator=(
  const AbstractWordReversibleSubstitutionModel& wrsm)
{
  AbstractReversibleSubstitutionModel::operator=(wrsm);
  new_alphabet_  = wrsm.new_alphabet_;
  VnestedPrefix_ = wrsm.VnestedPrefix_;
  rate_          = wrsm.rate_;
  p_             = wrsm.p_;

  unsigned int i;
  unsigned int num = wrsm.VSubMod_.size();

  if (wrsm.new_alphabet_)
    alphabet_ = new WordAlphabet(*(dynamic_cast<const WordAlphabet*>(wrsm.getAlphabet())));

  SubstitutionModel* pSM = 0;
  if ((num > 1) & (wrsm.VSubMod_[0] == wrsm.VSubMod_[1]))
    pSM = wrsm.VSubMod_[0]->clone();

  for (i = 0; i < num; i++)
  {
    VSubMod_[i] =  (pSM ? pSM : wrsm.VSubMod_[i]->clone());
  }

  return *this;
}

AbstractWordReversibleSubstitutionModel::~AbstractWordReversibleSubstitutionModel()
{
  if ((VSubMod_.size() > 1) && (VSubMod_[0] == VSubMod_[1]))
  {
    if (VSubMod_[0])
      delete VSubMod_[0];
  }
  else
    for (unsigned int i = 0; i < VSubMod_.size(); i++)
    {
      if (VSubMod_[i])
        delete VSubMod_[i];
    }
  if (new_alphabet_)
    delete alphabet_;
}

unsigned int AbstractWordReversibleSubstitutionModel::getNumberOfStates() const
{
  return getAlphabet()->getSize();
}

Alphabet* AbstractWordReversibleSubstitutionModel::extractAlph(const vector<SubstitutionModel*>& modelVector)
{
   unsigned int i;

   vector<const Alphabet*> vAlph;

  for (i = 0; i < modelVector.size(); i++)
  {
   vAlph.push_back(modelVector[i]->getAlphabet());
  }

  return new WordAlphabet(vAlph);
}

void AbstractWordReversibleSubstitutionModel::setNamespace(const std::string& prefix)
{
   AbstractReversibleSubstitutionModel::setNamespace(prefix);

  if (VSubMod_.size() < 2 || VSubMod_[0] == VSubMod_[1])
  {
   string t = "";
    for (unsigned int i = 0; i < VSubMod_.size(); i++)
    {
      t += TextTools::toString(i+1);
    }
    VSubMod_[0]->setNamespace(prefix + t + "_" + VnestedPrefix_[0]);
  }
  else
  {
    for (unsigned int i = 0; i < VSubMod_.size(); i++)
    {
      VSubMod_[i]->setNamespace(prefix + TextTools::toString(i+1) + "_" + VnestedPrefix_[i]);
    }
  }
}

/******************************************************************************/

void AbstractWordReversibleSubstitutionModel::updateMatrices()
{
  //First we update position specific models. This need to be done here and not
  //in fireParameterChanged, has some parameter aliases might have been defined
  //and need to be resolved first.
  if (VSubMod_.size() < 2 || VSubMod_[0] == VSubMod_[1])
    VSubMod_[0]->matchParametersValues(getParameters());
  else
    for (unsigned int i = 0; i < VSubMod_.size(); i++)
    {
      VSubMod_[i]->matchParametersValues(getParameters());
    }

  unsigned int nbmod = VSubMod_.size();
  unsigned int salph = getNumberOfStates();

  // Generator

  if (enableEigenDecomposition())
  {
    unsigned int i, j, n, l, k, m;

    vector<unsigned int> vsize;

    for (k = 0; k < nbmod; k++)
    {
      vsize.push_back(VSubMod_[k]->getNumberOfStates());
    }

    RowMatrix<double> gk, exch;

    m = 1;

    for (k = nbmod; k > 0; k--)
    {
      gk = VSubMod_[k - 1]->getGenerator();
      exch = (dynamic_cast<AbstractReversibleSubstitutionModel*>(VSubMod_[k - 1]))->getExchangeabilityMatrix();
      for (i = 0; i < vsize[k - 1]; i++)
      {
        for (j = 0; j < vsize[k - 1]; j++)
        {
          if (i != j)
          {
            n = 0;
            while (n < salph)
            { // loop on prefix
              for (l = 0; l < m; l++)
              { // loop on suffix
                generator_(n + i * m + l, n + j * m + l) = gk(i,j) * rate_[k - 1];
                exchangeability_(n + i * m + l, n + j * m + l) = exch(i,j) * rate_[k - 1];
              }
              n += m * vsize[k - 1];
            }
          }
        }
      }
      m *= vsize[k - 1];
    }
  }

  // modification of generator_ and freq_

  completeMatrices();

  // at that point generator_ and freq_ are done for models without
  // enableEigenDecomposition

  // Eigen values:


  if (enableEigenDecomposition())
  {
    unsigned int i, j;
    double x;

    unsigned int nbStop;
    Vdouble vi;

    for (i = 0; i < salph; i++)
    {
      x = 0;
      for (j = 0; j < salph; j++)
      {
        if (j != i)
          x += generator_(i, j);
      }
      generator_(i, i) = -x;
    }

    //02/03/10 Julien: this should be avoided, we have to find a way to avoid particular cases like this...
    if (AlphabetTools::isCodonAlphabet(getAlphabet()))
    {
      int gi = 0, gj = 0;

      const CodonAlphabet* pca = dynamic_cast<const CodonAlphabet*>(getAlphabet());

      RowMatrix<double> gk;

      nbStop = pca->numberOfStopCodons();
      gk.resize(salph - nbStop, salph - nbStop);
      for (i = 0; i < salph; i++)
      {
        if (!pca->isStop(i))
        {
          gj = 0;
          for (j = 0; j < salph; j++)
          {
            if (!pca->isStop(j))
            {
              gk(i - gi, j - gj) = generator_(i, j);
            }
            else
              gj++;
          }
        }
        else
          gi++;
      }

      EigenValue<double> ev(gk);
      eigenValues_ = ev.getRealEigenValues();
      vi = ev.getImagEigenValues();

      for (i = 0; i < nbStop; i++)
      {
        eigenValues_.push_back(0);
      }

      RowMatrix<double> rev = ev.getV();
      rightEigenVectors_.resize(salph, salph);
      gi = 0;
      for (i = 0; i < salph; i++)
      {
        if (pca->isStop(i))
        {
          gi++;
          for (j = 0; j < salph; j++)
          {
            rightEigenVectors_(i,j) = 0;
          }
          rightEigenVectors_(i,salph - nbStop + gi - 1) = 1;
        }
        else
        {
          for (j = 0; j < salph - nbStop; j++)
          {
            rightEigenVectors_(i, j) = rev(i - gi,j);
          }
          for (j = salph - nbStop; j < salph; j++)
          {
            rightEigenVectors_(i,j) = 0;
          }
        }
      }
    }
    else
    {
      EigenValue<double> ev(generator_);
      eigenValues_ = ev.getRealEigenValues();
      vi = ev.getImagEigenValues();
      rightEigenVectors_ = ev.getV();
      nbStop = 0;
    }

    MatrixTools::inv(rightEigenVectors_, leftEigenVectors_);

    // looking for the 0 eigenvector

    unsigned int nulleigen = 0;
    while (nulleigen < salph - nbStop)
    {
      if (abs(eigenValues_[nulleigen]) < 0.000001 && abs(vi[nulleigen]) < 0.000001)
        break;
      else
        nulleigen++;
    }

    if (nulleigen >= salph - nbStop)
    {
      cerr << "AbstractWordReversibleSubstitutionModel::updateMatrices : Problem in eigenspace of " << getName() << endl;
      exit(0);
    }

    for (i = 0; i < salph; i++)
    {
      freq_[i] = leftEigenVectors_(nulleigen, i);
    }

    x = 0;
    for (i = 0; i < salph; i++)
    {
      x += freq_[i];
    }

    for (i = 0; i < salph; i++)
    {
      freq_[i] /= x;
    }

    // normalization

    x = 0;
    for (i = 0; i < salph; i++)
    {
      x += freq_[i] * generator_(i, i);
    }

    MatrixTools::scale(generator_, -1. / x);

    for (i = 0; i < salph; i++)
    {
      eigenValues_[i] /= -x;
    }
  }
}


void AbstractWordReversibleSubstitutionModel::setFreq(std::map<int, double>& freqs)
{
  map<int, double> tmpFreq;
  unsigned int nbmod = VSubMod_.size();

  unsigned int i, j, s, k, d, size;

  d = size = getNumberOfStates();

  for (i = 0; i < nbmod; i++)
  {
    tmpFreq.clear();
    s = VSubMod_[i]->getAlphabet()->getSize();
    d /= s;
    for (j = 0; j < s; j++)
    {
      tmpFreq[j] = 0;
    }
    for (k = 0; k < size; k++)
    {
      tmpFreq[(k / d) % s] += freqs[k];
    }
    VSubMod_[i]->setFreq(tmpFreq);
  }

  updateMatrices();
}
