//
// File: OneChangeRegisterTransitionModel.cpp
// Created by: Laurent Gueguen
// Created on: samedi 24 octobre 2015, � 18h 50
//

/*
   Copyright or � or Copr. Bio++ Development Team, (November 16, 2004)

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

#include "OneChangeRegisterTransitionModel.h"
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Seq/Alphabet/IntegerAlphabet.h>

using namespace bpp;
using namespace std;

OneChangeRegisterTransitionModel::OneChangeRegisterTransitionModel(const SubstitutionModel& originalModel, const SubstitutionRegister& reg, size_t numReg) :
  AbstractParameterAliasable("OneChange."),
  AbstractFromSubstitutionModelTransitionModel(originalModel, "OneChange."),
  noChangedStates_(getNumberOfStates(), getNumberOfStates()),
  modelChanged_(),
  registerName_(reg.getName()),
  vNumRegs_(vector<size_t>(1,numReg))
{
  if ((numReg<=0) || (numReg>reg.getNumberOfSubstitutionTypes()))
    throw IndexOutOfBoundsException("OneChangeRegisterTransitionModel::OneChangeRegisterTransitionModel : wrong number for register category", numReg, 1, reg.getNumberOfSubstitutionTypes());

  // new alphabet to handle pit state
  shared_ptr<IntegerAlphabet> ialph(new IntegerAlphabet(getAlphabet()->getSize()));
    
  modelChanged_.reset(new AnonymousSubstitutionModel(ialph.get(), std::shared_ptr<const StateMap>(new CanonicalStateMap(ialph.get(),false))));
  modelChanged_->setScalable(false);
  
  for (size_t i=0;i<size_;i++)
  {
    vector<uint>& chS_i=noChangedStates_.getRow(i);
    
    for (size_t j=0;j<size_;j++)
      chS_i[j]=(reg.getType(i,j)!=numReg);
  }

  updateMatrices();
}

OneChangeRegisterTransitionModel::OneChangeRegisterTransitionModel(const SubstitutionModel& originalModel, const SubstitutionRegister& reg, vector<size_t> vNumRegs) :
  AbstractParameterAliasable("OneChange."),
  AbstractFromSubstitutionModelTransitionModel(originalModel, "OneChange."),
  noChangedStates_(getNumberOfStates(),getNumberOfStates()),
  modelChanged_(),
  registerName_(reg.getName()),
  vNumRegs_(vNumRegs)
{
  for (const auto& numReg : vNumRegs_)
    if ((numReg<=0) || (numReg>reg.getNumberOfSubstitutionTypes()))
      throw IndexOutOfBoundsException("OneChangeRegisterTransitionModel::OneChangeRegisterTransitionModel : wrong number for register category", numReg, 1, reg.getNumberOfSubstitutionTypes());

  shared_ptr<IntegerAlphabet> ialph(new IntegerAlphabet(getAlphabet()->getSize()));
    
  modelChanged_.reset(new AnonymousSubstitutionModel(ialph.get(), std::shared_ptr<const StateMap>(new CanonicalStateMap(ialph.get(),false))));
  modelChanged_->setScalable(false);

  for (size_t i=0;i<size_;i++)
  {
    vector<uint>& chS_i=noChangedStates_.getRow(i);

    for (size_t j=0;j<size_;j++)
      chS_i[j]=(std::find(vNumRegs_.begin(),vNumRegs_.end(),reg.getType(i,j))==vNumRegs_.end());
  }
  
  updateMatrices();  
}

/******************************************************************************/

void OneChangeRegisterTransitionModel::updateMatrices()
{
  const RowMatrix<double>& gen=getSubstitutionModel().getGenerator();
  
  for (size_t i = 0; i < size_; i++)
  {
    double si=0;
    const vector<double>& gen_i=gen.getRow(i);
    const vector<uint>& chS_i=noChangedStates_.getRow(i);
    
    for (size_t j = 0; j < size_; j++)
    {
      if (chS_i[j] || i==j)
        modelChanged_->setGenerator()(i,j)=gen_i[j];
      else
      {
        modelChanged_->setGenerator()(i,j)=0;
        si+=gen_i[j];
      }
    }
    modelChanged_->setGenerator()(i,size_)=si;
  }

  modelChanged_->updateMatrices();
}

/******************************************************************************/

double OneChangeRegisterTransitionModel::Pij_t(size_t i, size_t j, double t) const
{
  double ch_ijt=modelChanged_->Pij_t(i,j,t);  
  double ch_ist=modelChanged_->Pij_t(i,size_,t);  

  if (t==0)
  {
    const Matrix<double>& gen=getSubstitutionModel().getGenerator();
    if (ch_ist!=0)
      return (gen(i,j)-ch_ijt)/ch_ist;
    else
      return (i==j)?1:0;
  }

  if (ch_ist==0) // no proba from i to pit
    return getModel().Pij_t(i,j,t);
  else
    return (getModel().Pij_t(i,j,t)-ch_ijt)/ch_ist;
}

double OneChangeRegisterTransitionModel::dPij_dt(size_t i, size_t j, double t) const
{
  const Matrix<double>& ch_t=modelChanged_->getPij_t(t);

  if (t==0)
  {
    const RowMatrix<double>& Qch=modelChanged_->getGenerator();
    double si=Qch(i,size_);

    if (si==0)
      return 0;

    RowMatrix<double> Qch2;

    MatrixTools::mult<double>(Qch,Qch,Qch2);
    double dsi=Qch2(i,size_);

    const RowMatrix<double>& Q=getSubstitutionModel().getGenerator();
    double q2ij(0);

    for (size_t k = 0; k < size_; ++k)
      q2ij+=Q(i,k)*Q(k,j);
    
    return ((-q2ij+Qch2(i,j))/si + (Q(i,j)-Qch(i,j))*dsi/(si*si))/2;
  }
  
  double si=ch_t(i,size_);
  if (si==0)
    return 0;

  const Matrix<double>& dch_dt=modelChanged_->getdPij_dt(t);

  return ((getModel().dPij_dt(i,j,t)-dch_dt(i,j))*si-dch_dt(i,size_)*(getModel().Pij_t(i,j,t)-ch_t(i,j)))/(si*si);
}

double OneChangeRegisterTransitionModel::d2Pij_dt2(size_t i, size_t j, double t) const
{
  if (t==0)
  {
    const RowMatrix<double>& Q=getSubstitutionModel().getGenerator();
    const RowMatrix<double>& Qch=modelChanged_->getGenerator();
    double si=Qch(i,size_);

    if (si==0)
      return 0;

    RowMatrix<double> Q2;
    MatrixTools::mult<double>(Q,Q,Q2);

    RowMatrix<double> Qch2, Qch3;
    MatrixTools::mult<double>(Qch,Qch,Qch2);
    MatrixTools::mult<double>(Qch,Qch2,Qch3);
    
    double dsi=Qch2(i,size_);
    double d2si=Qch3(i,size_);

    double q3ij(0);

    for (size_t k = 0; k < size_; ++k)
      q3ij+=Q2(i,k)*Q(k,j);
    
    return ((Q(i,j)-Qch(i,j))*(si*d2si/3-dsi*dsi/2)+(Q2(i,j)-Qch2(i,j))*si*dsi/2-(q3ij-Qch3(i,j))*si*si/3)/(si*si*si);
  }
  
  const Matrix<double>& ch_t=modelChanged_->getPij_t(t);
  const Matrix<double>& dch_dt=modelChanged_->getdPij_dt(t);
  const Matrix<double>& d2ch_dt2=modelChanged_->getd2Pij_dt2(t);
  
  double d2u=getModel().d2Pij_dt2(i,j,t)-d2ch_dt2(i,j);
  double du=getModel().dPij_dt(i,j,t)-dch_dt(i,j);
  double u=getModel().Pij_t(i,j,t)-ch_t(i,j);
  
  double si=ch_t(i,size_);
  if (si==0)
    return 0;

  double si2=si*si;
  double dsi=dch_dt(i,size_);
  double d2si=d2ch_dt2(i,size_);
  
  return (d2u*si-d2si*u)/si-2*dsi*du/si2+2*dsi*dsi*u/(si2*si);
}


const Matrix<double>& OneChangeRegisterTransitionModel::getPij_t(double t) const
{
  if (t==0)
  {
    const RowMatrix<double>& Q=getSubstitutionModel().getGenerator();
    const RowMatrix<double>& Qch=modelChanged_->getGenerator();

    for (size_t i=0;i<size_;i++)
    {
      vector<double>& pi_t=pij_t.getRow(i);

      double si=Qch(i,size_);
      if (si!=0)
      {
        const vector<double>& qi=Q.getRow(i);
        const vector<double>& qchi=Qch.getRow(i);

        for (size_t j=0; j<size_;j++)
          pi_t[j]=(qi[j]-qchi[j])/si;
      }
      else
        for (size_t j=0; j<size_;j++)
          pi_t[j]= (i==j)?1:0;
    }
    return pij_t;
  }
  
  const RowMatrix<double>& orig_t=getModel().getPij_t(t);
  const RowMatrix<double>& ch_t=modelChanged_->getPij_t(t);
  
  for (unsigned int i = 0; i < size_; ++i) {
    vector<double>& pi_t=pij_t.getRow(i);
    const vector<double>& origi_t=orig_t.getRow(i);
    const vector<double>& chi_t=ch_t.getRow(i);

    double si=chi_t[size_];
    if (si==0)
      for (auto& x : pi_t)
        x=0;
    else      
      for (unsigned int j = 0; j < size_; ++j) 
        pi_t[j]=(origi_t[j]-chi_t[j])/si;
  }

  return pij_t;
}


const Matrix<double>& OneChangeRegisterTransitionModel::getdPij_dt(double t) const
{
  if (t==0)
  {
    const RowMatrix<double>& Q=getSubstitutionModel().getGenerator();
    const RowMatrix<double>& Qch=modelChanged_->getGenerator();
    
    RowMatrix<double> Qch2;
    MatrixTools::mult<double>(Qch,Qch,Qch2);
    
    RowMatrix<double> Q2;
    MatrixTools::mult<double>(Q,Q,Q2);
    
    for (size_t i=0;i<size_;i++)
    {
      vector<double>& dpi_t=dpij_t.getRow(i);
      
      const vector<double>& qchi=Qch.getRow(i);
      double si=qchi[size_];
      
      if (si!=0)
      {
        const vector<double>& qi=Q.getRow(i);
        const vector<double>& q2i=Q2.getRow(i);
        const vector<double>& q2chi=Qch2.getRow(i);
        double dsi=q2chi[size_];
      
        for (size_t j=0; j<size_;j++)
          dpi_t[j]= ((-q2i[j]+q2chi[j])/si + (qi[j]-qchi[j])*dsi/(si*si))/2;
      }
      else
        for (auto& x : dpi_t)
          x=0;
    }
    return dpij_t;
  }
  
  const RowMatrix<double>& orig_t=getModel().getPij_t(t);
  const RowMatrix<double>& ch_t=modelChanged_->getPij_t(t);
  const RowMatrix<double>& dorig_dt=getModel().getdPij_dt(t);
  const RowMatrix<double>& dch_dt=modelChanged_->getdPij_dt(t);

  for (unsigned int i = 0; i < size_; ++i) {
    vector<double>& dpi_dt=dpij_t.getRow(i);
    const vector<double>& chi_t=ch_t.getRow(i);
    
    double si=chi_t[size_];
    if (si==0)
    {
      for (auto& x : dpi_dt)
        x=0;
      continue;
    }
    
    const vector<double>& origi_t=orig_t.getRow(i);
    const vector<double>& dorigi_dt=dorig_dt.getRow(i);
    const vector<double>& dchi_dt=dch_dt.getRow(i);

    double dsi=dchi_dt[size_];

    for (unsigned int j = 0; j < size_; ++j) 
      dpi_dt[j]= ((dorigi_dt[j]-dchi_dt[j])*si-dsi*(origi_t[j]-chi_t[j]))/(si*si);
  }
  
  return dpij_t;
}


const Matrix<double>& OneChangeRegisterTransitionModel::getd2Pij_dt2(double t) const
{
  if (t==0)
  {
    const RowMatrix<double>& Q=getSubstitutionModel().getGenerator();
    const RowMatrix<double>& Qch=modelChanged_->getGenerator();
    
    RowMatrix<double> Qch2, Qch3;
    MatrixTools::mult<double>(Qch,Qch,Qch2);
    MatrixTools::mult<double>(Qch,Qch2,Qch3);
    
    RowMatrix<double> Q2, Q3;
    MatrixTools::mult<double>(Q,Q,Q2);
    MatrixTools::mult<double>(Q,Q2,Q3);
    
    for (size_t i=0;i<size_;i++)
    {
      vector<double>& d2pi_t=d2pij_t.getRow(i);
      
      const vector<double>& qchi=Qch.getRow(i);
      double si=qchi[size_];
      if (si!=0)
      {
        const vector<double>& qi=Q.getRow(i);
        const vector<double>& q2i=Q2.getRow(i);
        const vector<double>& q2chi=Qch2.getRow(i);
        const vector<double>& q3i=Q3.getRow(i);
        const vector<double>& q3chi=Qch3.getRow(i);
        double dsi=Qch2(i,size_);
        double d2si=Qch3(i,size_);
        
        for (size_t j=0; j<size_;j++)
          d2pi_t[j]= -((qi[j]-qchi[j])*(si*d2si/3-dsi*dsi/2)+(q2i[j]-q2chi[j])*si*dsi/2-(q3i[j]-q3chi[j])*si*si/3)/(si*si*si);
      }
      else
        for (auto& x : d2pi_t)
          x=0;
    }
    return d2pij_t;
  }
  
  const RowMatrix<double>& orig_t=getModel().getPij_t(t);
  const RowMatrix<double>& ch_t=modelChanged_->getPij_t(t);
  const RowMatrix<double>& dorig_dt=getModel().getdPij_dt(t);
  const RowMatrix<double>& dch_dt=modelChanged_->getdPij_dt(t);
  const RowMatrix<double>& d2orig_dt2=getModel().getd2Pij_dt2(t);
  const RowMatrix<double>& d2ch_dt2=modelChanged_->getd2Pij_dt2(t);

  for (unsigned int i = 0; i < size_; ++i) {
    vector<double>& d2pi_dt2=d2pij_t.getRow(i);
    const vector<double>& chi_t=ch_t.getRow(i);

    double si=chi_t[size_];
    if (si==0)
    {
      for (auto& x : d2pi_dt2)
        x=0;
      continue;
    }

    const vector<double>& origi_t=orig_t.getRow(i);
    const vector<double>& dorigi_dt=dorig_dt.getRow(i);
    const vector<double>& dchi_dt=dch_dt.getRow(i);
    const vector<double>& d2origi_dt2=d2orig_dt2.getRow(i);
    const vector<double>& d2chi_dt2=d2ch_dt2.getRow(i);
    
    double dsi=dchi_dt[size_];
    double d2si=d2chi_dt2[size_];

    double si2=si*si;
  
    for (unsigned int j = 0; j < size_; ++j)
    {
      double d2u=d2origi_dt2[j]-d2chi_dt2[j];
      double du=dorigi_dt[j]-dchi_dt[j];
      double u=origi_t[j]-chi_t[j];
  
      d2pi_dt2[j]=-((d2u*si-d2si*u-2*dsi*du)/si2+2*dsi*dsi*u/(si2*si));
    }
  }
  
  return d2pij_t;
}


