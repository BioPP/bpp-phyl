//
// File: RegisteredParameter.h
// Created by: Laurent Guéguen
// Created on: lundi 15 avril 2013, à 15h 08
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004,
  2005, 2006).

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

#ifndef _REGISTEREDPARAMETER_H_
#define _REGISTEREDPARAMETER_H_

#include <Bpp/Numeric/Parameter.h>
#include "../../Mapping/SubstitutionRegister.h"

//From the STL:

namespace bpp
{

  /**
   * @brief In the case a parameter is linked with a substitution register,
   *  this virtual class defines all that is necessary for it.
   */
  
  class RegisteredParameter:
    public Parameter
  {
  protected:
    const Alphabet* alphabet_;
    
    std::auto_ptr<SubstitutionRegister> SubReg_;

    // value used as a null hypothesis one
    double null_;

    void setSubstitutionRegister(SubstitutionRegister* pSubReg) {
      SubReg_.reset(pSubReg);
    }
    
  public:
    RegisteredParameter(const Alphabet* alpha): Parameter(), alphabet_(alpha), SubReg_(0), null_(1) {}

    RegisteredParameter(const Alphabet* alpha, const std::string& name, double value, Constraint* constraint, bool attachConstraint, double precision=0) : Parameter(name, value, constraint, attachConstraint, precision), alphabet_(alpha), SubReg_(0), null_(1) {}

    RegisteredParameter(const Alphabet* alpha, const std::string& name, double value, const Constraint* constraint = 0, double precision=0) : Parameter(name, value, constraint, precision), alphabet_(alpha), SubReg_(0), null_(1) {}

    RegisteredParameter(const RegisteredParameter& rp): Parameter(rp),
                                                        alphabet_(rp.alphabet_),
                                                        SubReg_(rp.SubReg_.get()->clone()),
                                                        null_(rp.null_) {}

    RegisteredParameter& operator=(const RegisteredParameter& rp)
    {
      Parameter::operator=(rp);
      alphabet_=rp.alphabet_;
      SubReg_=std::auto_ptr<SubstitutionRegister>(rp.SubReg_.get()->clone());
      null_=rp.null_;
      return *this;
    }

    virtual ~RegisteredParameter() {};
    
#ifndef NO_VIRTUAL_COV
    virtual RegisteredParameter* clone() const = 0;
#endif

    
  public:

    const Alphabet* getAlphabet() const
    {
      return alphabet_;
    }

    /**
     *@brief Returns the value for the null model.
     *
     */
    
    double getNull() const
    {
      return null_;
    }
    
    /**
     * @brief Method to compute an estimate of the value of the
     * parameter from values assigned to each register.
     * 
     * @param A vector of the values of the registers
     * @return The estimate parameter value
     */
    
    virtual double computeEstimate(std::vector<double> values) const = 0;

    const SubstitutionRegister& getSubstitutionRegister() const
    {
      return *SubReg_;
    }

  };

}// end of namespace bpp

#endif //_REGISTEREDPARAMETER_H_

