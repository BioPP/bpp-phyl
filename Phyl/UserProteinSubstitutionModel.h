//
// File: UserProteinSubstitutionModel.h
// Created by: Julien Dutheil
// Created on: Wed Aug 26 16:27 2005
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

#ifndef _USERPROTEINSUBSTITUTIONMODEL_H_
#define _USERPROTEINSUBSTITUTIONMODEL_H_

#include "ProteinSubstitutionModel.h"

// From SeqLib:
#include <Seq/ProteicAlphabet.h>

// From the STL:
#include <string>

using namespace std;

namespace bpp
{

/**
 * @brief Build an empirical protein substitution model from a file.
 * 
 * The file must follow PAML's format, and contain the exchangeabilities components (\f$S_{i,j}\f$)
 * and all equilibrium frequencies (\f$\pi_{i}\f$).
 * The generator is build so that \f$Q_{i,j} = \pi_i . S_{i,j}\f$, and is normalized
 * so that \f$\sum_i Q_{i,i} \times \pi_i = -1\f$.
 */
class UserProteinSubstitutionModel:
  public ProteinSubstitutionModel
{
  protected:
    const string _path;
  
  public:
    UserProteinSubstitutionModel(const ProteicAlphabet * alpha, const string & path, const string& prefix);

    virtual ~UserProteinSubstitutionModel() {}

#ifndef NO_VIRTUAL_COV
    UserProteinSubstitutionModel*
#else
    Clonable*
#endif
    clone() const { return new UserProteinSubstitutionModel(*this); }
      
  public:
    string getName() const;
    const string & getPath() const { return _path; }

  protected:
    void readFromFile();

};

} //end of namespace bpp.

#endif //_USERPROTEINSUBSTITUTIONMODEL_H_

