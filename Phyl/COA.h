//
// File: COACOA.h
// Created by: Bastien Boussau
// Created on: Tue May 18 15:23:20 2010
//

/*
   Copyright or � or Copr. CNRS, (November 16, 2004)
   COA
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


#ifndef _COA_H_
#define _COA_H_

#include "ProteinSubstitutionModel.h"
#include "AbstractSubstitutionModel.h"


// From SeqLib:
#include <Seq/ProteicAlphabet.h>

using namespace std;


namespace bpp
{
/**
 * @brief A empirico-parametric amino-acid substitution model.
 *
 * An empirical exchangeability matrix is used, and equilibrium frequencies are estimated in an economic way, using hyperparameters obtained with a PCA analysis of
 * proteins from the three domains of life.
 * This way really is economic mainly in non-homogeneous cases, when several matrices are used on a single tree topology.
 * _parameters is a list of parameters ordered as follows: first the axes probabilities, and then the sets of 20 equilibrium frequencies corresponding to these axes.
 */

class COA :
  public ProteinSubstitutionModel,
  public AbstractReversibleSubstitutionModel
{
protected:
  string _baseModel;
  static IncludingInterval _WeightConstraint; // = IncludingInterval(-1.0, 1.0);
  RowMatrix<double>* _P;

public:
  COA(const ProteicAlphabet* alpha,
      const string baseModel = "JTT92",
      vector<double> axisWeights = vector<double>( 19, 0.0 ),
      vector<double> equilibriumFrequencies = vector<double>( 380, 0.05 ));

  COA(const COA& coa);

  COA& operator=(const COA& coa);

  virtual ~COA() {}

#ifndef NO_VIRTUAL_COV
  COA*
#else
  Clonable*
#endif
  clone() const { return new COA(*this); }

public:
  string getName() const { return "COA"; }

protected:
  void readFromFile();
  void normalizeProbabilities();
  void computeEquilibriumFrequencies();
  void updateMatrices();
  vector<double> prodMatrixVector(RowMatrix<double>& P, vector<double>& V);
};
} // end of namespace bpp.
#endif  // _COA_H_
