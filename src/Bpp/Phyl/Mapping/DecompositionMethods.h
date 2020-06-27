//
// File: DecompositionMethods.h
// Created by: Laurent Guéguen
// Created on: mardi 18 juillet 2017, à 22h 42
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004, 2005, 2006)

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

#ifndef _DECOMPOSITION_METHODS_H_
#define _DECOMPOSITION_METHODS_H_

#include <Bpp/Numeric/Matrix/Matrix.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>

#include "../Model/SubstitutionModel.h"
#include "SubstitutionRegister.h"

namespace bpp
{

/**
 * @brief Methods useful for analytical substitution count and rewards
 * using the eigen decomposition method.
 *
 * The code is adapted from the original R code by Paula Tataru and
 * Asger Hobolth.
 *
 * @author Julien Dutheil, Laurent Guéguen
 */

  class DecompositionMethods
  {
  protected:
    const SubstitutionModel* model_;
    size_t nbStates_;
    size_t nbTypes_;
    mutable RowMatrix<double> jMat_, jIMat_;

    /*
     * @brief Real and imaginary eigenvectors, for non-reversible
     * computation
     */
    
    ColMatrix<double> rightEigenVectors_, rightIEigenVectors_;
    RowMatrix<double> leftEigenVectors_, leftIEigenVectors_;

    /*
     * @brief computation matrices
     *
     */
    
    std::vector< RowMatrix<double> > bMatrices_, insideProducts_, insideIProducts_;
    
  public:
    DecompositionMethods(const SubstitutionModel* model, SubstitutionRegister* reg);

    DecompositionMethods(SubstitutionRegister* reg);
    
    DecompositionMethods(const SubstitutionModel* model);

    DecompositionMethods(const StateMap& statemap);

    DecompositionMethods(const DecompositionMethods& dm) :
      model_(dm.model_),
      nbStates_(dm.nbStates_),
      nbTypes_(dm.nbTypes_),
      jMat_(dm.jMat_),
      jIMat_(dm.jIMat_),
      rightEigenVectors_(dm.rightEigenVectors_),
      rightIEigenVectors_(dm.rightIEigenVectors_),
      leftEigenVectors_(dm.leftEigenVectors_),
      leftIEigenVectors_(dm.leftIEigenVectors_),
      bMatrices_(dm.bMatrices_),
      insideProducts_(dm.insideProducts_),
      insideIProducts_(dm.insideIProducts_)
    {}				
    
    DecompositionMethods& operator=(const DecompositionMethods& dm)
    {
      model_          = dm.model_;
      nbStates_       = dm.nbStates_;
      nbTypes_        = dm.nbTypes_;
      jMat_           = dm.jMat_;
      jIMat_          = dm.jIMat_;
      
      rightEigenVectors_ = dm.rightEigenVectors_;
      rightIEigenVectors_ = dm.rightIEigenVectors_;
      leftEigenVectors_ = dm.leftEigenVectors_;
      leftIEigenVectors_ = dm.leftIEigenVectors_;
      bMatrices_      = dm.bMatrices_;
      insideProducts_ = dm.insideProducts_;
      insideIProducts_ =  dm.insideIProducts_;
      
      return *this;
    }				
    
    DecompositionMethods* clone() const { return new DecompositionMethods(*this); }

    virtual ~DecompositionMethods() {};


    /**
     * @brief Set the substitution model.
     *
     * @param model A pointer toward the substitution model to use.
     */
    
    void setSubstitutionModel(const SubstitutionModel* model);


  protected:

    void initStates_();
    
    void initBMatrices_();

    void computeProducts_();

    /*
     * @brief Perform the computation of the conditional expectations
     *
     */

    void computeExpectations(RowMatrix<double>& mapping, double length) const;

    void computeExpectations(std::vector< RowMatrix<double> >& mappings, double leangth) const;

    /**
     * @brief Compute the integral part of the computation
     *
     */
    
    void jFunction_(const std::vector<double>& lambda, double t, RowMatrix<double>& result) const;

    /**
     * @brief Compute the integral part of the computation, in complex numbers
     *
     */
    
    void jFunction_(const std::vector<double>& lambda, const std::vector<double>& ilambda, double t, RowMatrix<double>& result, RowMatrix<double>& iresult) const;
    
  };

} //end of namespace bpp.

#endif // _DECOMPOSITION_METHODS_H_

