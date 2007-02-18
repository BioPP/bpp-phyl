//
// File: DistanceMatrix.h
// Created on: Wed jun 08 10:39 2005
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

#ifndef _DISTANCEMATRIX_H_
#define _DISTANCEMATRIX_H_

// From the STL:
#include <vector>
#include <string>

// From Utils
#include <Utils/Exceptions.h>
#include <NumCalc/VectorExceptions.h> //DimensionException

// From NumCalc
#include <NumCalc/Matrix.h>

/**
 * @brief A Matrix class to store phylogenetic distances.
 */
class DistanceMatrix: public virtual RowMatrix<double> {

	private:
		vector<string> _names;

	public:

    /**
     * @brief Build a new distance matrix with specified names.
     *
     * The dimension of the matrix will be equal to the number of names
     *
     * @param names The names to use.
     */
		DistanceMatrix(const vector<string> & names): RowMatrix<double>(names.size(), names.size()), _names(names)
		{
			reset();
		}

		/**
     * @brief Build a new distance matrix with specified size.
     *
     * Row names will be named 'Taxon 0', 'Taxon 1', and so on.
     *
     * @param n The size of the matrix.
     */
    DistanceMatrix(unsigned int n): RowMatrix<double>(n, n), _names(n)
		{
			for(unsigned int i = 0; i < n; i++) _names[i] = "Taxon " + i;
		}

		virtual ~DistanceMatrix() {}

		DistanceMatrix(const DistanceMatrix & dist): RowMatrix<double>(dist), _names(dist._names)	{}

		DistanceMatrix & operator=(const DistanceMatrix & dist)
		{
			unsigned int n = dist.size();
			resize(n, n);
			for(unsigned int i = 0; i < n; i++)
      {
				for(unsigned int j = 0; j < n; j++)
        {
					operator()(i, j) = dist(i, j);
				}
			}
			_names = dist._names;
			return *this;
		}
		
	public:

    /**i
     * @brief Reset the distance matrix: all distances are set to 0.
     */
		void reset()
		{
			unsigned int n = size();
			for(unsigned int i = 0; i < n; i++) {
				for(unsigned int j = 0; j < n; j++) {
					operator()(i, j) = 0;
				}
			}
		}
		
    /**
     * @return The dimension of the matrix.
     */
		unsigned int size() const { return _names.size(); }

    /**
     * @return The names associated to the matrix.
     */
		vector<string> getNames() const { return _names; }

    /**
     * @return The ith name.
     * @param i Name index.
     * @throw IndexOutOfBoundsException If i is not a valid index.
     */
		string getName(unsigned int i) const throw (IndexOutOfBoundsException)
    { 
      if(i >= size()) throw IndexOutOfBoundsException("DistanceMatrix::getName. Invalid indice.", i, 0, size());
      return _names[i];
    }
    
    /**
     * @brief Set the ith name.
     * 
     * @param i Name index.
     * @param name The new name.
     * @throw IndexOutOfBoundsException If i is not a valid index.
     */
		void setName(unsigned int i, const string & name) throw (IndexOutOfBoundsException)
		{
			if(i >= size()) throw IndexOutOfBoundsException("DistanceMatrix::setName. Invalid indice.", i, 0, size());
			_names[i] = name;
		}

    /**
     * @brief Set the names associated to the matrix.
     * 
     * @param names Matrix names.
     * @throw DimensionException If 'names' have not the same size as the matrix.
     */
		void setNames(const vector<string> & names) throw (DimensionException)
		{
			if(names.size() != _names.size()) throw DimensionException("DistanceMatrix::setNames. Invalid number of names.", names.size(), _names.size());
			_names = names;
		}

    /**
     * @brief Get the index of a given name.
     *
     * @param name The name to look for.
     * @return The position of the name.
     * @throw Exception If no names are attached to this matrix, or if the name was not found.
     */
    unsigned int getNameIndex(const string & name) const throw (Exception);

    /**
     * @brief Access by name.
     *
     * @param iName Name 1 (row)
     * @param jName Name 2 (column)
     * @return A reference toward the specified distance.
     * @throw Exception if the matrix has no name of if one of the name do not match existing names.
     */
    virtual const double & operator()(const string & iName, const string & jName) const throw (Exception)
    {
      unsigned int i = getNameIndex(iName);
      unsigned int j = getNameIndex(jName);
      return operator()(i,j);
    }

    /**
     * @brief Access by name.
     *
     * @param iName Name 1 (row)
     * @param jName Name 2 (column)
     * @return A reference toward the specified distance.
     * @throw Exception if the matrix has no name of if one of the name do not match existing names.
     */
    virtual double & operator()(const string & iName, const string & jName) throw (Exception)
    {
      unsigned int i = getNameIndex(iName);
      unsigned int j = getNameIndex(jName);
      return operator()(i,j);
    }

    virtual const double & operator()(unsigned int i, unsigned int j) const
    {
      return RowMatrix<double>::operator()(i, j);
    }
    virtual double & operator()(unsigned int i, unsigned int j)
    {
      return RowMatrix<double>::operator()(i, j);
    }
    virtual const double & operator()(int i, int j) const
    {
      return RowMatrix<double>::operator()(i, j);
    }
    virtual double & operator()(int i, int j)
    {
      return RowMatrix<double>::operator()(i, j);
    }
};

#endif //_DISTANCEMATRIX_H_

