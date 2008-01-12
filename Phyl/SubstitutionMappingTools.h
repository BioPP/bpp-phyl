//
// File: SubstitutionMappingTools.h
// Created by: Julien Dutheil
// Created on: Wed Apr 5 13:04 2006
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004, 2005, 2006)

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

#ifndef _SUBSTITUTIONMAPPINGTOOLS_H_
#define _SUBSTITUTIONMAPPINGTOOLS_H_

#include "ProbabilisticSubstitutionMapping.h"
#include "SubstitutionCount.h"
#include "DRTreeLikelihood.h"

namespace bpp
{

/**
 * @brief Provide methods to compute substitution mappings.
 *
 * For now, only 4 methods are implemented, and provide probabilistic substitution mappings.
 *
 * See:
 * Dutheil J, Pupko T, Jean-Marie A, Galtier N.
 * A model-based approach for detecting coevolving positions in a molecule.
 * Mol Biol Evol. 2005 Sep;22(9):1919-28. Epub 2005 Jun 8.
 */
class SubstitutionMappingTools
{
  public:
    SubstitutionMappingTools() {}
    virtual ~SubstitutionMappingTools() {}
    
	public:
		
		/**
		 * @brief Compute the substitutions vectors for a particular dataset using the
		 * double-recurvive likelihood computation.
		 *
		 * @param drtl              A DRTreeLikelihood object.
		 * @param substitutionCount The SubstitutionCount to use.
		 * @param verbose           Print info to screen.
		 * @return A vector of substitutions vectors (one for each site).
     * @throw Exception If the likelihood object is not initialized.
		 */
		static ProbabilisticSubstitutionMapping * computeSubstitutionVectors(
			const DRTreeLikelihood & drtl,
			SubstitutionCount & substitutionCount,
			bool verbose = true) throw (Exception);
		
		/**
		 * @brief Compute the substitutions vectors for a particular dataset using the
		 * double-recurvive likelihood computation.
		 *
		 * In this method, substitution counts are computed using the pair of ancestral
		 * states with maximum likelihood.
		 * This is a kind of joint-pair ancestral reconstruction, as in Galtier and Boursot (1998).
		 * This reconstruction possibly takes into account several rate classes, and
		 * substitution counts are averaged over all rate classes, weighted by their conditional
		 * likelihood.
		 *
		 * This function is mainly for testing purpose (see Dutheil et al. 2005).
		 * For practical use, consider using the 'getSubstitutionVectors' method instead.
		 *
		 * @param drtl              A DRTreeLikelihood object.
		 * @param substitutionCount The substitutionsCount to use.
		 * @param verbose           Print info to screen.
		 * @return A vector of substitutions vectors (one for each site).
     * @throw Exception If the likelihood object is not initialized.
		 */
		static ProbabilisticSubstitutionMapping * computeSubstitutionVectorsNoAveraging(
			const DRTreeLikelihood & drtl,
			SubstitutionCount & substitutionCount,
			bool verbose = true) throw (Exception);
		
		/**
		 * @brief Compute the substitutions vectors for a particular dataset using the
		 * double-recurvive likelihood computation.
		 *
		 * In this method, all ancestral states are estimated using marginal likelihoods,
		 * putatively intregated over several rate classes.
		 * For each branch, the number of substitution given marginal states is used.
		 * This method, used with a SimpleSubstitutionCount objet is equivalent to
		 * Tufféry and Darlu's (2000) computation of substitution vectors.
		 *
		 * Use with another substitution count objet is in most cases irrelevent.
		 * 
		 * @param drtl              A DRTreeLikelihood object.
		 * @param substitutionCount The substitutionsCount to use.
		 * @param verbose           Print info to screen.
		 * @return A vector of substitutions vectors (one for each site).
     * @throw Exception If the likelihood object is not initialized.
		 */
		static ProbabilisticSubstitutionMapping * computeSubstitutionVectorsNoAveragingMarginal(
			const DRTreeLikelihood & drtl,
			SubstitutionCount & substitutionCount,
			bool verbose = true) throw (Exception);
		
		/**
		 * @brief Compute the substitutions vectors for a particular dataset using the
		 * double-recurvive likelihood computation.
		 *
		 * The marginal probability is used for weighting, i.e. the product of probabilities for the pair.
		 *
		 * This function is mainly for testing purpose (see Dutheil et al. 2005).
		 * For practical use, consider using the 'getSubstitutionVectors' method instead.
		 *
		 * @param drtl              A DRTreeLikelihood object.
		 * @param substitutionCount The substitutionsCount to use.
		 * @param verbose           Print info to screen.
		 * @return A vector of substitutions vectors (one for each site).
     * @throw Exception If the likelihood object is not initialized.
		 */
		static ProbabilisticSubstitutionMapping * computeSubstitutionVectorsMarginal(
			const DRTreeLikelihood & drtl,
			SubstitutionCount & substitutionCount,
			bool verbose = true) throw (Exception);
	
		/**
		 * @brief Write the substitutions vectors to a stream.
		 *
		 * @param substitutions The substitutions vectors to write.
		 * @param sites         The dataset associated to the vectors
		 * (needed to know the position of each site in the dataset).
		 * @param out           The output stream where to write the vectors.
		 * @throw IOException If an output error happens.
		 */
		static void writeToStream(
			const ProbabilisticSubstitutionMapping & substitutions,
			const SiteContainer & sites,
			ostream & out)
			throw (IOException);
	
		/**
		 * @brief Read the substitutions vectors from a stream.
		 *
		 * @param in            The input stream where to read the vectors.
		 * @param substitutions The mapping object to fill.
		 * @throw IOException If an input error happens.
		 */
		static void readFromStream(istream & in, ProbabilisticSubstitutionMapping & substitutions)
			throw (IOException);
    
};

} //end of namespace bpp.

#endif //_SUBSTITUTIONMAPPINGTOOLS_H_

