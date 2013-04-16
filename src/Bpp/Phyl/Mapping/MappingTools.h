//
// File: SubstitutionMappingTools.h
// Created by: Julien Dutheil
// Created on: Wed Apr 5 13:04 2006
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

#ifndef _MAPPINGTOOLS_H_
#define _MAPPINGTOOLS_H_

#include "ProbabilisticSubstitutionMapping.h"
#include "SubstitutionCount.h"
#include "../Likelihood/DRTreeLikelihood.h"

namespace bpp
{

  /**
   * @brief Provide general methods for count and estimation using substitution mappings.
   *
   * @author Julien Dutheil, Laurent Guéguen
   */
  
  class MappingTools
  {
  public:
    MappingTools() {}
    virtual ~MappingTools() {}
    
  public:
		
    /*
     * @brief Returns the counts on each branch.
     *
     * @param drtl              A DRTreeLikelihood object.
     * @param ids               The numbers of the nodes of the tree
     * @param model             The model on which the SubstitutionCount is built
     * @param reg               the Substitution Register
     * @param threshold         value above which counts are considered saturated
     *                                        (default: -1 means no threshold).
     * @return A vector of substitutions vectors (one per branch per type).
     */

    static std::vector< std::vector<double> > getCountsPerBranch(
                                                                 DRTreeLikelihood& drtl,
                                                                 const std::vector<int>& ids,
                                                                 SubstitutionModel* model,
                                                                 const SubstitutionRegister& reg,
                                                                 double threshold = -1);

    /*
     * @brief Returns the normalization factors due to the null model
     * on each branch, for each register
     *
     *
     * @param drtl              A DRTreeLikelihood object.
     * @param ids               The numbers of the nodes of the tree
     * @param nullModel         The model on which the SubstitutionCount is built
     * @param reg               the Substitution Register
     * @return A vector of normalization vectors (one per branch per type).
     */
    

    static std::vector< std::vector<double> > getNormalizationsPerBranch(
                                                                         DRTreeLikelihood& drtl,
                                                                         const std::vector<int>& ids,
                                                                         const SubstitutionModel* nullModel,
                                                                         const SubstitutionRegister& reg);

    /*
     * @brief Returns the normalization factors due to the set of null
     * models on each branch, for each register.
     *
     *
     * @param drtl              A DRTreeLikelihood object.
     * @param ids               The numbers of the nodes of the tree
     * @param nullModelSet      The model on which the SubstitutionCount is built
     * @param reg               the Substitution Register
     * @return A vector of normalization vectors (one per branch per type).
     */

    static std::vector< std::vector<double> > getNormalizationsPerBranch(
                                                                         DRTreeLikelihood& drtl,
                                                                         const std::vector<int>& ids,
                                                                         const SubstitutionModelSet* nullModelSet,
                                                                         const SubstitutionRegister& reg);

    
    /* *
     * @brief Returns the counts normalized by a null model
     *
     * @param drtl              A DRTreeLikelihood object.
     * @param ids               The numbers of the nodes of the tree
     * @param model             The model on which the SubstitutionCount is built
     * @param nullModel         The model on which the SubstitutionCount is built
     * @param reg               the Substitution Register
     * @param threshold         value above which counts are considered saturated
     *                                        (default: -1 means no threshold).
     */

    static std::vector< std::vector<double> >  getNormalizedCountsPerBranch(
                                                                            DRTreeLikelihood& drtl,
                                                                            const std::vector<int>& ids,
                                                                            SubstitutionModel* model,
                                                                            SubstitutionModel* nullModel,
                                                                            const SubstitutionRegister& reg,
                                                                            double threshold = -1);
    /* *
     * @brief Returns the counts normalized by a null model set
     *
     * @param drtl              A DRTreeLikelihood object.
     * @param ids               The numbers of the nodes of the tree
     * @param model             The model on which the SubstitutionCount is built
     * @param nullModelSet      The model on which the SubstitutionCount is built
     * @param reg               the Substitution Register
     * @param threshold         value above which counts are considered saturated
     *                                        (default: -1 means no threshold).
     */


    static std::vector< std::vector<double> >  getNormalizedCountsPerBranch(
                                                                            DRTreeLikelihood& drtl,
                                                                            const std::vector<int>& ids,
                                                                            SubstitutionModelSet* modelSet,
                                                                            SubstitutionModelSet* nullModelSet,
                                                                            const SubstitutionRegister& reg,
                                                                            double threshold = -1);

    /*
     * @brief Returns the counts relative to the frequency of the
     * states in case of non-stationarity.
     *
     * @param drtl              A DRTreeLikelihood object.
     * @param ids               The numbers of the nodes of the tree
     * @param model             The model on which the SubstitutionCount is built
     * @param reg               the Substitution Register
     * @param stationarity      if false, a correction is made if the SubstitutionRegister
     *                             is a CategorySubstitutionRegister
     * @param threshold         value above which counts are considered saturated
     *                                        (default: -1 means no threshold).
     *
     */

    static std::vector< std::vector<double> >  getRelativeCountsPerBranch(
                                                                          DRTreeLikelihood& drtl,
                                                                          const std::vector<int>& ids,
                                                                          SubstitutionModel* model,
                                                                          const SubstitutionRegister& reg,
                                                                          bool stationarity = true,
                                                                          double threshold = -1);

    /*
     * @brief Output the sum of the counts par branch per site, in a
     * file.
     *
     * @param filename          The name of the output file
     * @param drtl              A DRTreeLikelihood object.
     * @param ids               The numbers of the nodes of the tree
     * @param model             The model on which the SubstitutionCount is built
     * @param reg               the Substitution Register
     *
     */

    static void outputTotalCountsPerBranchPerSite(
                                                  std::string& filename,
                                                  DRTreeLikelihood& drtl,
                                                  const std::vector<int>& ids,
                                                  SubstitutionModel* model,
                                                  const SubstitutionRegister& reg);

  };
} //end of namespace bpp.

#endif //_MAPPINGTOOLS_H_

