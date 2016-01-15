//
// File: SENCA.h
// Created by: Fanny Pouyet 
// Created on: mars 2012
//
 
/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

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
 
#ifndef _SENCA_H_
#define _SENCA_H_
 
#include "AbstractCodonFitnessSubstitutionModel.h"
#include "AbstractCodonSubstitutionModel.h"
#include "AbstractCodonDistanceSubstitutionModel.h"

namespace bpp
{
  /**
   * @brief Class for non-synonymous and synonymous substitution models
   * on codons with parameterized equilibrium frequencies and
   * nucleotidic basic models.
   *
   * @author Fanny Pouyet, Laurent Guéguen
   *
   * See description in AbstractCodonDistanceSubstitutionModel
   * class, AbstractCodonFitnessSubstitutionModel class.
   *
   * Only substitutions with one letter changed are accepted. </p>
   *
   * The ratio between non-synonymous and synonymous substitutions 
   * rates is @f$\beta@f$ with positive parameter \c "beta".
   *
   * The fitness of a codon is a value between 0 and 1 defining the
   * relative advantage of a codon, compared to others. If a codon
   * @f$i@f$ has a fitness @f$\phi_i@f$ and another one (@f$j@f$)
   * has a fitness @f$\phi_j@f$, the substitution rate from codon
   * @f$i@f$ to codon @f$j@f$ is multiplied by
   * \f[-\frac{ \log(\frac{\phi_i}{\phi_j})}{1-\frac{\phi_i}{\phi_j}}\f]
   *
   * The set of fitnesses is implemented through a Codon
   * FrequenciesSet object. The parameters are named \c
   * "fit_NameOfTheParameterInTheFrequenciesSet".
   *
   * Reference:
   * -  Yang Z. and Nielsen R. (2008), _Molecular Biology and Evolution_ 25(3):568--579.
   */
  
  class SENCA :
    public virtual SubstitutionModel,
    public AbstractCodonSubstitutionModel,
    public AbstractCodonDistanceSubstitutionModel,
    public AbstractCodonFitnessSubstitutionModel
  {
  public:
    SENCA(
        const GeneticCode* gCode,
        NucleotideSubstitutionModel* pmod,
        FrequenciesSet* pfit,
        const AlphabetIndex2* pdist = 0);
    SENCA(
        const GeneticCode* gCode,
        NucleotideSubstitutionModel* pmod1,
        NucleotideSubstitutionModel* pmod2,
        NucleotideSubstitutionModel* pmod3,
        FrequenciesSet* pfit,
        const AlphabetIndex2* pdist = 0);

    virtual ~SENCA() {}

    SENCA* clone() const
    {
      return new SENCA(*this);
    }

  public:
    void fireParameterChanged(const ParameterList& parameterlist);

    std::string getName() const;

    double getCodonsMulRate(size_t i, size_t j) const;

    void setNamespace(const std::string&);

    /*
     * @brief set the fitness of the model from
     * given frequencies, such that the equilibrium frequencies of the
     * model matches at best the given ones.
     * 
     * Matching is done in two steps : first, frequencies from nucleotide substitution model are
     * matched at best, then the resulting discrepancy (in terms of
     * ratios between the given one and the one computed by the pahse
     * frequencies) is given for matching to the fitness.
     *
     * @ param frequencies  the frequencies to match on.
     */
    void setFreq(std::map<int,double>& frequencies);

  };

} // end of namespace bpp.

#endif

