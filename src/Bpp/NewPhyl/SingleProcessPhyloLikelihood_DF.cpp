//
// File: SingleProcessPhyloLikelihood_DF.cpp
// Authors: François Gindraud, Laurent Guéguen
// Creation: lundi 27 mai 2019, à 06h 35
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.

  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".

  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.

  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#include "SingleProcessPhyloLikelihood_DF.h"

using namespace bpp;
using namespace dataflow;
Vdouble SingleProcessPhyloLikelihood_DF::getLikelihoodPerSite() const
{
  auto vLik=getLikelihoodCalculation()->getSiteLikelihoods()->getTargetValue();
  Vdouble v(vLik.size());
  
  Eigen::VectorXd::Map(&v[0], v.size()) = vLik;
  return v;
}
      

VVdouble SingleProcessPhyloLikelihood_DF::getPosteriorProbabilitiesPerClass() const
{
  auto rates=getLikelihoodCalculation()->getSubstitutionProcess().getRateDistribution();
  auto nbS=getLikelihoodCalculation()->getNumberOfSites();
  VVdouble vv(nbS);
  
  if (!rates || rates->getNumberOfCategories()==1)
  {
    for (auto& v:vv)
      v.resize(1,1);
  }
  else
  {
    auto vvLik=getLikelihoodCalculation()->getSiteLikelihoodsForAllClasses();
    for (size_t i=0;i<nbS;i++)
    {
      vv[i].resize(vvLik.rows());
      Eigen::VectorXd::Map(&vv[i][0], vv[i].size()) = vvLik.col(i)/vvLik.col(i).sum();
    }
  }
  return vv;
}
      
Vdouble SingleProcessPhyloLikelihood_DF::getPosteriorProbabilitiesForSitePerClass(size_t pos) const
{
  auto rates=getLikelihoodCalculation()->getSubstitutionProcess().getRateDistribution();
  
  if (!rates || rates->getNumberOfCategories()==1)
    return Vdouble(1,1);
  else
  {
    Vdouble vv(rates->getNumberOfCategories());
    for (size_t i=0;i<vv.size();i++)
      vv[i]=(getLikelihoodCalculation()->getSiteLikelihoodsForAClass(i))(pos);
      
    vv/=VectorTools::sum(vv);
    return vv;
  }
}
