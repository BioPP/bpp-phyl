#ifndef _RASTOOLS_H_
#define _RASTOOLS_H_

#include "DiscreteRatesAcrossSites.h"

// From NumCalc:
#include <NumCalc/DiscreteDistribution.h>

/**
 * @brief Tools to deal with Rates Across Sites (RAS) models.
 */
class RASTools {

	public:

		/**
		 * @brief Get the rate distribution estimated from a dataset.
		 *
		 * This methods takes an objet implementing the DiscreteRatesAcroossSites
		 * interface as input and use the posterior probabilities of rates for each site
		 * to generate the corresponding distribution.
		 *
		 * @param treeLikelihood A Likelihood calculation implmenting the RAS model interface.
		 * @return The posterior distribution of rate classes.
		 */
		static DiscreteDistribution * getPosteriorRateDistribution(
				const DiscreteRatesAcrossSites & treeLikelihood);
};

#endif // _RASTOOLS_H_

