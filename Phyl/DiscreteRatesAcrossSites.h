#ifndef _DISCRETERATESACROSSSITES_H_
#define _DISCRETERATESACROSSSITES_H_

// From the NumCalc library:
#include <NumCalc/DiscreteDistribution.h>
#include <NumCalc/ParameterList.h>

/**
 * @brief Interface for RAS implementation.
 *
 * This interface provides methods for dealing with RAS models.
 */
class DiscreteRatesAcrossSites
{
	public:
		DiscreteRatesAcrossSites() {}
		~DiscreteRatesAcrossSites() {}

	public:

		/**
		 * @brief Get the rate distribution used for the computation.
		 *
		 * @return A const pointer toward the rate distribution of this instance.
		 */
		virtual const DiscreteDistribution * getRateDistribution() const = 0;

		/**
		 * @brief Get the rate distribution used for the computation.
		 *
		 * @return A pointer toward the rate distribution of this instance.
		 */
		virtual DiscreteDistribution * getRateDistribution() = 0;

		/**
		 * @brief Get the likelihood for a site knowing its rate class.
		 *
		 * @param site      The site index.
		 * @param rateClass The rate class index.
		 * @return The likelihood for the specified site and rate class.
		 */
		virtual double getLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const = 0;
		
		/**
		 * @brief Get the logarithm of the likelihood for a site knowing its rate class.
		 *
		 * @param site      The site index.
		 * @param rateClass The rate class index.
		 * @return The logarithm of the likelihood for the specified site and rate class.
		 */
		virtual double getLogLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const = 0;
	
		/**
		 * @brief Get the likelihood for each site and each rate class.
		 *
		 * @return A two-dimension vector with all likelihoods.
		 */
		virtual VVdouble getLikelihoodForEachSiteForEachRate() const = 0;
		
		/**
		 * @brief Get the logarithm of the likelihood for each site and each rate class.
		 *
		 * @return A two-dimension vector with all log likelihoods:
		 * <code>V[i][j] =</code> likelihood of site i and rate class j.
		 */
		virtual VVdouble getLogLikelihoodForEachSiteForEachRate() const = 0;
		
		/**
		 * @brief Get the posterior probability for each site of belonging to a
		 * particular rate class.
		 *
		 * @return A two-dimension vector with all posterior probabilities:
		 * <code>V[i][j] =</code> probablity for site i of belonging to rate class j.
		 */
		virtual VVdouble getPosteriorProbabilitiesOfEachRate() const = 0;
		
		/**
		 * @brief Get the posterior rate class (the one with maximum posterior
		 * probability) for each site.
		 *
		 * @return A vector with all rate classes indexes.
		 */
		virtual Vint getRateClassWithMaxPostProbOfEachSite() const = 0;

		/**
		 * @brief Get the posterior rate (the one with maximum posterior
		 * probability) for each site.
		 *
		 * @return A vector with all rate classes indexes.
		 */
		virtual Vdouble getRateWithMaxPostProbOfEachSite() const = 0;
	
		/**
		 * @brief Get the posterior rate, i.e. averaged over all classes
		 * and weighted with posterior probabilities, for each site.
		 *
		 * @return A vector with all rates.
		 */
		virtual Vdouble getPosteriorRateOfEachSite() const = 0;

		/**
		 * @brief Get the parameters associated to the rate distirbution.
		 *
		 * @return A ParameterList object with all rate distribution parameters.
		 */
		virtual ParameterList getRateDistributionParameters() const = 0;
		
};

#endif //_DISCRETERATESACROSSSITES_H_
