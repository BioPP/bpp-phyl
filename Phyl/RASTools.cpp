#include "RASTools.h"

// From NumCalc:
#include <NumCalc/SimpleDiscreteDistribution.h>

// From the STL

#include <map>

using namespace std;

DiscreteDistribution * RASTools::getPosteriorRateDistribution(
		const DiscreteRatesAcrossSites & treeLikelihood)
{
	
	// Get all posterior rate classes for each sites:
	Vint classes = treeLikelihood.getRateClassWithMaxPostProbOfEachSite();
	map<unsigned int, unsigned int> counts;
	for(unsigned int i = 0; i < classes.size(); i++)
		counts[classes[i]]++;
	
	// Now compute the distribution:
	const DiscreteDistribution * rDist = treeLikelihood.getRateDistribution();
	map<double, double> distribution;
	for(map<unsigned int, unsigned int>::iterator i = counts.begin(); i != counts.end(); i++) {
		distribution[rDist -> getCategory(i -> first)] = (double)i -> second / (double)classes.size();
	}

	// Build a new distribution and return it:
	return new SimpleDiscreteDistribution(distribution);
}

