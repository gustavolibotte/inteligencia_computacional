#define USE_MULTI_THREAD

#pragma once

#define isfinite(a) (!(fabs((a)) > DBL_MAX) && !((a) != (a)))

namespace optimization {
	/* EPSILON -  Smallest value */
	const double EPSILON = 10e-6;
	/* DEN - Alternative value for denominator equals to zero */
	const double DEN = 1.0;
	/* NC - Crossover distribuition index */
	const double NC = 20;
	/* NM - Mutation distribuition index */
	const double NM = 20;
	double HAMMER_BEST[] = {
		1.0753,
		0.126606,
		0.757567,
		2.28036,
		2.41894,
		4.81758,
		100,
		201
	};//1.2722

	double Kappel[] = {
		1.0628331407,
		0.1437784285,
		0.762,
		2.3794476715,
		2.5730748603,
		5.0,
		100,
		201,
	};
	double Gustavo[] = {
		1.072151,
		0.151212,
		0.758954,
		2.278674,
		2.502342,
		4.688076,
		100,
		201
	};
}