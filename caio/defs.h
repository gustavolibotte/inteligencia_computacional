#define USE_MULTI_THREAD

#pragma once

#define isfinite(a) (!(fabs((a)) > DBL_MAX) && !((a) != (a)))

namespace optimization {
	const double EPSILON = 10e-6;
	const double DEN = 1.0;
	double HAMMER_BEST[] = {
		0.437629,
		0.0643525,
		0.27849,
		2.48756,
		4.44427,
		4.46137,
		101.377,
		201.935
	}; 
	double HAMMER_BEST2[] = {
		0.644076,
		0.194823,
		0.583695,
		2.8064,
		4.40238,
		4.51823,
		100,
		201,
	};
	double HAMMER_BEST3[] = {
		0.661133,
		0.214906,
		0.597326,
		2.86149,
		4.96617,
		4.4171,
		100,
		201,
	};

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
	double paper_best[] = {
		0.6280,
		0.1604,
		0.6808,
		2.7087,
		3.0394,
		4.7638,
		101,
		201
	};
}