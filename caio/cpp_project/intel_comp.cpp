// intel_comp.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include <vector>
#include <functional>
#include <algorithm>
#include <Windows.h>
#include <iostream>
#include "individual.h"
#include "bounds.h"
#include "population.h"
#include "functionObjective.h"
#include "algorithm.h"
#include <time.h>
#include "EVAL_HAMMER.h"
#include "Timer.h"
#include "hammer.h"
#include "util.h"
#include "batch.h"

//#define BATCH
#define POINT_EVAL

using namespace optimization;
int _tmain(int argc, _TCHAR* argv[])
{
	Timer t1 = Timer();
	Timer t2 = Timer();

	int p_size = 300;
	int max_it = 100;
	double cross = 1.0;
	double mut = 0.05;
	double cross_f = 1.0;
	double mut_f = 0.05;
	Bounds b = Bounds(8,4);
	b.lower_bounds[0] = 0.508; b.upper_bounds[0] = 1.27;
	b.lower_bounds[1] = 0.0254; b.upper_bounds[1] = 0.254;
	b.lower_bounds[2] = 0.0254; b.upper_bounds[2] = 0.762;

	/*b.lower_bounds[0] = 0.2; b.upper_bounds[0] = 0.5;
	b.lower_bounds[1] = 0.01; b.upper_bounds[1] = 0.1;
	b.lower_bounds[2] = 0.01; b.upper_bounds[2] = 0.3;*/

	b.lower_bounds[3] = 2.0; b.upper_bounds[3] = 5.0;
	b.lower_bounds[4] = 2.0; b.upper_bounds[4] = 5.0;
	b.lower_bounds[5] = 2.0; b.upper_bounds[5] = 5.0;
	b.lower_bounds[6] = 100.0; b.upper_bounds[6] = 101.99;
	b.lower_bounds[7] = 200.0; b.upper_bounds[7] = 202.99;
	

	b.lower_fit[0] = b.lower_fit[1] = b.lower_fit[2] = b.lower_fit[3] = 10e10;
	b.upper_fit[0] = b.upper_fit[1] = b.upper_fit[2] = b.upper_fit[3] = -10e10;

#ifndef BATCH

#ifndef POINT_EVAL
	t1.tic();
	Population p = Population(p_size, &b);
	MOEA<HAMMER_MULT<6>> algo = MOEA<HAMMER_MULT<6>>();

	std::vector<Individual> v1 = algo.reduceSpace(p, cross, mut, max_it);
	t1.toc();
	t2.tic();

	std::vector<Individual> v2 = algo.solveReduced(p, v1, cross, mut, max_it);
	 
	std::cout << "First TPs" << std::endl;
	for(int i = 0; i < v1.size(); i++) {
		v1[i].print();
	}

	std::cout << "\n\nLast TPs" << std::endl;
	for(int i = 0; i < v2.size(); i++) {
		v2[i].print();
	}
	t2.toc();
	t1.print();
	t2.print();
#else
	HammerInterface hm = HammerInterface(L"hammer\\");
	hm.eval2(HAMMER_BEST);
	std::cout << "Kef: " << hm.kef << "\tFluxo: " << hm.fluxo << "\tFP: " << hm.FP << std::endl;
	std::cout << "Fit: " << hm.f << std::endl;
#endif
#else
	BatchRunMOEA<HAMMER_MULT<6>> batch = BatchRunMOEA<HAMMER_MULT<6>>();
	batch.run("batch_result.txt", b, p_size, 100, max_it, cross, cross_f, 4,
					mut, mut_f, 3);
#endif
	system("PAUSE");
	return 0;
}

