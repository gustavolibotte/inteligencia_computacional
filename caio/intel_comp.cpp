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
#include "util.h";

#define EVOLVE
#define USE_HAMMER

using namespace optimization;
int _tmain(int argc, _TCHAR* argv[])
{
	Timer t1 = Timer();
	Timer t2 = Timer();
#ifdef EVOLVE
	int p_size = 40;
	int max_it = 2;
	double cross = 1.0;
	double mut = 1.0/p_size;
#ifdef USE_HAMMER
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
	b.lower_bounds[6] = 100.0; b.upper_bounds[6] = 100.99;
	b.lower_bounds[7] = 201.0; b.upper_bounds[7] = 201.99;
	
	b.lower_fit[0] = b.lower_fit[1] = b.lower_fit[2] = b.lower_fit[3] = 10e10;
	t1.tic();
	Population p = Population(p_size, &b);
	MOEA<HAMMER_MULT<4>> algo = MOEA<HAMMER_MULT<4>>();
#else
	Bounds b = Bounds(3,3);
	b.lower_bounds[0] = 0.0; b.upper_bounds[0] = 1.0;
	b.lower_bounds[1] = 0.0; b.upper_bounds[1] = 1.0;
	b.lower_bounds[2] = 0.0; b.upper_bounds[2] = 1.0;
	
	b.lower_fit[0] = b.lower_fit[1] = b.lower_fit[2] = 10e10;
	t2.tic();
	Population p = Population(p_size, &b);
	MOEA<FTZL2> algo = MOEA<FTZL2>();
#endif
	std::vector<Individual> v1 = algo.reduceSpace(p, cross, mut, max_it);

	std::cout << "First TPs" << std::endl;
	for(int i = 0; i < v1.size(); i++) {
		v1[i].print();
	}
	t1.toc();
	t2.tic();

	std::vector<Individual> v2 = algo.solveReduced(p, v1, cross, mut, max_it);
	 

	std::cout << "Last TPs" << std::endl;
	for(int i = 0; i < v2.size(); i++) {
		v2[i].print();
	}
	t2.toc();
	t1.print();
	t2.print();
#else
	t1.tic();
	for(int i=0; i < 40; i++)
		eval2(Kappel);

	t1.toc();
	t2.tic();
	HammerInterface hm = HammerInterface(L"hammer\\");
	for(int i=0; i < 40; i++)
		hm.eval2(Kappel);

	system("cls");
	t2.toc();
	t1.print();
	t2.print();
	std::cout << "Kef: " << hm.kef << "\tFluxo: " << hm.fluxo << "\tFP: " << hm.FP << std::endl;
	std::cout << "Fit: " << hm.f << std::endl;
#endif
	system("PAUSE");
	return 0;
}

