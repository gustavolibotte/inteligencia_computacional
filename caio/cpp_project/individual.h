#include <vector>
#include "defs.h"
#include <assert.h>
#include <random>
#include <iostream>

#pragma once
namespace optimization {

	std::random_device rd;
	std::knuth_b e2(rd());
	std::knuth_b e1(rd());
	std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
	/*
	minstd_rand0	
	minstd_rand	
	mt19937	
	mt19937_64	
	ranlux24_base	
	ranlux48_base
	ranlux24
	ranlux48
	*knuth_b
	*/

	class Individual {
		struct Update {
			bool fit;
		};
	public:
		std::vector<double> x;
		std::vector<double> fit;
		std::vector<double> grid;
		Update needUpdate;
		
		double asf;

		int lc_class;
		double lambda_c;
		bool is_bounded;
		int rank;
		double dist;

		Individual(int i_size, int f_size) {
			x = std::vector<double>(i_size);
			fit = std::vector<double>(f_size);
			grid = std::vector<double>(f_size);
			needUpdate.fit = true;
			is_bounded = false;
		}

		int size() const {
			return x.size();
		}

		double const & operator [](int i) const {
			return x[i];
		}

		double & operator [](int i) {
			return x[i];
		}

		void print() {
			std::cout << "x:\t";
			for(int i = 0; i < x.size(); i++) {
				std::cout << x[i] << ", ";
			}
			std::cout << std::endl;

			std::cout << "Fit:\t";
			for(int i = 0; i < fit.size(); i++) {
				std::cout << fit[i] << ", ";
			}
			std::cout << std::endl;
			std::cout << "ASF:\t" << asf << std::endl;
		}
	};

	auto dominate = [](const Individual &i1, const Individual &i2) -> bool {
						for(int i = 0; i < i1.fit.size(); i++) {
							if(i1.fit[i] > i2.fit[i]) {
								return false;
							}
						}
						return true;
					};
	auto cmp_ind_by_c_asf = [](const Individual &i1, const Individual &i2) -> bool { 
						if(i1.lc_class == i2.lc_class) {
							return i1.asf < i2.asf;
						} else {
							return i1.lc_class < i2.lc_class;
						}
					};
	auto cmp_ind_by_c_lc = [](const Individual &i1, const Individual &i2) -> bool { 
						if(i1.lc_class == i2.lc_class) {
							return i1.lambda_c < i2.lambda_c;
						} else {
							return i1.lc_class < i2.lc_class;
						}
					};
	auto cmp_ind_by_c = [](const Individual &i1, const Individual &i2) -> bool { 
						return i1.lc_class < i2.lc_class;
					};

	auto cmp_ind_by_b = [](const Individual &i1, const Individual &i2) -> bool { 
						return (i1.is_bounded && !i2.is_bounded) ? true : false;
					};
	auto cmp_ind_by_r = [](const Individual &i1, const Individual &i2) -> bool { 
						return i1.rank < i2.rank;
					};
	auto cmp_ind_by_dist_asc = [](const Individual &i1, const Individual &i2) -> bool { 
						return i1.dist < i2.dist;
					};
	auto cmp_ind_by_dist_desc = [](const Individual &i1, const Individual &i2) -> bool { 
						return i1.dist > i2.dist;
					};
	auto cmp_ind_by_fit = [](const Individual &i1, const Individual &i2) -> bool { 
						double sum1 = 0, sum2 = 0;
						for(int i = 0; i < i1.fit.size(); i++) {
							sum1 += i1.fit[i];
							sum2 += i2.fit[i];
						}
						return sum1 < sum2;
					};

	double rand_uniform() {
		//return (double)rand()/(double)RAND_MAX;
		
		 return uniform_dist(e2);
	}

	double cross_fact() {
		double r = rand_uniform();
		double b = pow((2.0 * r), (1.0 / (NC + 1.0)));

		double b2 = b > 1.0 ? (r == 1.0 ? 0 : 1 / pow(2 * (1-r), 1.0 / (NC + 1.0))) : b;
		assert(isfinite(b2));
		return b2;
	}

	double mut_fact() {
		double r = rand_uniform();
		if(r > 0.5) {
			r = pow((2.0 * r), (1.0 / (NM + 1.0)) - 1.0);
		} else {
			r = pow((2*(1-r)), (1.0 / (NM + 1.0)));
		}
		assert(isfinite(r));
		return r;
	}

	void bin_crossover_p(Individual & p1, Individual & p2, Individual & c, double rate) {
		for(int i = 0; i < c.size(); i++) {
			c[i] = 0.5 * ((1 - rate) * p1[i] + (1 + rate) * p2[i]);
		}
	}

	void bin_crossover_m(Individual & p1, Individual & p2, Individual & c, double rate) {
		for(int i = 0; i < c.size(); i++) {
			c[i] = 0.5 * ((1 + rate) * p1[i] + (1 - rate) * p2[i]);
		}
	}

	void bin_mutation(Individual & p1, Individual & p2, Individual & c, double rate) {
		for(int i = 0; i < c.size(); i++) {
			double range = abs(p1[i] - p2[i]);
			c[i] = c[i] + range * rate;
		}
	}
}