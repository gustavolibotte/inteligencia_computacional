#include <vector>
#include "bounds.h"
#include "individual.h"
#include <algorithm>
#include <limits>
#include "defs.h"

#pragma once
namespace optimization {

	enum PopulationDivision {
		ALL,
		FIRST_HALF,
		SECOND_HALF
	};

	auto lambda = [](int i, int j) -> double {
		return i == j ? 1.0 : EPSILON;
	};
	
	class Population {
		Bounds * bounds_;
		std::vector<Individual> pop;
		std::vector<std::vector<double>> dist;
		int p_size;
	public:

		Population(int pop_size, Bounds * bounds) : bounds_(bounds), p_size(pop_size) {
			pop = std::vector<Individual>(pop_size*2, Individual((*bounds_).lower_bounds.size(), (*bounds_).lower_fit.size()));
			for(int i = 0; i < pop_size; i++) {
				for(int j = 0; j < pop[i].size(); j++) {
					pop[i][j] = (*bounds_).lower_bounds[j] + ((*bounds_).upper_bounds[j] - (*bounds_).lower_bounds[j]) * rand_uniform();
				}
			}
			dist = std::vector<std::vector<double>>(pop_size*2, std::vector<double>(pop_size*2));
		}

		void sort_by_c_asf(PopulationDivision p) {
			switch(p) {
				case PopulationDivision::FIRST_HALF:
					std::sort(pop.begin(), pop.begin() + pop.size() / 2, cmp_ind_by_c_asf);
					break;
				case PopulationDivision::SECOND_HALF:
					std::sort(pop.begin() + pop.size() / 2, pop.end(), cmp_ind_by_c_asf);
					break;
				case PopulationDivision::ALL:
				default:
					std::sort(pop.begin(), pop.end(), cmp_ind_by_c_asf);
					break;
			}
		}

		void sort_by_c_lc(PopulationDivision p) {
			switch(p) {
				case PopulationDivision::FIRST_HALF:
					std::sort(pop.begin(), pop.begin() + pop.size() / 2, cmp_ind_by_c_lc);
					break;
				case PopulationDivision::SECOND_HALF:
					std::sort(pop.begin() + pop.size() / 2, pop.end(), cmp_ind_by_c_lc);
					break;
				case PopulationDivision::ALL:
				default:
					std::sort(pop.begin(), pop.end(), cmp_ind_by_c_lc);
					break;
			}
		}

		void sort_by_c(PopulationDivision p) {
			switch(p) {
				case PopulationDivision::FIRST_HALF:
					std::sort(pop.begin(), pop.begin() + pop.size() / 2, cmp_ind_by_c);
					break;
				case PopulationDivision::SECOND_HALF:
					std::sort(pop.begin() + pop.size() / 2, pop.end(), cmp_ind_by_c);
					break;
				case PopulationDivision::ALL:
				default:
					std::sort(pop.begin(), pop.end(), cmp_ind_by_c);
					break;
			}
		}

		void sort_by_b(PopulationDivision p) {
			switch(p) {
				case PopulationDivision::FIRST_HALF:
					std::sort(pop.begin(), pop.begin() + pop.size() / 2, cmp_ind_by_b);
					break;
				case PopulationDivision::SECOND_HALF:
					std::sort(pop.begin() + pop.size() / 2, pop.end(), cmp_ind_by_b);
					break;
				case PopulationDivision::ALL:
				default:
					std::sort(pop.begin(), pop.end(), cmp_ind_by_b);
					break;
			}
		}

		void sort_by_fit(PopulationDivision p) {
			switch(p) {
				case PopulationDivision::FIRST_HALF:
					std::sort(pop.begin(), pop.begin() + pop.size() / 2, cmp_ind_by_fit);
					break;
				case PopulationDivision::SECOND_HALF:
					std::sort(pop.begin() + pop.size() / 2, pop.end(), cmp_ind_by_fit);
					break;
				case PopulationDivision::ALL:
				default:
					std::sort(pop.begin(), pop.end(), cmp_ind_by_fit);
					break;
			}
		}

		void sort_by_rank(int begin, int end) {
			std::sort(pop.begin() + begin, pop.begin() + end, cmp_ind_by_r);
		}
		
		void sort_by_dist_asc(int begin, int end) {
			std::sort(pop.begin() + begin, pop.begin() + end, cmp_ind_by_dist_asc);
		}

		void sort_by_dist_desc(int begin, int end) {
			std::sort(pop.begin() + begin, pop.begin() + end, cmp_ind_by_dist_desc);
		}

		void update_z_min(PopulationDivision pd) {
			int begin = 0, end = 0;
			getIndex(pd, begin, end);
			assert(begin != end);

			for(int i = 0; i < pop[begin].fit.size(); i++) {
				(*bounds_).lower_fit[i] = pop[begin].fit[i];
			}

			for(int i = begin + 1; i < end; i++) {
				for(int j = 0; j < pop[i].fit.size(); j++) {
					(*bounds_).lower_fit[j] = pop[i].fit[j] < (*bounds_).lower_fit[j] ? pop[i].fit[j] : (*bounds_).lower_fit[j];
				}
			}
		}

		void update_z_max(int begin, int end) {
			for(int i = 0; i < pop[begin].fit.size(); i++) {
				(*bounds_).upper_fit[i] = pop[begin].fit[i];
			}

			for(int i = begin + 1; i < end; i++) {
				for(int j = 0; j < pop[i].fit.size(); j++) {
					(*bounds_).upper_fit[j] = pop[i].fit[j] > (*bounds_).upper_fit[j] ? pop[i].fit[j] : (*bounds_).upper_fit[j];
				}
			}
		}

		void ensureBounds(PopulationDivision p) {
			switch(p) {
				case PopulationDivision::FIRST_HALF:
					(*bounds_).ensureBounds(pop.begin(), pop.begin() + pop.size() / 2);
					break;
				case PopulationDivision::SECOND_HALF:
					(*bounds_).ensureBounds(pop.begin() + pop.size() / 2, pop.end());
					break;
				case PopulationDivision::ALL:
				default:
					(*bounds_).ensureBounds(pop.begin(), pop.end());
					break;
			}
		}

		void getIndex(PopulationDivision pd, int & begin, int & end) {
			switch(pd) {
				case PopulationDivision::FIRST_HALF:
					begin = 0;
					end = pop.size() / 2;
					break;
				case PopulationDivision::SECOND_HALF:
					begin = pop.size() / 2;
					end = pop.size();
					break;
				case PopulationDivision::ALL:
				default:
					begin = 0;
					end = pop.size();
					break;
			}
		}

		void update_asf(PopulationDivision pd) {
			int begin = 0, end = 0;
			getIndex(pd, begin, end);
			assert(begin != end);
			for(int i = begin; i < end; i++) {
				double max = -(std::numeric_limits<double>::max)();
				for(int j = 0; j < (*bounds_).lower_fit.size(); j++) {
					double tmp = (pop[i].fit[j] - (*bounds_).lower_fit[j]) /  lambda(j, pop[i].lc_class);
					max = tmp > max ? tmp : max;
				}
				assert(isfinite(max));
				pop[i].asf = max;
			}
		}

		void updateBounded(PopulationDivision pd) {
			int begin = 0, end = 0;
			getIndex(pd, begin, end);
			assert(begin != end);
			for(int i = begin; i < end; i++) {
				pop[i].is_bounded = true;
				for(int j = 0; j < pop[i].fit.size(); j++) {
					if(pop[i].fit[j] > (*bounds_).upper_fit[j]) {
						pop[1].is_bounded = false;
						break;
					}
				}
			}
		}

		void restart(std::vector<Individual> tps) {

			std::vector<double> range((*bounds_).lower_bounds.size());
			for(int i = 0; i < (*bounds_).lower_bounds.size(); i++) {
				range[i] = ((*bounds_).upper_bounds[i] - (*bounds_).lower_bounds[i]) / 2.0;
			}
			int i = 0;
			while(i < p_size) {
				for(int k = 0; k < tps.size(); k++) {
					for(int j = 0; j < pop[i].size(); j++) {
						pop[i][j] = tps[k][j] + range[j] * (rand_uniform() - 0.5);
					}
					pop[i].needUpdate.fit = true;
					i++;
					if(i >= p_size) {
						break;
					}
				}
			}
			ensureBounds(PopulationDivision::FIRST_HALF);
		}

		int s_in(PopulationDivision pd) {
			int begin = 0, end = 0;
			getIndex(pd, begin, end);
			assert(begin != end);
			int ret = 0;
			for(int i = begin; i < end; i++) {
				if(pop[i].is_bounded) {
					ret++;
				}
			}
			return ret;
		}

		int nd(int begin, int end) {
			int ret = 0;
			for(int i = begin; i < end; i++) {
				if(pop[i].rank == 0) {
					ret++;
				}
			}
			return ret;
		}

		void rank_by_dom(int begin, int end) {
			for(int i = begin; i < end; i++) {
				pop[i].rank = 0;
				for(int j = begin; j < end; j++) {
					if(i != j) {
						if(dominate(pop[j], pop[i]) && !dominate(pop[i], pop[j])) {
							pop[i].rank++;
						}
					}
				}
			}
			sort_by_rank(begin, end);
		}

		void update_grid(int begin, int end) {
			for(int i = begin; i < end; i++) {
				for(int j = 0; j < pop[i].grid.size(); j++) {
					double den = ((*bounds_).upper_fit[j] - (*bounds_).lower_fit[j]);
					den = den == 0.0 ? 10e-6 : den;
					pop[i].grid[j] = p_size * ((pop[i].fit[j] - (*bounds_).lower_fit[j]) / den);
				}
			}
		}

		void update_dist_grid(int begin, int end) {
			for(int i = begin; i < end; i++) {
				for(int j = i+1; j < end; j++) {
					dist[i][j] = 0;
					for(int k = 0; k < pop[i].grid.size(); k++) {
						dist[i][j] += (pop[i].grid[k] - pop[j].grid[k]) * (pop[i].grid[k] - pop[j].grid[k]);
					}
					dist[i][j] = sqrt(dist[i][j]);
					assert(isfinite(dist[i][j]));
					dist[j][i] = dist[i][j];
				}
			}
		}

		void update_min_dist(int begin, std::vector<Individual> tps) {
			for(int k = begin; k < pop.size(); k++) {
				pop[k].dist = (std::numeric_limits<double>::max)();
				for(int i = 0; i < tps.size(); i++) {
					double dist = 0;
					for(int l = 0; l < tps[i].size(); l++) {// for tp[i]
						dist += (tps[i][l] -  pop[i][l]) * (tps[i][l] -  pop[i][l]);
					}
					pop[k].dist = dist < pop[k].dist ? dist : pop[k].dist;
					for(int j = i+1; j < tps.size(); j++) {
						dist = 0;
						for(int l = 0; l < tps[i].size(); l++) {// for middle tp[i] tp[j]
							dist += ((tps[i][l] + tps[j][l]) / 2.0  -  pop[i][l]) * ((tps[i][l] + tps[j][l]) / 2.0 -  pop[i][l]);
						}
						pop[k].dist = dist < pop[k].dist ? dist : pop[k].dist;
					}
				}
			}
		}

		void update_max_dist(int from_begin, int from_end, int to_begin, int to_end) {
			for(int i = from_begin; i < from_end; i++) {
				pop[i].dist = -(std::numeric_limits<double>::max)();
				for(int j = to_begin; j < to_end; j++) {
					double dist = 0;
					for(int k = 0; k < pop[i].size(); k++) {
						dist += (pop[i][k] - pop[j][k]) * (pop[i][k] - pop[j][k]);
					}
					pop[i].dist = pop[i].dist < dist ? dist : pop[i].dist;
				}
			}
		}

		void sort_by_dist_grid(int begin, int end) {
			std::vector<bool> mark = std::vector<bool>(end - begin, false);
			for(int i = 0; i < (end - begin - p_size) ; i++) {
				int idx = -1;
				pop[i].rank = 0;
				double min = (std::numeric_limits<double>::max)();
				for(int j = begin; j < end; j++) {
					if(mark[j-begin]) {
						continue;
					}
					for(int k = begin; k < end; k++) {
						if(i==j || mark[k-begin]) {
							continue;
						}
						if(dist[j][k] < min) {
							idx = i;
							min = dist[j][k];
						}
					}
				}
				mark[idx] = true;
				pop[idx].rank = 1;
			}
			sort_by_rank(begin, end);
		}

		Individual & operator [](int i) {
			return pop[i];
		}

		std::vector<Individual> & operator *() {
			return pop;
		}

		int realSize() {
			return pop.size();
		}

		int initialSize() {
			return p_size;
		}

		Bounds * getBounds() {
			return bounds_;
		}
	};
}