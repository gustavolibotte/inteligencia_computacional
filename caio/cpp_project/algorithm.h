#include "population.h"
#include "individual.h"
#include "bounds.h"
#include <vector>
#include <list>
#include <limits>
#include "defs.h"
#include <random>

#pragma once

namespace optimization {

	template<class fObj>
	class MOEA {
		fObj myfObj;
		void updateFit(Population & p, PopulationDivision pd);
		void classify(Population & p, PopulationDivision pd);
		void balance(Population & p, PopulationDivision pd);
		std::vector<int> findPartitions(Population & p, PopulationDivision pd);
		void genOffsprings(Population & p, std::vector<int> partitions, double c_rate, double m_rate);
		void reorderHalving(Population & p, std::vector<int> partition);
		std::vector<int> findTargetPoints(Population & p, PopulationDivision pd);

		//////////////////////////////////////////////////////

	public:

		MOEA() {
			//myfObj = fObj();
		}

		void run(Population & p);
		std::vector<Individual> reduceSpace(Population & p, double cross_rate, double mut_rate, int max_it);
		std::vector<Individual> solveReduced(Population & p, std::vector<Individual> tps, double cross_rate, double mut_rate, int max_it);
	};

	template<class fObj>
	void MOEA<fObj>::updateFit(Population & p, PopulationDivision pd) {
		int begin = 0, end = 0;
		p.getIndex(pd, begin, end);
		assert(begin != end);
#ifdef USE_MULTI_THREAD
		myfObj.eval(*p, begin, end);
#else
		for(int i = begin; i < end; i++) {
			if(p[i].needUpdate.fit) {
				myfObj.eval(p[i].x, p[i].fit);
				p[i].needUpdate.fit = false;
			}
		}
#endif
	}

	template<class fObj>
	void MOEA<fObj>::classify(Population & p, PopulationDivision pd) {

		int begin = 0, end = 0;
		p.getIndex(pd, begin, end);
		assert(begin != end);

		Bounds * b = p.getBounds();

		for(int i = begin; i < end; i++) {

			int iClass = -1;
			double min = (std::numeric_limits<double>::max)();

			for(int j = 0; j < (*b).lower_fit.size(); j++) {

				double max = -(std::numeric_limits<double>::max)();
				double den = ((*b).lower_fit[j] == p[i].fit[j]) ? DEN : p[i].fit[j] - (*b).lower_fit[j];

				for(int k = 0; k < (*b).lower_fit.size(); k++) {
					if(j == k) {
						max = EPSILON;
					} else {
						double lc_ik = (p[i].fit[k] - (*b).lower_fit[k]) / den;
						max = max > lc_ik ? max : lc_ik;
					}
				}
				if(max < min) {
					min = max;
					iClass = j;
				}
			}
			p[i].lc_class = iClass;
			p[i].lambda_c = min;
		}
	}

	template<class fObj>
	void MOEA<fObj>::balance(Population & p, PopulationDivision pd) {
		int num_classes = p[0].fit.size();
		int begin = 0, end = 0;

		p.getIndex(pd, begin, end);
		assert(begin != end);

		int min_ind_per_class = floor((end - begin) / num_classes);

		std::vector<int> class_counter(num_classes);
		std::vector<int> worst_idx(num_classes);

		p.sort_by_c_lc(pd);

		for(int i = begin; i < end; i++) {
			class_counter[p[i].lc_class]++;
			worst_idx[p[i].lc_class] = i;
		}

		//transfer 1 by 1 from max to min
		while(true) {
			int min = findMin(class_counter);
			if(class_counter[min] >= min_ind_per_class) {
				break;
			}
			int max = findMax(class_counter);
			int i = worst_idx[max];
			p[i].lc_class = min;
			worst_idx[max]--;
			class_counter[max]--;
			class_counter[min]++;
		}
		p.sort_by_c(pd);
	}
	template<class fObj>
	std::vector<int> MOEA<fObj>::findPartitions(Population & p, PopulationDivision pd) {
		assert(p.realSize() > 0);
		assert(p[0].fit.size() > 0);

		int begin = 0, end = 0;
		p.getIndex(pd, begin, end);
		assert(begin != end);
		int cur = p[begin].lc_class;
		int idx = 0;
		std::vector<int> part(p[0].fit.size(), end);
		for(int i = begin + 1; i < end; i++) {
			if(p[i].lc_class > cur) {
				part[idx] = i;
				idx++;
				cur = p[i].lc_class;
			}
		}
		return part;
	}
	template<class fObj>
	void MOEA<fObj>::genOffsprings(Population & p, std::vector<int> partitions, double c_rate, double m_rate) {
		int begin = 0, end = 0, mid = 0;
		
		end = p.realSize();
		mid = end / 2;
		bool up = true;
		double cross, mut;
		int i1, i2;
		for(int i = 0; i < partitions.size(); i++) {
			std::uniform_int_distribution<int> ui(begin, partitions[i]);
			for(int j = begin; j < partitions[i]; j++) {
				if(rand_uniform() <= c_rate) {
					if(up) {
						cross = cross_fact();
						i1 = ui(e1);
						i2 = ui(e1);
						up = false;
						bin_crossover_p(p[i1], p[i2], p[mid], cross);
					} else {
						up = true;
						bin_crossover_m(p[i1], p[i2], p[mid], cross);
					}
					if(rand_uniform() <= m_rate) {
						mut = mut_fact();
						bin_mutation(p[i1], p[i2], p[mid], mut);
					}
				} else {
					p.restartIndividual(mid);
				}
				p[mid].needUpdate.fit = true;
				p[mid].lc_class = i;
				mid++;
			}
			begin = partitions[i];
		}
	}

	template<class fObj>
	void MOEA<fObj>::reorderHalving(Population & p, std::vector<int> partitions) {
		int last = 0;
		int c = partitions.size() + 1;
		for(int i = 0; i < partitions.size(); i++) {
			for(int j = last + ((partitions[i] - last) / 2); j < partitions[i]; j++) {
				p[j].lc_class = c;
			}
			last = partitions[i];
		}
		p.sort_by_c(PopulationDivision::ALL);
	}

	template<class fObj>
	std::vector<int> MOEA<fObj>::findTargetPoints(Population & p, PopulationDivision pd) {
		int begin = 0, end = 0;
		p.getIndex(pd, begin, end);
		p.update_z_min(PopulationDivision::FIRST_HALF);
		std::vector<int> v = std::vector<int>(p[0].fit.size());
		Bounds * b = p.getBounds();
		for(int i = 0; i < v.size(); i++) {//for each tp
			double min = (std::numeric_limits<double>::max)();
			int idx = -1;
			for(int j = begin; j < end; j++) {//for each individual
				double max = (p[j].fit[0] - (*b).lower_fit[0]) /  lambda(i, 0);
				for(int k = 1; k < v.size(); k++) {//for each fit
					double tmp = (p[j].fit[k] - (*b).lower_fit[k]) /  lambda(i, k);
					max = tmp > max ? tmp : max;
				}
				if(max < min) {
					idx = j;
					min = max;
				}
			}
			assert(idx != -1);
			v[i] = idx;
		}
		return v;
	}

	template<class fObj>
	std::vector<Individual> MOEA<fObj>::reduceSpace(Population & p, double cross_rate, double mut_rate, int max_it) {
			updateFit(p, PopulationDivision::FIRST_HALF);
			p.update_z_min(PopulationDivision::FIRST_HALF);
			classify(p, PopulationDivision::FIRST_HALF);
			balance(p, PopulationDivision::FIRST_HALF);
			for(int i = 0; i < max_it; i++) {
				if((i*100)% max_it == 0) {
					std::cout << (i*100)/ max_it << " % Complete" << std::endl;
				}
				std::vector<int> partitions = findPartitions(p, PopulationDivision::FIRST_HALF);
				genOffsprings(p, partitions, cross_rate, mut_rate);
				p.ensureBounds(PopulationDivision::SECOND_HALF);
				updateFit(p, PopulationDivision::SECOND_HALF);
				p.update_z_min(PopulationDivision::ALL);
				p.update_asf(PopulationDivision::ALL);
				p.sort_by_c_asf(PopulationDivision::ALL);
				partitions = findPartitions(p, PopulationDivision::ALL);
				reorderHalving(p, partitions);
				p.update_z_min(PopulationDivision::FIRST_HALF);
				classify(p, PopulationDivision::FIRST_HALF);
				balance(p, PopulationDivision::FIRST_HALF);
			}
			
			std::vector<int> tps = findTargetPoints(p, PopulationDivision::FIRST_HALF);
			std::vector<Individual> ret = std::vector<Individual>();
			for(int i = 0; i < tps.size(); i++) {
				ret.push_back(p[tps[i]]);
			}
			return ret;
		}
	/////////////////////////////////////

	template<class fObj>
	std::vector<Individual> MOEA<fObj>::solveReduced(Population & p, std::vector<Individual> tps, double cross_rate, double mut_rate, int max_it) {
		p.restart(tps);
		updateFit(p, PopulationDivision::FIRST_HALF);
		for(int i = 0; i < max_it; i++) {
			if((i*100)% max_it == 0) {
				std::cout << (i*100)/ max_it << " % Complete" << std::endl;
			}
			p.update_z_max(0, p.initialSize());
			genOffsprings(p, std::vector<int>(1, p.initialSize()), cross_rate, mut_rate);
			p.ensureBounds(PopulationDivision::SECOND_HALF);
			updateFit(p, PopulationDivision::SECOND_HALF);
			p.updateBounded(PopulationDivision::ALL);
			p.sort_by_b(PopulationDivision::ALL);
			int s_in = p.s_in(PopulationDivision::ALL);
			if(s_in > p.initialSize()) {
				p.rank_by_dom(0, s_in);
				int nd = p.nd(0, s_in);
				if(nd > p.initialSize()) {
					p.update_z_min(PopulationDivision::ALL);
					p.update_z_max(0, nd);
					p.update_grid(0, nd);
					p.update_dist_grid(0, nd);
					p.sort_by_dist_grid(0, nd);
				} else {
					p.update_max_dist(nd, s_in, 0, nd);
					p.sort_by_dist_asc(nd, s_in);
				}
			} else {
				p.update_min_dist(s_in, tps);
				p.sort_by_dist_desc(s_in, p.realSize());
			}
			p.update_z_min(PopulationDivision::FIRST_HALF);
			std::vector<int> tps_idx = findTargetPoints(p, PopulationDivision::FIRST_HALF);
			tps.clear();
			for(int l = 0; l < tps_idx.size(); l++) {
				tps.push_back(p[tps_idx[l]]);
			}
		}
		p.sort_by_fit(PopulationDivision::FIRST_HALF);
		return tps;
	}

	/////////////////////////////////////

	int findMin(std::vector<int> & vec) {
		assert(vec.size() > 0);
		int idx = 0;
		for(int i = 1; i < vec.size(); i++) {
			idx = vec[i] < vec[idx] ? i : idx;
		}
		return idx;
	}

	int findMax(std::vector<int> & vec) {
		assert(vec.size() > 0);
		int idx = 0;
		for(int i = 1; i < vec.size(); i++) {
			idx = vec[i] > vec[idx] ? i : idx;
		}
		return idx;
	}
}