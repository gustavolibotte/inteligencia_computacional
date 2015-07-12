#include <vector>
#include "individual.h"
#include <assert.h>
#include <cmath>
#include "defs.h"

#pragma once
namespace optimization {

	class Bounds {

	public:
		Bounds(int d_dim, int i_dim) {
			lower_bounds = std::vector<double>(d_dim);
			upper_bounds = std::vector<double>(d_dim);
			lower_fit = std::vector<double>(i_dim);
			upper_fit = std::vector<double>(i_dim);
		}
		std::vector<double> lower_bounds;
		std::vector<double> upper_bounds;
		std::vector<double> lower_fit;
		std::vector<double> upper_fit;

		void ensureBounds(std::vector<Individual>::iterator begin, std::vector<Individual>::iterator end) {
			assert(begin != end);
			assert(lower_bounds.size() == (*begin).size() || upper_bounds.size() == (*begin).size());
			std::vector<double> range(lower_bounds.size());
			for(int i = 0; i < lower_bounds.size(); i++) {
				range[i] = upper_bounds[i] - lower_bounds[i];
			}

			for(auto p = begin; p !=end; p++) {
				for(int i = 0; i < (*p).size(); i++) {
					assert(isfinite((*p)[i]));
					if((*p)[i] < lower_bounds[i] || (*p)[i] > upper_bounds[i]) {
						double replace = fmod((*p)[i] - lower_bounds[i], range[i]);
						replace = replace > 0 ? replace: range[i] + replace;
						(*p)[i] = replace + lower_bounds[i];
						assert((*p)[i] >= lower_bounds[i] && (*p)[i] <= upper_bounds[i]);
					}
					assert(isfinite((*p)[i]));
				}
			}
		}
	};
}