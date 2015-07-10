#include <vector>
#include <assert.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include "EVAL_HAMMER.h"
#include "individual.h"
#include <thread>
#include "hammer.h"
#include "util.h"

#pragma once

namespace optimization {

	/*
	class FTZL2 {
		std::vector<double> c;
		std::vector<double> s;

	public:
		 void eval(std::vector<double> & in, std::vector<double> & out) {
			assert(out.size() > 0);
			assert(out.size() <= in.size());
			c.resize(out.size());
			s.resize(out.size());
			double sum = 0;
			for(int i = out.size() - 1; i < in.size(); i++) {
				sum += (in[i] - 0.5) * (in[i] - 0.5);
			}
			c[0] = s[0] = 1;
			for(int i = 1; i < out.size(); i++) {
				c[i] = c[i-1] * cos(in[i-1] * 0.5 * M_PI);
				s[i] = sin(in[out.size()-1-i] * 0.5 * M_PI);
			}
			for(int i = 0; i < out.size(); i++) {
				out[i] = (1 + sum) * s[i] * c[out.size()-1-i];
			}
		}
	};

	class HAMMER {
	public:
		void eval(std::vector<double> & in, std::vector<double> & out) {
			::eval(in);
			out[0] = dkef > 0.01 ? 5*dkef : 0.0;
			out[1] = dfluxo > 0.01 ? 5*dfluxo : 0.0;
			out[2] = FP;
			out[3] = ktest < kef ? 5*dktest/0.03 : 0.0;
		}
	};*/

	template<int NUM>
	class HAMMER_MULT {
	
		HammerInterface * hm[NUM];
		int begin, end;
		std::vector<Individual> * pop;

		void run_thread(int idx);
	public:
		HAMMER_MULT(std::wstring & path = std::wstring(L"hammer"));
		void eval(std::vector<Individual> & pop, int begin, int end);
	};

	template<int NUM>
	void HAMMER_MULT<NUM>::run_thread(int idx) {
		for(int i = begin + idx; i < end; i+= NUM) {
			hm[idx]->eval((*pop)[i].x);
			
			(*pop)[i].fit[0] = hm[idx]->dkef > 0.01 ? 5*hm[idx]->dkef : 0.0;
			(*pop)[i].fit[1] = hm[idx]->dfluxo > 0.01 ? 5*hm[idx]->dfluxo : 0.0;
			(*pop)[i].fit[2] = hm[idx]->FP;
			(*pop)[i].fit[3] = hm[idx]->ktest < hm[idx]->kef ? 5*hm[idx]->dktest/0.03 : 0.0;
		}
	}

	template<int NUM>
	HAMMER_MULT<NUM>::HAMMER_MULT(std::wstring & path = std::wstring(L"hammer")) {
		for(int i = 0; i < NUM; i++) {
			std::wstring ws = std::to_wstring(i);
			_wmkdir(ws.c_str());
			CopyDirTo(path, ws);
			hm[i] = new HammerInterface(ws);
		}
	}

	template<int NUM>
	void HAMMER_MULT<NUM>::eval(std::vector<Individual> & _pop, int _begin, int _end) {
		begin = _begin;
		end = _end;
		pop = &_pop;
		std::thread threads [NUM];
		for(int i = 0; i < NUM; i++) {
			threads[i] = std::thread(&HAMMER_MULT<NUM>::run_thread, this, i);
		}
		for(int i = 0; i < NUM; i++) {
			threads[i].join();
		}
	}
}