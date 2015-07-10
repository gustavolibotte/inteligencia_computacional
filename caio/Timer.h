#include "defs.h"
#include <Windows.h>
#include <iostream>

#pragma once

class Timer {
	long t;
	long l;
public:
	Timer() {}

	void tic() {
		t = GetTickCount();
	}

	void toc() {
		 l = GetTickCount() - t;
	}

	void print() {
		int ms = l % 1000;
		l /= 1000;
		int sec =  l % 60;
		l /= 60;
		int min = l % 60;
		l /= 60;
		int hour = l % 24;
		std::cout << "Tempo total: " << hour << ":" << min << ":" << sec << ".:." << ms << std::endl;
	}
};