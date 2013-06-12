#ifndef Polynomial_h
#define Polynomial_h

#define BIN_PRECISION 0.00000000000000000001
#define MAX_BIN_ITS 100

#include <vector>

class Polynomial {
	private:
		double * coefs;
		int numCoefs;
		double binarySearch(double min, double max, bool inc, int its);
	public:
		Polynomial(double * cfs, int numCfs);
		std::vector<double> getRoots();
		double getValue(double x);
};

#endif
