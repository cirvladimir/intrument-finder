#ifndef data_h
#define data_h

#include "matrix.h"

#include <vector>

using namespace std;

class Data
{
	private:
		int dim;
		double * avgVec;
		vector<Eigenpair> egPrs;
		vector<double *> basis;
	public:
		Data(int dimension);
		void load(double * points, int numPts);
		vector<Eigenpair> getPairs();
		void compOrthoBasis(double minIVal);
		vector<double *> getOrthoBasis();
		void getProjection(double * vec, double * proj);
		int match(double ** vecs, int numVs, double * v);
};

#endif
