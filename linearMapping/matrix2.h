#ifndef Matrix_h
#define Matrix_h

#define ZERO_THRESH 0.00000000001
#define PRINT_SIZE 6

#include <vector>

using namespace std;

class Eigenpair
{
	public:
		Eigenpair(double val, double * vec);
		double value;
		double * vector;
};

class Matrix 
{
	private:
	 	double ** mat;
	 	int rows, cols;
	 	Matrix();
	 	void getEigenvector(double val, double * iVec);
	public:
		Matrix(int rs, int cls);
		void set(int row, int col, double val);
		double get(int row, int col);
		void rowReduce();
		void add(int row, int col, double val);
		void print();
		vector<Eigenpair> getEigenpairs();
		Matrix copy();
};

#endif
