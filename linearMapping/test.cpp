#include <iostream>
#include <cmath>
#include <cstdlib>
#include "data.h"
#include "polynomial.h"

#define NUM_POINTS 10
#define _USE_MATH_DEFINES

using namespace std;

class Point
{
	public:
		Point (double, double);
		double x, y;
};

Point::Point(double _x, double _y)
{
	x = _x;
	y = _y;
}

void printArr(double * ar, int dim)
{
	for (int i = 0; i <dim; i++)
	{
		cout << ar[i] << ", ";
	}
	cout << endl;
}


int main() 
{
	int w = 3;
	int h = 3;
	Matrix m(w, h);
	for (int i = 0; i < w; i ++)
	{
		for (int j = i; j < h; j++)
		{
			double nextV = rand() % 100 + 1;
			m.set(i, j, nextV);
			m.set(j, i, nextV);
		}
	}
	
	m.print();
	//m.multByRot(1, 2, 0.7144496363095);
	//m.simpMat(1, 2);
	//m.print();
	cout << "working" << endl;
	vector<Eigenpair> prs = m.getEigenpairs();
	for (int i = 0; i < prs.size(); i++)
	{
		cout << "value: " << prs.at(i).value;
		cout << " vector: ";
		printArr(prs.at(i).vector, 2);
	}
	cout << "done" << endl;
	
	/*
	int numPts = 40;
	int numDims = 5;
	double dimRange[5][2] = {{0, 10}, {0, 1}, {0, 10}, {0, 1}, {0, 1}};
	double pts[numPts * numDims];
	for (int i = 0; i < numPts; i++)
	{
		for (int j = 0; j < numDims; j++)
		{
			double randD = (rand() % 100000) * 1.0 / 100000;
			pts[i * numDims + j] = (dimRange[j][1] - dimRange[j][0]) * randD + dimRange[j][0];
			//cout << pts[i * numDims + j] << endl;
		}
	}
	Data dt(numDims);
	dt.load(pts, numPts);
	vector<Eigenpair> egPrs = dt.getPairs();
  for (int i = 0; i < egPrs.size(); i++)
  {
		cout << "eigenvalue: " << egPrs.at(i).value << endl;
		cout << "eigenvector: ";
		for (int j = 0; j < numDims; j++)
		{
			cout << egPrs.at(i).vector[j];
			if (j != numDims - 1) 
				cout << ", ";
		}
		cout << endl;
	}
	
	double avgEig = 0;
	for (int i = 0; i < egPrs.size(); i++)
	{
		avgEig += egPrs.at(i).value / egPrs.size();
	}
	
	dt.compOrthoBasis(avgEig);
	vector<double *> bs = dt.getOrthoBasis();
	for (int i = 0; i < bs.size(); i++)
	{
		for (int j = 0; j < numDims; j++)
		{
			cout << bs.at(i)[j] << ", ";
		}
		cout << endl;
	}
	
	double sampVecs[4][5] = { 
		{ 10, 20, 10, 10, 10},
		{7, .5, .5, -4, .5},
		{20, 20, -20, 20, 20},
		{-10, 2, 2, 3, 41}};
		
	double * sV[4] = { sampVecs[0], sampVecs[1], sampVecs[2], sampVecs[3] };
	
	double testVec[5] = {5.69204, 0.488711, 6.36872, 0.518744, 0.579536};
	
	//int bestM = dt.match(sV, 4, testVec);
	//cout << "best match " << bestM << endl;
	
	double * proj = new double[bs.size()];
	dt.getProjection(testVec, proj);
	for (int i = 0; i < bs.size(); i++)
	{
		cout << proj[i] << ", ";
	}
	cout << endl;
	
	/*
	
	double pl[7] = {1.842880000213704e-10, -0.059413640000000004, 1};
	Polynomial ps(pl, 3);
	vector<double> roots = ps.getRoots();
	for (int i = 0; i < roots.size(); i++)
		cout << roots.at(i) << ", ";
	cout << endl;
	
	double prepMat[3][3] = { { 4, 6, 0}, {4, 0, 1}, {2, 3, 0} };
	
	double * mat[3] = {prepMat[0], prepMat[1], prepMat[2]};
	
	printM(mat, 3, 3);
	cout << "-----" << endl;
	
	rowReduce(mat, 3, 3);
	
	printM(mat, 3, 3);
	
	double nullVec[3];
	findNullVector(mat, 3, 3, nullVec);
	double * vecSet[1] = {nullVec};
	cout << "eigen values: ";
	printM(vecSet, 1, 3);*/
	/*
	
	Point ps [NUM_POINTS] = {Point(-1, -1), Point(1, 1), Point(-2, -2), Point(2, 2), 
		Point(-4, -3), Point(4, 3), Point(-4, -5), Point(4, 5), Point(-5, -5), Point(5, 5)};
	double mat[2][2];
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			mat[i][j] = 0;
		}
	} 
	for (int i = 0; i < NUM_POINTS; i++)
	{
		mat[0][0] += ps[i].x * ps[i].x;
		mat[0][1] += ps[i].x * ps[i].y;
		mat[1][0] += ps[i].y * ps[i].x;
		mat[1][1] += ps[i].y * ps[i].y;
	}
	
	cout << "matrix: " << endl;
	double * ttt[2] = { mat[0], mat[1] };
	printM(ttt, 2, 2);
	
	double vals[2];
	solveQuad(1, -mat[0][0] - mat[1][1], mat[0][0] * mat[1][1] - 
		mat[0][1] * mat[1][0], vals);
	//cout << "eigenvalues: " << endl;
	//cout << vals[0] << ", " << vals[1] << endl;
	//cout << "eigen vectors: " << endl;
	for (int i = 0; i < 2; i++)
	{
		cout << "eigen value " << vals[i] << " vector: ";
		double eigVec[2];
		getEigenvector(mat, vals[i], eigVec);
		cout << eigVec[0] << ", " << eigVec[1] << endl;
	}
	cout << "hello" << endl;
	
	double polCfs[3] = { -1, 0, 1 };
	Polynomial pol(polCfs, 3);
	vector<double> roots = pol.getRoots();
	for (vector<double>::iterator it = roots.begin(); it != roots.end(); ++it)
    cout << *it << endl;*/
    
  /*Matrix mat(3, 3);
  mat.set(0, 0, 1);
  mat.set(0, 1, 2);
  mat.set(0, 2, 3);
  mat.set(1, 0, 2);
  mat.set(1, 1, 3);
  mat.set(1, 2, 4);
  mat.set(2, 0, 3);
  mat.set(2, 1, 4);
  mat.set(2, 2, 5);
  vector<Eigenpair> egPrs = mat.getEigenpairs();
  for (int i = 0; i < egPrs.size(); i++)
  {
		cout << "eigenvalue: " << egPrs.at(i).value << endl;
		cout << "eigenvector: " << egPrs.at(i).vector[0] << ", " <<
			egPrs.at(i).vector[1] << endl;
	}*/
	
}
