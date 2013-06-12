#include "data.h"
#include <iostream>
#include <cmath>

using namespace std;

Data::Data(int d)
{
	dim = d;
}

vector<Eigenpair> Data::getPairs()
{
	return egPrs;
}
void printAr(double * ar, int dim)
{
	for (int i = 0; i <dim; i++)
	{
		cout << ar[i] << ", ";
	}
	cout << endl;
}

void Data::load(double * points, int numPts)
{
	egPrs.empty();
	//get avg for each coord
	avgVec = new double[dim];
	for (int i = 0; i < dim; i++)
		avgVec[i] = 0;
	for (int ptInd = 0; ptInd < numPts; ptInd++)
	{
		for (int dmInd = 0; dmInd < dim; dmInd++)
		{
			avgVec[dmInd] += *(points + ptInd * dim + dmInd) / numPts;
		}
	}
	//printAr(avgVec, dim);
		
	//compute covariance matrix
	Matrix covM(dim, dim);
	for (int ptInd = 0; ptInd < numPts; ptInd++)
	{
		for (int i = 0; i < dim; i++)
		{
			//cout << points[ptInd * dim + i] << endl;
			for (int j = 0; j < dim; j++)
			{
				//cout << (*(points + ptInd * dim + i) - avgVec[i]) *
				//	(*(points + ptInd * dim + j) - avgVec[j]) << " ";
				covM.add(i, j, (*(points + ptInd * dim + i) - avgVec[i]) *
					(*(points + ptInd * dim + j) - avgVec[j]));
			}
		}
	}
	
	//covM.set(0, 0, 0.0575962);        covM.set(0, 1,  0.0102312  );       
	//covM.set(1, 0, 0.0102312   );  covM.set(1, 1,    0.00181744);
	
	//covM.print();
	
	//return eigen pairs for covariance matrix
	egPrs = covM.getEigenpairs();
}

vector<double *> Data::getOrthoBasis()
{
	return basis;
}

double dot(double * v1, double * v2, int dim)
{
	double ret = 0;
	for (int i = 0; i < dim; i++)
	{
		ret += *(v1++) * *(v2++);
	}
	return ret;
}

void normVect(double * v, int dim, double * outp)
{
	double len = sqrt(dot(v, v, dim));
	for (int i = 0; i < dim; i++)
	{
		outp[i] = v[i] / len;
	}
}

void Data::compOrthoBasis(double minIVal)
{
	basis.empty();
	for (int i = 0; i < egPrs.size(); i++)
	{
		if (egPrs.at(i).value >= minIVal)
		{
			double * bs = new double[dim];
			normVect(egPrs.at(i).vector, dim, bs);
			//printAr(bs, dim);
			basis.push_back(bs);
		}
	}
}

double getDist(double * v1, double * v2, int dim)
{
	double diff[dim];
	for (int i = 0; i < dim; i++)
	{
		diff[i] = v1[i] - v2[i];
	}
	return sqrt(dot(diff, diff, dim));
}

double getAng(double * v1, double * v2, int dim)
{
	double * origin = new double[dim];
	for (int i = 0; i < dim; i++)
	{
		origin[i] = 0;
	}
	double v1Norm[dim];
	double v2Norm[dim];
	normVect(v1, dim, v1Norm);
	normVect(v2, dim, v2Norm);
	double cosine = dot(v1Norm, v2Norm, dim);
	return abs(cosine - 1);
}


int Data::match(double ** vecs, int numVs, double * v)
{
	if (numVs == 0)
	{
		return -1;
	}
	else
	{
		double vProj[basis.size()];
		getProjection(v, vProj);
		double tProj[basis.size()];
		getProjection(vecs[0], tProj);
		double minDist = getAng(vProj, tProj, basis.size());//getDist(vProj, tProj, basis.size());
		int minInd = 0;
		for (int i = 1; i < numVs; i++)
		{
			getProjection(vecs[i], tProj);
			double curDist = getAng(vProj, tProj, basis.size());//getDist(vProj, tProj, basis.size());
			if (curDist < minDist)
			{
				minDist = curDist;
				minInd = i;
			}
		}
		return minInd;
	}
}

void Data::getProjection(double * vec, double * proj)
{
	double vecMAvg[dim];
	for (int i = 0; i < dim; i++)
	{
		vecMAvg[i] = vec[i] - avgVec[i];
	}
	for (int bInd = 0; bInd < basis.size(); bInd++)
	{
		proj[bInd] = dot(basis[bInd], vecMAvg, dim);
	}
}
