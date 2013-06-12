#include "polynomial.h"
#include <cmath>
#include <iostream>

using namespace std;

double Polynomial::binarySearch(double min, double max, bool inc, int its)
{
	double mid = (min + max) / 2;
	//cout << "min: (" << min << ", " << getValue(min) << "), max: (" <<
	//	max << ", " << getValue(max) << ")" << mid<< endl;
	double val = getValue(mid);
	if (its > MAX_BIN_ITS)
		return mid;
	else if (abs(val) < BIN_PRECISION)
		return mid;
	else if ((val > 0) != inc) 
		return binarySearch(mid, max, inc, its + 1);
	else
		return binarySearch(min, mid, inc, its + 1);
}

Polynomial::Polynomial(double * cfs, int numCfs)
{
	coefs = cfs;
	numCoefs = numCfs;
}

double Polynomial::getValue(double x)
{
	double val = 0;
	for (int pow = 0; pow < numCoefs; pow++)
	{
		val += coefs[pow] * std::pow(x, pow);
	}
	return val;
}

vector<double> Polynomial::getRoots()
{
	vector<double> roots;
	if (numCoefs == 1)
	{
		if (*coefs == 0)
			roots.push_back(0);
	}
	else if (numCoefs == 2)
	{
		if (*(coefs + 1) != 0)
		{
			roots.push_back(-*coefs / *(coefs + 1));
		}
		else
		{
			if (*coefs == 0)
				roots.push_back(0);
		}
	}
	else
	{
		//find critical points
		double dxdy[numCoefs - 1];
		for (int pow = 1; pow < numCoefs; pow++)
		{
			dxdy[pow - 1] = pow * coefs[pow];
		} 
		vector<double> critPts = Polynomial(dxdy, numCoefs - 1).getRoots();
		double val1 = getValue(critPts.at(0));
		
		//check on - end
		bool incAtNegInf = (coefs[numCoefs - 1] > 0) == ((numCoefs % 2) == 1);
		if (incAtNegInf == (val1 <= 0))
		{
			double endCrit = critPts.at(0);
			double nextVal = -1;
			while (getValue(endCrit + nextVal) * val1 > 0)
				nextVal *= 2;
			roots.push_back(binarySearch(endCrit + nextVal, endCrit, !incAtNegInf, 0));
		}
		
		//check between critical points
		for (int critInd = 0; critInd < critPts.size() - 1; critInd++)
		{
			double val2 = getValue(critPts.at(critInd + 1));
			if (val1 * val2 <= 0)
			{//points on opposite sides
				roots.push_back(binarySearch(critPts.at(critInd), critPts.at(critInd + 1), val2 > val1, 0)); 
			}
			val1 = val2;
		}
		
		//check on + end
		bool incAtPlusInf = coefs[numCoefs - 1] > 0;
		if (incAtPlusInf == (val1 <= 0))
		{
			double endCrit = critPts.at(critPts.size() - 1);
			double nextVal = 1;
			while (getValue(endCrit + nextVal) * val1 > 0)
				nextVal *= 2;
			roots.push_back(binarySearch(endCrit, endCrit + nextVal, incAtPlusInf, 0));
		}
	}
	return roots;
}
