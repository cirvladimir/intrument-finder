#include "matrix.h"
#include "polynomial.h"
#include <iostream>
#include <cmath>
#include <sstream>
#include <string>

#define _USE_MATH_DEFINES

using namespace std;

void moveDownZeros(double ** mat, int maxRows, int colNum, int maxCols)
{
	int dRow = 0;
	while ((dRow < maxRows) && (*(*(mat + dRow) + colNum) == 0))
	{
		dRow++;
	}
	if (dRow == maxRows)
	{
		//If colNum < maxCols - 1, recursion for next col
		if (colNum < maxCols - 1)
			moveDownZeros(mat, maxRows, colNum + 1, maxCols);
	}
	else
	{
		//swap dRow with row zero
		double * tempRow = *mat;
		*mat = *(mat + dRow);
		*(mat + dRow) = tempRow;
		//if maxRows > 2 and colNum < maxCols - 1, recursion for next col + row
		if ((maxRows > 2) && (colNum < maxCols - 1))
			moveDownZeros(mat + 1, maxRows - 1, colNum + 1, maxCols);
	}
}

//row r1 -> r1 - a * r2 
void subtMultRow(double ** mat, int r1, int r2, double a, int cols)
{
	for (int i = 0; i < cols; i++)
	{
		*(*(mat + r1) + i) -= a * *(*(mat + r2) + i);
	}
}

void reduceDir(double ** mat, int dRow, int times, int col, int cols)
{
	if (*(*mat + col) != 0)
	{
		for (int i = 1; i <= times; i++)
		{
			subtMultRow(mat, i * dRow, 0, *(*(mat + i * dRow) + col) /
				*(*mat + col), cols);
		}
	}
}

//r1 -> r1 * a
void multRow (double ** mat, double a, int cols)
{
	for (int i = 0; i < cols; i++)
	{
		*(*mat + i) *= a;
	}
}


Eigenpair::Eigenpair(double val, double * vec)
{
	value = val;
	vector = vec;
}

Matrix::Matrix() {}

Matrix::Matrix(int rs, int cls)
{
	rows = rs;
	cols = cls;
	//double * tmp[rows];
	mat = new double * [rows];
	for (int row = 0; row < rows; row++)
	{
		mat[row] = new double[cols];
		for (int col = 0; col < cols; col++)
		{
			mat[row][col] = 0;
		}
	}
}

Matrix Matrix::copy()
{
	Matrix cop;
	cop.rows = rows;
	cop.cols = cols;
	cop.mat = new double * [rows];
	for (int row = 0; row < rows; row++)
	{
		cop.mat[row] = new double[cols];
		for (int col = 0; col < cols; col++)
		{
			cop.mat[row][col] = mat[row][col];
		}
	}
	return cop;
}

void Matrix::set(int r, int c, double val)
{
	mat[r][c] = val;
}

double Matrix::get(int r, int c)
{
	return mat[r][c];
}

void Matrix::add(int r, int c, double val)
{
	mat[r][c] += val;
}

void Matrix::print()
{
	cout << "------" << endl;
	for (int r = 0; r < rows; r++)
	{
		for (int c = 0; c < cols; c++)
		{
			stringstream strStr;
			strStr << round(pow(10, PRINT_SIZE + 2) * mat[r][c]) / pow(10, PRINT_SIZE + 2);
			string str = strStr.str();
			if (str.size() > PRINT_SIZE) {
				str = str.substr(0, PRINT_SIZE);
				if (str.at(PRINT_SIZE - 1) == '.')
					str = str.substr(0, PRINT_SIZE - 1);
			}
			cout << str;
			for (int i = 0; i < PRINT_SIZE + 1 - str.size(); i++)
			{
				cout << " ";
			}
		}	
		cout << endl;
	}
	cout << "------" << endl;
}

void Matrix::rowReduce()
{
	int minDim = min(rows, cols);
	int nextCol = 0;
	int curRow = 0;
	//get into eschelon form
	while (curRow < minDim - 1)
	{
		//make sure all zeroes go down
		moveDownZeros(mat, rows, 0, cols);
		while ((nextCol < cols) && (*(*(mat + curRow) + nextCol) == 0))
			nextCol++;
		if (nextCol == cols)
			break;
		reduceDir(mat + curRow, 1, rows - curRow - 1, nextCol, cols);
		curRow++;
	}
	curRow--;
	
	//finish row-reduced eschelon
	while (curRow >= 0)
	{
		while ((nextCol > 0) && (*(*mat + nextCol) == 0))
			nextCol--;
		multRow(mat + curRow, 1 / *(*(mat + curRow) + nextCol), cols);
		reduceDir(mat + curRow, -1, curRow, nextCol, cols);
		curRow--;
	}
}


//finds null vector of row-reduced eschelon-form matrix=
void findNullVector(double ** mat, int rows, int cols, double * vec)
{
	//go through each column, keeping track of current pivot row
	//if value at row, col is 1, then col is pivot column
	//	store - summ of all other entries in row into the cell, inc pivot row
	//else store 1
	int row = 0;
	for (int col = 0; col < cols; col++)
	{
		if (row >= rows)
		{
			*vec = 1;
		}
		else if (abs(*(*(mat + row) + col) - 1) < ZERO_THRESH)
		{
			*vec = 0;
			//loop through all columns in pivot row, get negative 
			for (int chkCol = col + 1; chkCol < cols; chkCol++)
			{
				*vec -= *(*(mat + row) + chkCol);
			}
			row++;
		}
		else
		{
			*vec = 1;
		}
		vec++;
	}
}

void printArrr(double * ar, int dim)
{
	for (int i = 0; i <dim; i++)
	{
		cout << ar[i] << ", ";
	}
	cout << endl;
}
void Matrix::getEigenvector(double val, double * iVec)
{
	Matrix copMat = copy();
	for (int i = 0; i < rows; i++)
	{
		copMat.add(i, i, -val);
	}
	copMat.rowReduce();
	findNullVector(copMat.mat, rows, cols, iVec);
}

vector<Eigenpair> Matrix::getEigenpairs()
{
	Matrix copMat = copy();
	for (int i = 0; i < rows; i++)
	{
		for (int j = i + 1; j < rows; j++)
		{
			double th;
			if (copMat.mat[i][i] == copMat.mat[j][j])
				th = M_PI / 4;
			else
				th = 0.5 * atan(2 * copMat.mat[i][j] / (copMat.mat[j][j] - copMat.mat[i][i]));
			double c = cos(th);
			double s = sin(th);
			double sii = copMat.mat[i][i];
			double sjj = copMat.mat[j][j];
			double sij = copMat.mat[i][j];
			copMat.mat[i][i] = c * c * sii - 2 * s * c * sij + s * s * sjj;
			copMat.mat[j][j] = c * c * sii + 2 * s * c * sij + s * s * sjj;
			copMat.mat[i][j] = 0;
			copMat.mat[j][i] = copMat.mat[i][j];
			for (int k = i + 1; k < cols; k++)
			{
				if ((k != i) && (k != j))
				{
					double sik = copMat.mat[i][k];
					double sjk = copMat.mat[j][k];
					copMat.mat[i][k] = c * sik - s * sjk;
					copMat.mat[k][i] = copMat.mat[i][k];
					copMat.mat[j][k] = s * sik + c * sjk;
					copMat.mat[k][j] = copMat.mat[j][k];
				}
			}
		}
	}
	vector<Eigenpair> prs;
	for (int i = 0; i < rows; i++)
	{
		double eigVal = copMat.mat[i][i];
		double * eigVec = new double[rows];
		getEigenvector(eigVal, eigVec);
		prs.push_back(Eigenpair(eigVal, eigVec));
	}
	return prs;	
}

