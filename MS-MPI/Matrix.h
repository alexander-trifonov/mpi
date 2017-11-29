#pragma once
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <cstdlib>
using namespace std;
int GenerateMatrix(double* &matrix, int rows, int cols, int min_range, int max_range)
{
	srand(time(NULL));
	matrix = (double*)malloc(cols*rows * sizeof(double));
	for (int i = 0; i < cols*rows; i++)
	{
		
		matrix[i] = rand() % max_range + min_range;
		//cout << i/((rows*cols) / 100) <<"%" << '\r'<<flush;
	}
	return 0;
}

void PrintMatrix(double* &matrix, int rows, int cols)
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
			cout << matrix[j + (cols*i)] << "\t";
		cout << endl;
	}
}
void PrintMatrixRank(double* &matrix, int rows, int cols, int rank)
{
	for (int i = 0; i < rows; i++)
	{
		cout << "r:"<<rank<<" "<<"i: "<<i<<"\t";
		for (int j = 0; j < cols; j++)
			cout << matrix[j + (cols*i)] << "\t";
		cout << endl;
	}
}

double* FindMin(double* &matrix, int Cols, int DataSize)
{
	double* min = (double*)malloc(DataSize * sizeof(double));
	for (int i = 0; i < DataSize; i++)
		min[i] = 65000;

	for (int i = 0; i < (DataSize); i++)
	{
		for (int j = 0; j < Cols; j++)
		{
			if (matrix[j + i*Cols] < min[i])
			{
				min[i] = matrix[j + i*Cols];
			}
		}
	}
	return min;
}

bool AreEqual(double* &array1, double* &array2, int Rows, int Cols)
{
	for (int i = 0; i < Rows; i++)
	{
		for (int j = 0; j < Cols; j++)
		{
			if (array1[i*Cols + j] != array2[i*Cols + j])
			{
				cout << array1[i*Cols + j] << " ?= " << array2[i*Cols + j] << endl;
				return false;
			}
		}

	}
	return true;
}

void CopyInto(double* &array1, double* &array2, int Rows, int Cols)
{
	for (int i = 0; i < Rows; i++)
	{
		for (int j = 0; j < Cols; j++)
			array2[i*Cols + j] = array1[i*Cols + j];
	}
}