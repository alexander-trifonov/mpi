#include "mpi.h"
#include <iostream>
#include <stdlib.h>
#include <time.h>
using namespace std;

int GenerateMatrix(double* &matrix, int rows, int cols, int min_range, int max_range)
{
	srand(time(NULL));
	matrix = (double*)malloc(cols*rows*sizeof(double));
	for (int i = 0; i < cols*rows; i++)
	{
	matrix[i] = rand() % max_range + min_range;
	}
	return 0;
}

void PrintMatrix(double* &matrix, int rows, int cols)
{
	for (int i = 0; i < cols; i++)
	{
	for (int j = 0; j < rows; j++)
	cout << matrix[j + (cols*i)] << " ";
	cout << endl;
	}
}

double* FindMin(double* &matrix, int Rows, int Cols, int RankSize)
{
	double* min = (double*)malloc(Rows/RankSize * sizeof(double));
	for (int i = 0; i < Rows/RankSize; i++)
		min[i] = 65000;

	for (int i = 0; i < (Rows / RankSize); i++)
	{
		for (int j = 0; j < Cols; j++)
		{
			if (matrix[j + i*Cols] < min[i])
				min[i] = matrix[j + i*Cols];
		}
	}
	return min;
}

void main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	int rank, RankSize;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &RankSize);

	
	double* matrix;
	int Rows = 6000;
	int Cols = 6000;
	if (rank == 0)
	{
		cout << "[[ Generating matrix ]]" << endl;
		/*
		cout << "Rows: ";
		cin >> Rows;
		cout << "Cols: ";
		cin >> Cols;
		*/
		GenerateMatrix(matrix, Rows, Cols, 0, 10);
		//PrintMatrix(matrix, Rows, Cols);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	double TimeStart;
	double TimeEnd;
	if (RankSize > 1)
	{
		double* min;
		double* minRecv = (double*)malloc(Rows * sizeof(double));
		if (rank == 0)
		{
			cout << "[[ MPI ]]" << endl;
			for (int i = 1; i < RankSize; i++)
			{
				MPI_Send(&matrix[i*(Cols*(Rows / RankSize))], Cols*(Rows / RankSize), MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
				TimeStart = MPI_Wtime();
				min = FindMin(matrix, Rows, Cols, RankSize);
			}
		}
		else
		{
			double *RecvBuf = (double*)malloc(Cols*(Rows / RankSize) * sizeof(double));
			MPI_Recv(RecvBuf, Cols*(Rows / RankSize), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			min = FindMin(RecvBuf, Rows, Cols, RankSize);
			free(RecvBuf);
		}
		MPI_Gather(min, Rows / RankSize, MPI_DOUBLE, minRecv, Rows / RankSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if (rank == 0)
		{
			TimeEnd = MPI_Wtime();
			/*
			cout << "Rank 0 received all mins: " << endl;
			for (int i = 0; i < Rows; i++)
				cout << i << ": " << minRecv[i] << endl;
				*/
			cout << "Time: " << TimeEnd - TimeStart << endl;
		}
		free(minRecv);
		free(min);

		if (rank == 0)
		{
			cout << "[[ Sequential ]]" << endl;
			TimeStart = MPI_Wtime();
			min = FindMin(matrix, Rows, Cols, 1);
			TimeEnd = MPI_Wtime();
			cout << "Time: " << TimeEnd - TimeStart << endl;
		}

	}
	MPI_Finalize();
}