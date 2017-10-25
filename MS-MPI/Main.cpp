#include "mpi.h"
#include <iostream>
#include <stdlib.h>
#include <time.h>
using namespace std;

int GenerateMatrix(double** &matrix, int rows, int cols, int min_range, int max_range)
{
	srand(time(NULL));
	/*
	matrix = (double*)malloc(cols*rows*sizeof(double));
	for (int i = 0; i < cols*rows; i++)
	{
		matrix[i] = rand() % max_range + min_range;
	}
	*/
	matrix = new double*[rows];
	for (int i = 0; i < rows; i++)
	{
		matrix[i] = new double[cols];
		for(int j = 0; j < cols; j++)
			matrix[i][j] = rand() % max_range + min_range;
	}
	return 0;
}

void PrintMatrix(double** &matrix, int rows, int cols)
{
	/*
	for (int i = 0; i < cols; i++)
	{
		for (int j = 0; j < rows; j++)
			cout << matrix[j + (cols*i)] << " ";
		cout << endl;
	}
	*/
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			cout << matrix[i][j] << " ";
		}
		cout << "\n";
	}
}

void main3(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	int rank, RankSize;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &RankSize);
	double** matrix;
	int Rows = 6;
	int Cols = 6;
	if (rank == 0)
	{
		GenerateMatrix(matrix, Rows, Cols, 0, 10);
		PrintMatrix(matrix, Rows, Cols);

		if (RankSize > 1)
		{
			//Create data blocks for future transport
			int odd = 0; //If odd, then add one row to main process
			int RowsInSingleBlock = Rows / RankSize;
			cout << "RowsInSingleBlock: " << RowsInSingleBlock << endl;
			if (RowsInSingleBlock % 2 != 0)
				odd = 1;
			RowsInSingleBlock += odd;

			double* Data_Buf;
			for (int i = 1; i < RankSize; i++)
			{
				//double* Data_Buf_Cont = (double*)malloc(RowsInSingleBlock * Cols * sizeof(double));
				Data_Buf = (double*)malloc(RowsInSingleBlock * Cols * sizeof(double*)); //may be 1 or more rows per data block for transfer
				for (int j = 0; j < RowsInSingleBlock; j++)
				{
					for (int k = 0; k < Cols; k++)
					{
						Data_Buf[k + j*Cols] = matrix[j+(i-1)*RowsInSingleBlock][k];
						cout << "Data_Buf["<< k + j*Cols <<"]: " << Data_Buf[k + j*Cols] << endl;
					}
				}
				cout << "Sending to.. " << i << endl;
				MPI_Send(&Data_Buf[0], Cols*RowsInSingleBlock, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
				int b = 1;
				//MPI_Send(&b, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
				cout << "Sedned to " << i << endl;
			}
		}
		cin.get();
		MPI_Finalize();
	}

	if (rank != 0)
	{
		double* Recv_Buf = (double*)malloc(1 * Cols * sizeof(double*));
		int a;
		cout << "I am rank: " << rank << endl;
		MPI_Recv(&Recv_Buf[0], Cols*(Rows / RankSize), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
		//MPI_Recv(&a, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		cout << "Rank: " << rank << " received: " << endl;
		for(int i = 0; i < (Rows / RankSize) * Cols; i++)
			cout << Recv_Buf[i] << " ";
		//PrintMatrix(Recv_Buf, Rows, Rows / RankSize);
		cout << "\n";
	}

		// Unfinished realization with pointers.
		/*
		//Unfinished realization with malloc and scatter
		double* matrix;
		int Rows = 3;
		int Cols = 3;
		double* data_buf = (double*)malloc(sizeof(double)*Cols);
		if (rank == 0)
		{
		GenerateMatrix(matrix, Rows, Cols, 0, 10);
		PrintMatrix(matrix, Rows, Cols);
		if (RankSize > 1)
		{
		MPI_Scatter(matrix, Cols, MPI_DOUBLE, data_buf, Cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		}
		}
		*/
	
}