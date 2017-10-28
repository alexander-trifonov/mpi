///*
//	Task 17: Searching mins in matrix rows;
//	Restrictions: # process%Cols = 0;
//*/
//#include "mpi.h"
//#include <iostream>
//#include <stdlib.h>
//#include <time.h>
//#include <cstdlib>
//using namespace std;
//
//int GenerateMatrix(double* &matrix, int rows, int cols, int min_range, int max_range)
//{
//	srand(time(NULL));
//	matrix = (double*)malloc(cols*rows * sizeof(double));
//	for (int i = 0; i < cols*rows; i++)
//	{
//		matrix[i] = rand() % max_range + min_range;
//		//cout << i/((rows*cols) / 100) <<"%" << '\r'<<flush;
//	}
//	return 0;
//}
//
//void PrintMatrix(double* &matrix, int rows, int cols)
//{
//	for (int i = 0; i < rows; i++)
//	{
//		for (int j = 0; j < cols; j++)
//			cout << matrix[j + (cols*i)] << " ";
//		cout << endl;
//	}
//}
//
//double* FindMin(double* &matrix, int Cols, int DataSize)
//{
//	double* min = (double*)malloc(DataSize * sizeof(double));
//	for (int i = 0; i < DataSize; i++)
//		min[i] = 65000;
//
//	for (int i = 0; i < (DataSize); i++)
//	{
//		for (int j = 0; j < Cols; j++)
//		{
//			if (matrix[j + i*Cols] < min[i])
//			{
//				min[i] = matrix[j + i*Cols];
//			}
//		}
//	}
//	return min;
//}
//
//bool AreEqual(double* &array1, double* &array2, int Rows)
//{
//	for (int i = 0; i < Rows; i++)
//	{
//		if (array1[i] != array2[i])
//			return false;
//	}
//	return true;
//}
//void main17(int argc, char** argv)
//{
//	MPI_Init(&argc, &argv);
//	int rank, RankSize;
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//	MPI_Comm_size(MPI_COMM_WORLD, &RankSize);
//	double TimeStart, TimeEnd;
//
//	double* matrix;	// a main matrix, init in rank 0;
//	int Rows = atoi(argv[1]);	// # rows in the matrix
//	int Cols = Rows;// # cols in the matrix
//	int DataSize;	// # rows to send each process
//	double *RecvBuf;// Matrix for FindMin function
//	double* minSeq;
//	if (rank == 0)
//	{
//		//[[ Matrix initialization ]]
//		cout << "[[ Generating matrix ]]" << endl;
//		cout << "[ Rows: " << Rows << " Cols: " << Cols << " ]" << endl;
//		GenerateMatrix(matrix, Rows, Cols, 0, 10);
//		//PrintMatrix(matrix, Rows, Cols);
//
//		//[[ Sequential ]]
//		cout << "[[ Sequential ]]" << endl;
//		TimeStart = MPI_Wtime();
//		minSeq = FindMin(matrix, Rows, Cols);
//		TimeEnd = MPI_Wtime();
//		//PrintMatrix(minSeq, Rows, 1); //Uncomment if we want to show mins;
//		cout << "Time: " << TimeEnd - TimeStart << endl;
//
//		//[[ MPI ]]
//		cout << "[[ MPI ]]" << endl;
//		TimeStart = MPI_Wtime();
//		DataSize = Rows / RankSize;
//	}
//
//	// [[ Shared code ]]
//	MPI_Bcast(&DataSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
//
//	// [[ Private code ]]
//	if(rank == 0)
//	{
//		for (int i = 1; i < RankSize; i++)
//		{
//			MPI_Send(&matrix[i*(Cols*DataSize)], Cols*DataSize, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
//		}
//		RecvBuf = matrix;
//	}
//	else
//	{
//		RecvBuf = (double*)malloc(Cols*(DataSize) * sizeof(double));
//		MPI_Recv(RecvBuf, Cols*(DataSize), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//	}
//	// [[ Shared code ]]
//	double* minMPI_datasize = FindMin(RecvBuf, Cols, DataSize); // Array.length = DataSize; Contains mins
//	double* minMPI = (double*)malloc(Rows * sizeof(double));
//	MPI_Gather(minMPI_datasize, DataSize, MPI_DOUBLE, minMPI, Rows / RankSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//	TimeEnd = MPI_Wtime();
//
//	// [[ Private code ]]
//	if (rank == 0)
//	{
//		//PrintMatrix(minRecv, Rows, 1); //Uncomment if we want to show mins;
//		cout << "Time: " << TimeEnd - TimeStart << endl;
//		if (!AreEqual(minSeq, minMPI, Rows))
//		{
//			cout << "Error: minSeq != minMPI" << endl;
//		}
//		else
//		{
//			cout << "Success: minSeq == minMPI" << endl;
//		}
//	}
//
//	free(minSeq);
//	free(minMPI);
//	free(minMPI_datasize);
//	MPI_Finalize();
//}