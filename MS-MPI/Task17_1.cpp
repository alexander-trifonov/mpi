///*
//	Task 17: Searching mins in matrix rows;
//	Restrictions: # process%Cols = 0;
//*/
//#include "mpi.h"
//#include <iostream>
//#include <stdlib.h>
//#include <time.h>
//#include <cstdlib>
//#include "Matrix.h"
//using namespace std;
//
//
//void main(int argc, char** argv)
//{
//	MPI_Init(&argc, &argv);
//	int rank, RankSize;
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//	MPI_Comm_size(MPI_COMM_WORLD, &RankSize);
//	double TimeStart, TimeEnd;
//
//	double* minMPI;
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
//		minSeq = FindMin(matrix, Cols, Rows);
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
//		minMPI = (double*)malloc(Rows * sizeof(double));//only for root
//	}
//	else
//	{
//		RecvBuf = (double*)malloc(Cols*(DataSize) * sizeof(double));
//		MPI_Recv(RecvBuf, Cols*(DataSize), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//	}
//	// [[ Shared code ]]
//	double* minMPI_datasize = FindMin(RecvBuf, Cols, DataSize); // Array.length = DataSize; Contains mins
//	MPI_Gather(minMPI_datasize, DataSize, MPI_DOUBLE, minMPI, Rows / RankSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//
//	// [[ Private code ]]
//	if (rank == 0)
//	{
//		TimeEnd = MPI_Wtime();//only for root
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
//		free(minMPI);
//		free(matrix);
//	}
//	else
//		free(RecvBuf);
//	free(minSeq);
//	free(minMPI_datasize);
//	MPI_Finalize();
//}



///*
//Task 17: Searching mins in matrix rows;
//Restrictions: # process%Cols = 0;
//*/
//#include "mpi.h"
//#include <iostream>
//#include <stdlib.h>
//#include <time.h>
//#include <cstdlib>
//#include "Matrix.h"
//using namespace std;
//
//
//
//void main(int argc, char** argv)
//{
//	//[ Matrix initialization ]
//	double* matrix;	// a main matrix, init in rank 0;
//	int Rows = atoi(argv[1]);	// # rows in the matrix
//	int Cols = Rows + 0;// # cols in the matrix
//
//	MPI_Init(&argc, &argv);
//	int rank, RankSize;
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//	MPI_Comm_size(MPI_COMM_WORLD, &RankSize);
//
//	double TimeStart, TimeEnd;
//	double* minMPI;
//	int DataSize;	// # rows process has
//	int* DataSizeArray; // # elems to send each process
//	int* RowSizeArray = new int[RankSize]; //# rows to send each process
//	int* Displ;
//	int* RowSum = new int[RankSize];
//	double *RecvBuf;// Matrix for FindMin function
//	double* minSeq;
//
//	if (rank == 0)
//	{
//		//[[ Matrix initialization ]]
//		cout << "[[ Generating matrix ]]" << endl;
//		cout << "[ Rows: " << Rows << " Cols: " << Cols << " ]" << endl;
//		GenerateMatrix(matrix, Rows, Cols, 0, 10);
//		PrintMatrix(matrix, Rows, Cols);
//
//		//[[ Sequential ]]
//		cout << "[[ Sequential ]]" << endl;
//		TimeStart = MPI_Wtime();
//		minSeq = FindMin(matrix, Cols, Rows);
//		TimeEnd = MPI_Wtime();
//		//PrintMatrix(minSeq, Rows, 1); //Uncomment if we want to show mins;
//		cout << "Time: " << TimeEnd - TimeStart << endl;
//
//		//[[ MPI ]]
//		cout << "[[ MPI ]]" << endl;
//		TimeStart = MPI_Wtime();
//		Displ = new int[RankSize];
//		DataSizeArray = new int[RankSize];
//		int tmp_x = (Rows / RankSize) * RankSize;
//		int sum = 0;
//		int rowsum = 0;
//		RowSum[0] = 0;
//		for (int i = 0; i < RankSize; i++)
//		{
//			if (tmp_x < Rows)
//			{
//				RowSizeArray[i] = (Rows / RankSize + 1);
//				tmp_x++;
//			}
//			else
//			{
//				RowSizeArray[i] = (Rows / RankSize);
//			}
//			DataSizeArray[i] = RowSizeArray[i] * Cols;
//			Displ[i] = sum;
//			RowSum[i] = rowsum;
//			sum += DataSizeArray[i];
//			rowsum += RowSizeArray[i];
//			DataSize = RowSizeArray[i];
//			if (i != 0)
//				MPI_Send(&DataSize, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
//		}
//		minMPI = (double*)malloc(Rows * sizeof(double));//root
//		DataSize = RowSizeArray[0];						//root
//	}
//	else
//	{
//		MPI_Recv(&DataSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//	}
//
//	// [[ Shared code ]]
//	RecvBuf = (double*)malloc(DataSize*Cols * sizeof(double));
//	//MPI_Scatter(matrix, DataSize*Cols, MPI_DOUBLE, RecvBuf, DataSize*Cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//	MPI_Scatterv(matrix, DataSizeArray, Displ, MPI_DOUBLE, RecvBuf, DataSize*Cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//
//	//std::cout << "RecvBuf[0]:"<<RecvBuf[0]<<" i:" << rank << std::endl;
//	PrintMatrix(RecvBuf, DataSize, Cols);
//
//
//	double* minMPI_datasize = FindMin(RecvBuf, Cols, DataSize); // Array.length = DataSize; Contains mins
//																//MPI_Gather(minMPI_datasize, DataSize, MPI_DOUBLE, minMPI, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//
//	MPI_Gatherv(minMPI_datasize, DataSize, MPI_DOUBLE, minMPI, RowSizeArray, RowSum, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//	// [[ Private code ]]
//	if (rank == 0)
//	{
//		TimeEnd = MPI_Wtime();//only for root
//							  //	PrintMatrix(minMPI, Rows, 1); //Uncomment if we want to show mins;
//		cout << "Time: " << TimeEnd - TimeStart << endl;
//		if (!AreEqual(minSeq, minMPI, Rows))
//		{
//			cout << "Error: minSeq != minMPI" << endl;
//		}
//		else
//		{
//			cout << "Success: minSeq == minMPI" << endl;
//		}
//		free(minMPI);
//		free(matrix);
//	}
//	else
//		free(RecvBuf);
//	free(minSeq);
//	//free(minMPI_datasize);
//	MPI_Finalize();
//}