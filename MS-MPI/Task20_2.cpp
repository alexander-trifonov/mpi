/*
Task 17: Searching mins in matrix rows;
Restrictions: # process%Cols = 0;
*/
#include "mpi.h"
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <cstdlib>
#include "Matrix.h"
using namespace std;
void ChooseMain(double* &matrix, int Rows_Start, int Cols, int Rows, int Column_Index, double &main, int &index)
{
	main = matrix[Column_Index];
	index = Rows_Start;
	for (int i = Rows_Start; i < Rows; i++)
		if (Column_Index > 0)
		{
			if (matrix[i*Cols + (Column_Index - 1)] == 0)
				if (matrix[i*Cols + Column_Index] > main)
				{
					main = matrix[i*Cols + Column_Index];
					index = i;
				}
		}
		else
			if (matrix[i*Cols + Column_Index] > main)
			{
				main = matrix[i*Cols + Column_Index];
				index = i;
			}
			
}

void GaussForward(double* &matrix, int Rows, int Cols)
{
	double main_value;
	int main_index;


	for (int k = 0; k < Rows; k++)
	{
		ChooseMain(matrix, k, Cols, Rows, k, main_value, main_index);
		for (int i = 0; i < Rows; i++)
		{
			if (i != main_index)
			{
				double tmp = matrix[i*Cols + k];
				for (int j = k; j < Cols; j++)
				{
					matrix[i*Cols + j] = matrix[i*Cols + j] - (tmp * (matrix[main_index*Cols + j] / matrix[main_index*Cols + k]));
				}
			}
			else
				for (int j = k; j < Cols; j++)
				{
					matrix[i*Cols + j] = matrix[i*Cols + j] / matrix[main_index*Cols + k];
				}
		}
	}
}
//void OP_MainIndex(double *, double *, int *, MPI_Datatype *);
void OP_MainIndex(double *in, double *inout, int *len, MPI_Datatype *dptr)
{
	if (double(in[0]) > double(inout[0]))
	{
		inout[0] = in[0]; // value
		inout[1] = in[1]; // row index
		inout[2] = in[2]; // rank
	}
}

MPI_Op MPI_MainIndex;
double recvmain[3] = { -65000.0, 0, 0 }; //global main: value, local_index, rank
void GaussForwardMPI(double* &matrix, int Rows, int Cols, int rank, int Seq)
{
	double main_value;
	int main_index;
	double* main_row = new double[Cols];
	double main[3]; // local main
	for (int k = 0; k < Cols-1; k++)
	{
		ChooseMain(matrix, 0, Cols, Rows, k, main_value, main_index);
		main[0] = main_value;
		main[1] = main_index;
		main[2] = rank;
		if(Seq != 0)
			MPI_Allreduce(main, recvmain, 3, MPI_DOUBLE, MPI_MainIndex, MPI_COMM_WORLD); //Who has main row?
		else
		{
			recvmain[0] = main[0];
			recvmain[1] = main[1];
			recvmain[2] = main[2];
		}
		//cout << "recvmain: " << recvmain[0] << " " << recvmain[1] << " " << recvmain[2] << endl;
		if (rank == recvmain[2]) // recvmain[2] rank has the main row, let's send it to everyvone.
			for (int i = 0; i < Cols; i++)
				main_row[i] = matrix[int(recvmain[1])*Cols + i];
		MPI_Bcast(main_row, Cols, MPI_DOUBLE, recvmain[2], MPI_COMM_WORLD);
		for (int i = 0; i < Rows; i++)
		{
			//cout << "!"<< rank << " == " << recvmain[2] << " " << i << " == " << recvmain[1] << endl;
			if (!((rank == recvmain[2]) && (i == recvmain[1]))) //isn't it main row?
			{//yes, isn't
				if (k > 0)
				{
					if (matrix[i*Cols + (k - 1)] == 0)
					{
						double tmp = matrix[i*Cols + k];
						for (int j = k; j < Cols; j++)
						{
							//cout << matrix[i*Cols + j] << " -> " << matrix[i*Cols + j] - (tmp * (main_row[j] / recvmain[0])) << " = " << matrix[i*Cols + j] << " - (" << tmp << " * (" << main_row[j] << " : " << recvmain[0] << "))\n";
							matrix[i*Cols + j] = matrix[i*Cols + j] - (tmp * (main_row[j] / recvmain[0]));
						}
					}
				}
				else
				{
					double tmp = matrix[i*Cols + k];
					for (int j = k; j < Cols; j++)
					{
					//	cout << matrix[i*Cols + j] << " -> " << matrix[i*Cols + j] - (tmp * (main_row[j] / recvmain[0])) << " = " << matrix[i*Cols + j] << " - (" << tmp << " * (" << main_row[j] << " : " << recvmain[0] << "))\n";
						matrix[i*Cols + j] = matrix[i*Cols + j] - (tmp * (main_row[j] / recvmain[0]));
					}
				}
			}
			else
				for (int j = k; j < Cols; j++)
				{
					matrix[i*Cols + j] = matrix[i*Cols + j] / recvmain[0];
				}
		}
		if(Seq != 0)
			MPI_Barrier(MPI_COMM_WORLD);
		//cout << " CHANGES: " << endl;
		//PrintMatrix(matrix, Rows, Cols);
	}
	delete[] main_row;
}

int DataDistr(int* &Displ, int &RankSize, int* &DataSizeArray, int Rows, int Cols, int DataSize )
{
	int root_datasize;
	Displ = new int[RankSize];
	DataSizeArray = new int[RankSize];
	int tmp_x = (Rows / RankSize) * RankSize;
	int sum = 0;
	for (int i = 0; i < RankSize; i++)
	{
		if (tmp_x < Rows)
		{
			DataSize = (Rows / RankSize + 1);
			tmp_x++;
		}
		else
			DataSize = (Rows / RankSize);
		DataSizeArray[i] = DataSize*Cols;
		Displ[i] = sum;
		sum += DataSizeArray[i];

		if (i != 0)
			MPI_Send(&DataSize, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		else
			root_datasize = DataSize;
	}
	return root_datasize;
}

void main(int argc, char** argv)
{
	//[ Matrix initialization ]
	double* matrix;	// a main matrix, init in rank 0;
	int Rows = atoi(argv[1]);	// # rows in the matrix
	int Cols = Rows + 1;// # cols in the matrix
	
	MPI_Init(&argc, &argv);
	int rank, RankSize;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &RankSize);
	MPI_Op_create((MPI_User_function*)OP_MainIndex, 0, &MPI_MainIndex);
	

	double TimeStart, TimeEnd;
	int DataSize;	// # rows process has
	int* DataSizeArray; // # elems to send each process
	int* Displ;
	double *RecvBuf;// Matrix for FindMin function
	double* matrixSeq;
	double* matrixMPI;


	if (rank == 0)
	{
		//[[ Matrix initialization ]]
		cout << "[[ Generating matrix ]]" << endl;
		cout << "[ Rows: " << Rows << " Cols: " << Cols << " ]" << endl;
		GenerateMatrix(matrix, Rows, Cols, 1, 10);
		PrintMatrix(matrix, Rows, Cols);

		//[[ Sequential ]]
		cout << "[[ Sequential ]]" << endl;
		TimeStart = MPI_Wtime();
		matrixSeq = new double[Cols*Rows];
		CopyInto(matrix, matrixSeq, Rows, Cols);
		GaussForwardMPI(matrixSeq, Rows, Cols, rank, 0); //0 means sequential
		TimeEnd = MPI_Wtime();
		cout << "Time: " << TimeEnd - TimeStart << endl;
		PrintMatrix(matrixSeq, Rows, Cols);
		

		//[[ MPI ]]
		cout << "[[ MPI ]]" << endl;
		TimeStart = MPI_Wtime();
		matrixMPI = new double[Cols*Rows];
		CopyInto(matrix, matrixMPI, Rows, Cols);
		DataSize = DataDistr(Displ, RankSize, DataSizeArray, Rows, Cols, DataSize);
	}
	else
	{
		MPI_Recv(&DataSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	
	// [[ Shared code ]]
	RecvBuf = (double*)malloc(DataSize*Cols * sizeof(double));
	//MPI_Scatter(matrix, DataSize*Cols, MPI_DOUBLE, RecvBuf, DataSize*Cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatterv(matrixMPI, DataSizeArray, Displ, MPI_DOUBLE, RecvBuf, DataSize*Cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	//std::cout << "RecvBuf[0]:"<<RecvBuf[0]<<" i:" << rank << std::endl;
	//PrintMatrix(RecvBuf, DataSize, Cols);

	GaussForwardMPI(RecvBuf, DataSize, Cols, rank, 1);
	MPI_Gatherv(RecvBuf, DataSize*Cols, MPI_DOUBLE, matrixMPI, DataSizeArray, Displ, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if (rank == 0)
	{
		TimeEnd = MPI_Wtime();
		PrintMatrix(matrixMPI, Rows, Cols);
		cout << "Time: " << TimeEnd - TimeStart << endl;
		free(matrix);
		delete[] Displ;
		delete[] DataSizeArray;
		delete[] matrixSeq;
		delete[] matrixMPI;
	}
	free(RecvBuf);
	MPI_Finalize();
}
