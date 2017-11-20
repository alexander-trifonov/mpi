/*
Task 17: Searching mins in matrix rows;
Restrictions: # process%Cols = 0;
*/
#include "mpi.h"
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <cstdlib>
#include <cmath>
#include "Matrix.h"

using namespace std;
void ChooseMain(double* &matrix, int Rows_Start, int Cols, int Rows, int Column_Index, double &main, int &index)
{
	main = -65000.0; //matrix[Column_Index] doesn't work
	index = Rows_Start;
	for (int i = Rows_Start; i < Rows; i++)
		if (Column_Index > 0)
		{
			if (matrix[i*Cols + (Column_Index - 1)] == 0.0) //Don't include already processed rows
			{
				if (matrix[i*Cols + Column_Index] > main)
				{
					main = matrix[i*Cols + Column_Index];
					index = i;
				}
			}
		}
		else
			if (matrix[i*Cols + Column_Index] > main)
			{
				main = matrix[i*Cols + Column_Index];
				index = i;
			}
			
}

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
//Collective function: require to call for each process together
//double* &matrix size may vary in each process
//int RankSize needs only because MPI_Bcast won't work properly
void GaussForwardMPI(double* &matrix, int Rows, int Cols, int rank, int Seq, int RankSize)
{
	double main_value;
	int main_index;
	double* main_row = (double*)malloc(Cols*sizeof(double));
	//recvmain[0] global max value, everyvone will see it
	//recvmain[1] global row index relative recvmain[2] rank
	//recvmain[2] global process rank, has max value and row index.
	//global means every process has same data
	double recvmain[3] = { -65000.0, 0, 0 };
	//main[0] local max value in current process
	//main[1] local row index with max value in current process
	//main[2] local current process rank
	//local means each process has own data
	double main[3];
	for (int k = 0; k < Cols-1; k++)
	{
		recvmain[0] = -65000.0; //Refreshing global max value for each iteration (else we will compare new value with old one). 
		if (Seq != 0)
			MPI_Barrier(MPI_COMM_WORLD);
		ChooseMain(matrix, 0, Cols, Rows, k, main_value, main_index);
		main[0] = main_value;
		main[1] = main_index;
		main[2] = rank;
		if (Seq != 0)
		{
			MPI_Allreduce(main, recvmain, 3, MPI_DOUBLE, MPI_MainIndex, MPI_COMM_WORLD); //Who has the main row?
		}
		else
		{
			recvmain[0] = main[0];
			recvmain[1] = main[1];
			recvmain[2] = main[2];
		}
		//cout << "recvmain: " << recvmain[0] << " " << recvmain[1] << " " << recvmain[2] << endl;
		if (rank == recvmain[2]) // recvmain[2] rank has the main row -> fill main row in main_row
			for (int i = 0; i < Cols; i++)
			{
				main_row[i] = (double)matrix[int(recvmain[1])*Cols + i];
			}
		if (Seq != 0)
		{
			/*
				Because for some reason MPI_Bcast doesn't work properly I had to write my Bcast.
				MPI_Bcast(main_row, Cols, MPI_DOUBLE, recvmain[2], MPI_COMM_WORLD) doesn't work:
				Sometimes it sends data with a big delay, so processes have old data each iteration
				and sometimes it sends fine.
			*/
			MPI_Barrier(MPI_COMM_WORLD);
			if (rank == recvmain[2])
			{
				for (int i = 0; i < RankSize; i++)
					if (i != rank)
						MPI_Send(main_row, Cols, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			}
			else
				MPI_Recv(main_row, Cols, MPI_DOUBLE, recvmain[2], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//std::cout << "Rank: " << rank << " received main row: ";
		//PrintMatrixRank(main_row, 1, Cols, rank);
		for (int i = 0; i < Rows; i++)
		{
			bool fine = true;
			if (!((rank == recvmain[2]) && (i == recvmain[1]))) //If this row isn't main, then
			{
				if (k > 0) //if this column isn't first then
				{
					if (matrix[i*Cols + (k - 1)] == 0) //if this row wasn't main previously. First check.
					{
						for (int t = 0; t < k; t++) //Additional check for ALL previous zeros
						{
							if (matrix[i*Cols + t] != 0)
							{
								fine = false;
								break;
							}
						}
						/*bool AllZeros = false;
						for (int t = 0; t < Cols - 1; t++ )
						{
							if (matrix[i*Cols + t] != 0)
							{
								AllZeros = false;
								break;
							}
						}
						if (AllZeros)
						{
							cout << "ROW FULL OF ZEROS" << endl;
							MPI_Finalize();
						}*/
						if (fine)
						{
							double tmp = matrix[i*Cols + k];
							for (int j = k; j < Cols; j++)
							{
								//cout << matrix[i*Cols + j] << " -> " << matrix[i*Cols + j] - (tmp * (main_row[j] / recvmain[0])) << " = " << matrix[i*Cols + j] << " - (" << tmp << " * (" << main_row[j] << " : " << recvmain[0] << "))\n";
								matrix[i*Cols + j] = matrix[i*Cols + j] - (tmp * (main_row[j] / recvmain[0]));
							}
						}							
					}
				}
				else
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
				for (int j = k; j < Cols; j++)
					matrix[i*Cols + j] = matrix[i*Cols + j] / recvmain[0];
			}
		}
		if (Seq != 0)
			MPI_Barrier(MPI_COMM_WORLD);
		//cout << " CHANGES: " << endl;
		//PrintMatrixRank(matrix, Rows, Cols,rank);
	}
	free(main_row);
}
double Determinant(double* &matrix, int Rows, int Cols, int RowStart)
{
	double result;
	if (Rows - RowStart == 2)
	{
		return matrix[RowStart*Cols] * matrix[(RowStart + 1)*Cols + 1] - matrix[RowStart*Cols + 1] * matrix[(RowStart + 1)*Cols];
	}
	for (int i = RowStart; i < Rows; i++)
	{
		return matrix[i*Cols] * (pow(-1, i))*Determinant(matrix, Rows, Cols, RowStart + 1);
	}
}
void GaussBackwardMPI(double* &matrix, double* &x, int Rows, int Cols, int rank, int Seq, int RankSize)
{
	//cout << rank << " rank has " << Rows << " rows and " << Cols << " cols" << endl;
	for (int j = Cols - 2; j >= 0; j--)
	{
		x[j] = -9999;
	}
	if(Seq!=0)
		MPI_Barrier(MPI_COMM_WORLD);
	//cout << rank << " rank passed barrier" << endl;
	for (int j = Cols - 2; j >= 0; j--)
	{
		bool main = false;
		for (int i = 0; i < Rows; i++)
		{
			if (matrix[i*Cols + j] == 1.0)
			{
				//cout << "matrix[" << int(i*Cols + j) << "]==1.0" << endl;
				//cout << "i:" << i << " j:" << j << endl;
				bool fine = true;
				//cout << rank << " rank entered ==1.0" << endl;
				if (j != 0)
				{
					for (int m = 0; m < j; m++)
					{
						//cout << matrix[i*Cols + m] << " m = "<<m<< endl;
						if (matrix[i*Cols + m] != 0)
						{
							//cout << "fine = false on m =" << m << " matrix[i*Cols + m]=" << matrix[i*Cols + m] << endl;
							fine = false;
							break;
						}
					}
				}
				//cout << rank << " passed fine module" << endl;
				if (fine)
				{
					//cout << rank << " main = true" << endl;
					main = true;
					x[j] = matrix[i*Cols + Cols - 1];
					//cout << "b: x[" << j << "] = " << x[j] << endl;
					for (int k = j + 1; k < Cols - 1; k++)
					{
						//cout << "x[k:" << k << "]=" << x[k] << "\tmatrix[i,k]=" << matrix[i*Cols + k] << endl;
						x[j] = x[j] - matrix[i*Cols + k] * x[k];
					}
					//cout << "x[" << j << "] = " << x[j] << endl;
					if (Seq != 0)
						for (int k = 0; k < RankSize; k++)
							if (k != rank)
							{
								MPI_Send(&x[j], 1, MPI_DOUBLE, k, 0, MPI_COMM_WORLD);
								//cout << "Sended to " << k << endl;
							}
				}
			}
		}
		//cout << rank << " rank exited loop" << endl;
		if (Seq != 0)
			if ((!main) && (x[j] == -9999))
			{
				//cout << rank << " rank is waiting" << endl;
				MPI_Recv(&x[j], 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				//cout << rank << " is no longer waiting" << endl;
				//cout << rank << " rank received x[" << j << "]=" << x[j] << endl;
			}
	}

	/*cout << "Backward ended: " << endl;
	for (int j = Cols - 2; j >= 0; j--)
	{
		cout << "x[" << j << "] = " << x[j] << endl;
	}*/

}
//Determine displ[], sendcounts[]
//Send DataSize each process
int DataDistr(int* &Displ, int* &sendcounts, int Rows, int Cols, int DataSize, int &RankSize )
{
	int root_datasize;
	Displ = new int[RankSize];
	sendcounts = new int[RankSize];
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
		sendcounts[i] = DataSize*Cols;
		Displ[i] = sum;
		sum += sendcounts[i];

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
	int* DataSizeArray; //root: # elems to send each process
	int* Displ;
	double *RecvBuf;
	double* matrixSeq;
	double* matrixMPI;
	double* xSeq;
	double* xMPI;

	if (rank == 0)
	{
		//[[ Matrix initialization ]]
		cout << "[[ Generating matrix ]]" << endl;
		cout << "[ Rows: " << Rows << " Cols: " << Cols << " ]" << endl;
		GenerateMatrix(matrix, Rows, Cols, 1, pow(10,RankSize));
		//PrintMatrix(matrix, Rows, Cols);

		//[[ Sequential ]]
		cout << "[[ Sequential ]]" << endl;
		TimeStart = MPI_Wtime();
		matrixSeq = new double[Cols*Rows];
		CopyInto(matrix, matrixSeq, Rows, Cols);
		GaussForwardMPI(matrixSeq, Rows, Cols, rank, 0, RankSize); //0 means sequential
		//PrintMatrix(matrixSeq, Rows, Cols);
		xSeq = new double[Cols];
		cout << "Determinant = " << Determinant(matrix, Rows, Cols, 0) << endl;
		GaussBackwardMPI(matrixSeq, xSeq, Rows, Cols, rank, 0, RankSize);
		TimeEnd = MPI_Wtime();
		cout << "Time: " << TimeEnd - TimeStart << endl;
		//PrintMatrix(xSeq, 1, Cols - 1);
		

		//[[ MPI ]]
		cout << "[[ MPI ]]" << endl;
		//PrintMatrix(matrix, Rows, Cols);
		TimeStart = MPI_Wtime();
		matrixMPI = new double[Cols*Rows];
		CopyInto(matrix, matrixMPI, Rows, Cols);
		DataSize = DataDistr(Displ, DataSizeArray, Rows, Cols, DataSize, RankSize);
	}
	else
	{
		MPI_Recv(&DataSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// [[ Shared code ]]
	MPI_Barrier(MPI_COMM_WORLD);
	RecvBuf = (double*)malloc(DataSize*Cols * sizeof(double));
	MPI_Scatterv(matrixMPI, DataSizeArray, Displ, MPI_DOUBLE, RecvBuf, DataSize*Cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	GaussForwardMPI(RecvBuf, DataSize, Cols, rank, 1, RankSize);
	MPI_Gatherv(RecvBuf, DataSize*Cols, MPI_DOUBLE, matrixMPI, DataSizeArray, Displ, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if (rank == 0)
	{
		//PrintMatrix(matrixMPI, Rows, Cols);
	}
	cout << "GausForward" << endl;
	xMPI = new double[Cols];
	GaussBackwardMPI(RecvBuf, xMPI, DataSize, Cols, rank, 1, RankSize);
	
	
	if (rank == 0)
	{
		TimeEnd = MPI_Wtime();
		cout << "MPI:" << endl;
		//PrintMatrix(matrixMPI, Rows, Cols);
		cout << "Time: " << TimeEnd - TimeStart << endl;
	//	PrintMatrix(xSeq, 1, Cols - 1);
	//	PrintMatrix(xMPI, 1, Cols - 1);
		if (AreEqual(xSeq, xMPI, 1, Cols-1))
			cout << "xSeq == xMPI" << endl;
		else
			cout << "xSeq != xMPI" << endl;
		cout << "Determinant = " << Determinant(matrix, Rows, Cols, 0) << endl;
		free(matrix);
		delete[] Displ;
		delete[] DataSizeArray;
		delete[] matrixSeq;
		delete[] matrixMPI;
		delete[] xSeq;
	}
	delete[] xMPI;
	free(RecvBuf);
	MPI_Finalize();
}
