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
bool IsProccesedRow(double* &matrix, int Cols, int Column_Index, int Row_Index)
{
	if(Column_Index>0)
		for (int j = Column_Index - 1; j >= 0; j--)
			if (matrix[Row_Index*Cols + j] != 0)
			{
				return true;
			}
	return false;
}
void ChooseMain(double* &matrix, int Rows_Start, int Cols, int Rows, int Column_Index, double &main, int &index)
{
	main = -DBL_MAX;
	for (int i = Rows_Start; i < Rows; i++)
	{
		if (!IsProccesedRow(matrix, Cols, Column_Index, i))//Если у процесса ВСЕ строки обработаны, то сюда он не зайдет. Соответственно, его main будет иметь значение предыдуще итерации. Соответственно, нужно избавляться от данных предущих итераций. Соответственно, необходимо обнулить main. 
		{
			if (matrix[i*Cols + Column_Index] > main)
			{
				main = matrix[i*Cols + Column_Index];
				index = i;
			}
		}
	}
}

//Collective function: require to call for each process together
//double* &matrix size may vary in each process
//int RankSize needs only because MPI_Bcast won't work properly
void GaussForwardMPI(double* &matrix, int Rows, int Cols, int rank, int Seq, int RankSize)
{
	double main_value;
	int main_index;
	double* main_row;
	double *recvmain;//root only
	if (rank == 0)
		recvmain = new double[3 * RankSize];
	//main[0] local max value in current process
	//main[1] local row index with max value in current process
	//main[2] local current process rank
	double* main = new double[3];
	for (int k = 0; k < Cols - 1; k++)
	{
		ChooseMain(matrix, 0, Cols, Rows, k, main_value, main_index);
		main[0] = main_value;
		main[1] = main_index;
		main[2] = rank;

		//Root will find what row is main
		MPI_Gather(main, 3, MPI_DOUBLE, recvmain, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if (rank == 0)
		{
			main[0] = recvmain[0];
			for (int i = 0; i < 3 * RankSize; i += 3)
			{
				if (recvmain[i] > main[0])
				{
					main[0] = recvmain[i];
					main[1] = recvmain[i + 1];
					main[2] = recvmain[i + 2];
				}
			}
		}
		MPI_Bcast(main, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
		//If current rank has main row, then normalize it
		if (rank == main[2])
		{
			main_row = &(matrix[int(main[1]) * Cols]);
			for (int j = k; j < Cols; j++)
				matrix[int(main[1])*Cols + j] = matrix[int(main[1])*Cols + j] / main[0];
		}
		else
			main_row = (double*)malloc(Cols * sizeof(double));
		MPI_Bcast(main_row, Cols, MPI_DOUBLE, main[2], MPI_COMM_WORLD);
	
		if (rank != main[2])
			for (int i = 0; i < Rows; i++)
			{
				if (!IsProccesedRow(matrix,Cols,k,i))
				{
					double tmp = matrix[i*Cols + k];
					for (int j = k; j < Cols; j++)
					{
						matrix[i*Cols + j] = matrix[i*Cols + j] - (tmp * (main_row[j]));
					}
				}
			}
		else
			for (int i = 0; i < Rows; i++)
			{
				if ((i != main[1]) && (!IsProccesedRow(matrix, Cols, k, i)))
				{
					double tmp = matrix[i*Cols + k];
					for (int j = k; j < Cols; j++)
					{
						matrix[i*Cols + j] = matrix[i*Cols + j] - (tmp * (main_row[j]));
					}
				}
			}
		if(rank != main[2])
			free(main_row);
	}
	if (rank == 0)
		delete[] recvmain;
	delete[] main;
}
//Collective function: require to call for each process together
//double* &matrix size may vary in each process
//int RankSize needs only because MPI_Bcast won't work properly
void GaussForward(double* &matrix, int Rows, int Cols)
{
	double main_value;
	int main_index;
	double* main_row;
	//main[0] local max value in current process
	//main[1] local row index with max value in current process
	//main[2] local current process rank
	double* main = new double[3];
	for (int k = 0; k < Cols - 1; k++)
	{
		ChooseMain(matrix, 0, Cols, Rows, k, main_value, main_index);
		main[0] = main_value;
		main[1] = main_index;

		//Root will find what row is main
		
		//If current rank has main row, then normalize it
		
		main_row = &(matrix[int(main[1]) * Cols]);
		for (int j = k; j < Cols; j++)
			matrix[int(main[1])*Cols + j] = matrix[int(main[1])*Cols + j] / main[0];
		

		for (int i = 0; i < Rows; i++)
		{
			if ((i != main[1]) && (!IsProccesedRow(matrix, Cols, k, i)))
			{
				double tmp = matrix[i*Cols + k];
				for (int j = k; j < Cols; j++)
				{
					matrix[i*Cols + j] = matrix[i*Cols + j] - (tmp * (main_row[j]));
				}
			}
		}
	}
	delete[] main;
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
void GaussBackwardMPI(double* &matrix, double* &x, int Rows, int Cols, int rank, int RankSize)
{
	for (int j = Cols - 2; j >= 0; j--)
	{
		MPI_Barrier(MPI_COMM_WORLD);//Иначе, процесс, который получил свой x[j], может отправить новый x[j+1] процессу, до которого еще не дошел цикл для отправки x[j].
		bool main = false;
		for (int i = 0; i < Rows; i++)
		{
			if ((matrix[i*Cols + j] == 1.0)&&(!IsProccesedRow(matrix,Cols,j,i)))
			{
				main = true;
				x[j] = matrix[i*Cols + Cols - 1];
				for (int k = j + 1; k < Cols - 1; k++)
				{
					//cout << "x[" << j << "]= " << x[j] << " - " << matrix[i*Cols + k] << " * " << x[k] << endl;
					//cout << "x[k:" << k << ":]=" << x[k] << endl;
					x[j] = x[j] - matrix[i*Cols + k] * x[k];
				}
				//cout << "x[" << j << "]=" << x[j] << endl;
				for (int k = 0; k < RankSize; k++)
					if (k != rank)
					{
						MPI_Send(&x[j], 1, MPI_DOUBLE, k, 0, MPI_COMM_WORLD);
						//cout << rank << " rank. Sending x[" << j << "]=" << x[j]<<" to "<<k << endl;
					}
			}
		}
		if ((!main))
		{
			MPI_Recv(&x[j], 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			//cout << rank << " rank. Received x[" << j << "]=" << x[j] << endl;
		}
	}
}
void GaussBackward(double* &matrix, double* &x, int Rows, int Cols)
{
	for (int j = Cols - 2; j >= 0; j--)
	{
		for (int i = 0; i < Rows; i++)
		{
			if ((matrix[i*Cols + j] == 1.0) && (!IsProccesedRow(matrix, Cols, j, i)))
			{
				x[j] = matrix[i*Cols + Cols - 1];
				for (int k = j + 1; k < Cols - 1; k++)
				{
					//cout << "x[" << j << "]= " << x[j] << " - " << matrix[i*Cols + k] << " * " << x[k] << endl;
					x[j] = x[j] - matrix[i*Cols + k] * x[k];
					//cout << "x[j] = x[j] - " << matrix[i*Cols + k] * x[k] << endl;
				}
				//cout << "x[" << j << "]=" << x[j] << endl;
			}
		}
	}
}
//Determine displ[], sendcounts[]
//Send DataSize each process
int DataDistr(int* &Displ, int* &sendcounts, int Rows, int Cols, int DataSize, int &RankSize )
{
	int root_datasize;
	Displ = new int[RankSize];
	sendcounts = new int[RankSize];
	int sum = 0;
	
	int mod = Rows % RankSize;
	for (int i = 0; i < RankSize; i++)
	{
		DataSize = Rows / RankSize + (mod > 0 ? 1 : 0);
		mod--;
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

	double TimeStart, TimeEnd;
	int DataSize;	// # rows process has
	int* DataSizeArray; //root: # elems to send each process
	int* Displ;
	double *RecvBuf;
	double* matrixSeq;
	double* matrixMPI;
	double* xSeq;
	double* xMPI;
	double TimeSeq;
	double TimeMPI;
	if (rank == 0)
	{
		//[[ Matrix initialization ]]
		cout << "[[ Generating matrix ]]" << endl;
		cout << "[ Rows: " << Rows << " Cols: " << Cols << " ]" << endl;
		GenerateMatrix(matrix, Rows, Cols, 1, Rows*10);
		cout << "Determinant = " << Determinant(matrix, Rows, Cols, 0) << endl;
		matrixSeq = new double[Cols*Rows];
		CopyInto(matrix, matrixSeq, Rows, Cols);
		matrixMPI = new double[Cols*Rows];
		CopyInto(matrix, matrixMPI, Rows, Cols);
		
		//PrintMatrix(matrix, Rows, Cols);

		//[[ Sequential ]]
		cout << "[[ Sequential ]]" << endl;
		TimeStart = MPI_Wtime();
		GaussForward(matrixSeq, Rows, Cols);
		//GaussForwardMPI(matrixSeq, Rows, Cols, rank, 0, RankSize); //0 means sequential
		//PrintMatrix(matrixSeq, Rows, Cols);
		xSeq = new double[Cols-1];
		//GaussBackwardMPI(matrixSeq, xSeq, Rows, Cols, rank, 0, RankSize);
		GaussBackward(matrixSeq, xSeq, Rows, Cols);
		TimeEnd = MPI_Wtime();
		TimeSeq = TimeEnd - TimeStart;
		cout << "Time: " << TimeSeq << endl;
		//PrintMatrix(xSeq, 1, Cols - 1);
		

		//[[ MPI ]]
		cout << "[[ MPI ]]" << endl;
		//PrintMatrix(matrix, Rows, Cols);
		TimeStart = MPI_Wtime();
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
		cout << "MPI:" << endl;
		//PrintMatrix(matrixMPI, Rows, Cols);
	}
	xMPI = new double[Cols-1];
	
	GaussBackwardMPI(RecvBuf, xMPI, DataSize, Cols, rank, RankSize);

	if (rank == 0)
	{
		TimeEnd = MPI_Wtime();
		cout << "MPI:" << endl;
		//PrintMatrix(matrixMPI, Rows, Cols);
		TimeMPI = TimeEnd - TimeStart;
		cout << "Time: " << TimeMPI << endl;
		cout << "Acceleration MPI/Seq = " << TimeMPI / TimeSeq << endl;
		//PrintMatrix(xSeq, 1, Cols - 1);
		//PrintMatrix(xMPI, 1, Cols - 1);
		if (AreEqual(xSeq, xMPI, 1, Cols-1))
			cout << "xSeq == xMPI" << endl;
		else
			cout << "xSeq != xMPI" << endl;
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
