//#include "mpi.h"
//#include <stdlib.h> //for _countof
//#include "stdio.h"
//
//void main1(int argc, char* argv[])
//{
//	MPI_Init(&argc, &argv);
//	int rank;
//	//MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//	if (rank == 0)
//	{
//		char hello[] = "Hello World!";
//		MPI_Send(hello, _countof(hello), MPI_CHAR, 1, 0, MPI_COMM_WORLD);
//	}
//	else if (rank == 1) 
//	{
//		char hello[13];
//		MPI_Recv(hello, _countof(hello), MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
//		printf("Received: %s", hello);
//	}
//	MPI_Finalize();
//}