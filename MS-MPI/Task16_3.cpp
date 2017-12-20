//Покрытие отрезков/Segments covering

//Возьмём все точки - концы отрезков(как левые, так и правые) и отсортируем их.При этом для каждой точки сохраним вместе с
//ней номер отрезка, а также то, каким концом его она является(левым или правым).Кроме того, отсортируем точки таким
//образом, что, если есть несколько точек с одной координатой, то сначала будут идти левые концы, и только потом - правые.
//Заведём стек, в котором будут храниться номера отрезков, рассматриваемых в данный момент; изначально стек пуст.
//Будем двигаться по точкам в отсортированном порядке.Если текущая точка - левый конец, то просто добавляем номер её 
//отрезка в стек.Если же она является правым концом, то проверяем, не был ли покрыт этот отрезок(для этого можно просто 
//завести массив булевых переменных).Если он уже был покрыт, то ничего не делаем и переходим к следующей точке(забегая
//вперёд, мы утверждаем, что в этом случае в стеке текущего отрезка уже нет).Если же он ещё не был покрыт, то мы
//добавляем текущую точку в ответ, и теперь мы хотим отметить для всех текущих отрезков, что они становятся покрытыми.
//Поскольку в стеке как раз хранятся номера непокрытых ещё отрезков, то будем доставать из стека по одному отрезку и 
//отмечать, что он уже покрыт, пока стек полностью не опустеет.По окончании работы алгоритма все отрезки будут покрыты,
//и притом наименьшим числом точек(повторимся, здесь важно требование, что при равенстве координат сначала идут левые
//концы, и только затем правые).
//
//Таким образом, весь алгоритм выполняется за O(N), не считая сортировки точек, а итоговая сложность алгоритма как 
//раз равна O(N log N).
//http://e-maxx.ru/algo/covering_segments

#include "mpi.h"
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <stack>
#include <time.h>
#include <cstdlib>
#include <algorithm>

using namespace std;
struct Point
{
	int x;
	int segment; //the point belongs to segment № %int segment%
	//0 for left segment end
	//1 for right segment end
	int side;

	Point() {};
	Point(int _x, int _segment, int _side)
	{
		x = _x;
		segment = _segment;
		side = _side;
	}
	Point(Point* _p)
	{
		x = _p->x;
		segment = _p->segment;
		side = _p->side;
	}
	Point& operator=(Point* _p)
	{
		x = _p->x;
		segment = _p->segment;
		side = _p->side;
		return *this;
	}
	void Print(bool full_information)
	{
		cout << x;
		if (full_information)
		{
			cout << "\t" << segment;
			cout << "\t" << (side == 0 ? "left" : "right");
		}
		cout << endl;
	}
};
//Проблема в том, что необходимо будет в дальнейшем сортировать по точкам, чтобы левые шли за правыми.
//Поэтому, целесообразнее будет отправлять процессам не сегменты, а сразу точки.
//Ну а потом каждый процесс просто отсортирует у себя.
//Иначе, если отправлять сегментами, то нужно будет доставать эти точки, перепаковывать куда-нибудь и там с ними работать.
struct Segment //Однако, сегменты отлично пригодятся для генерации множества точек.
{
	pair<Point, Point> points;

	Segment() {};
	Segment(Point left, Point right)//Points have to have identical Point.segment and different Point.side fields.
	{
		points.first = left;
		points.second = right;
	}
};
Segment* Generate(int size, int min, int max)
{
	srand(time(NULL));
	Segment* line = new Segment[size];
	for (int i = 0; i < size; i++)
	{
		Point left(int((rand() % (max - min)) + min), i, 0);
		Point right(int(left.x + rand() % (max - min) + min + 1), i, 1);
		line[i].points.first = left;
		line[i].points.second = right;
	}
	return line;
}

int FindByValue(int x, vector<int> &table)
{
	for (int i = 0; i < table.size(); i++)
		if (table[i] == x) return i;
	return -1;
}
int FindByValue(int x, vector<Point> &table)
{
	for (int i = 0; i < table.size(); i++)
		if (table[i].x == x) return i;
	return -1;
}
vector<int> Cover(Point* &line, int segments_count)
{
	vector<int> result;
	vector<int> marked;//Contains all № of marked segments;
	stack<int> unmarked;//Gives faster marking alghorithm
	for (int i = 0; i < segments_count*2; i++)
	{
		if (line[i].side == 0)//If point is a left segment end
			unmarked.push(line[i].segment);
		else//If point is a right segment end
		{
			bool IsMarked = false;//Was this segment marked before?
			for (int j = 0; j < marked.size(); j++)
			{
				if (marked[j] == line[i].segment)
				{
					IsMarked = true;
					break;
				}
			}
			if (!IsMarked)
			{
				if (FindByValue(line[i].x, result) != -1) continue;
				result.push_back(line[i].x);//Add point to the solution
				for (int j = 0; j < unmarked.size(); j++)
				{
					marked.push_back(unmarked.top());//And mark the segment
					unmarked.pop();
				}
			}
		}
	}

	return result;
}

int comparePoint(const void* a, const void* b)
{
	if (((Point*)(a))->x <	((Point*)(b))->x) return -1;
	if (((Point*)(a))->x >	((Point*)(b))->x) return 1;
	if (((Point*)(a))->x == ((Point*)(b))->x) return 0;
}
int compareInt(const void* a, const void* b)
{
	if (*(int*)a	<	*(int*)b) return -1;
	if (*(int*)a	>	*(int*)b) return 1;
	if (*(int*)a	==	*(int*)b) return 0;
}
//size - how many segments will generate
//total_segments - stores 
int DataDistr(Point* &line, vector<Point>& global_line, int* &sendcounts, int* &Displ, int &size, int &RankSize, int &max, int &min)
{
	int root_SegmentsCount;
	int SegmentsCount;
	int segment_id = 0;
	int sum = 0;
	int mod = (size) % RankSize;
	Displ = new int[RankSize];
	sendcounts = new int[RankSize];//amount of points

	//Calculate root's number of segmetns:
	SegmentsCount = (size) / RankSize + (mod > 0 ? 1 : 0);
	sendcounts[0] = SegmentsCount*2*3;
	Displ[0] = sum;
	mod--;
	root_SegmentsCount = SegmentsCount;//buffer
	line = new Point[SegmentsCount * 2];
	for (int i = 0; i < SegmentsCount * 2 - 1; i = i + 2)
	{
		Point left(int((rand() % (max - min)) + min), segment_id, 0);
		Point right(int(left.x + rand() % (max - min) + min + 1), segment_id, 1);
		global_line.push_back(left);
		global_line.push_back(right);
		line[i] = left;
		line[i + 1] = right;
		segment_id++;
	}
	sum += sendcounts[0];

	//Calculate other number of segments:
	for (int i = 1; i < RankSize; i++)
	{
		SegmentsCount = (size) / RankSize + (mod > 0 ? 1 : 0);
		sendcounts[i] = SegmentsCount*2*3;
		Displ[i] = sum;
		mod--;
		MPI_Send(&SegmentsCount, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		for (int j = 0; j < SegmentsCount; j++)
		{
			Point left(int((rand() % (max - min)) + min), segment_id, 0);
			Point right(int(left.x + rand() % (max - min) + min + 1), segment_id, 1);
			global_line.push_back(left);
			global_line.push_back(right);
			MPI_Send(&left, 3, MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Send(&right, 3, MPI_INT, i, 0, MPI_COMM_WORLD);
			segment_id++;
		}
		sum += sendcounts[i];
	}
	return root_SegmentsCount;
}

void main(int argc, char** argv)
{
	int rank;
	int RankSize;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &RankSize);
	int* sendcounts;
	int* Displ;
	Point* line;
	vector<Point> global_line;//for sync and duplicate deleting
	int SegmentsCount;
	int size;
	int* gather_solution;
	vector<int> Final_Solution;
	if (rank == 0)
	{
		int max, min;
		int root_database;
		max = 10;
		min = 0;
		size = atoi(argv[1]);//Segments count
		//global_line.resize(size * 2);
		srand(time(NULL));

		SegmentsCount = DataDistr(line, global_line, sendcounts, Displ, size, RankSize, max, min);
		//int segments_count = 0;//Two points have to belong to same segment. This var is for generating method below;
		//int mod = (size) % RankSize;
		////Calculate root's number of segmetns:
		//SegmentsCount = (size) / RankSize + (mod > 0 ? 1 : 0);
		//mod--;
		//root_database = SegmentsCount;//buffer
		//line = new Point[SegmentsCount*2];
		//for (int i = 0; i < SegmentsCount*2-1; i = i+2)
		//{
		//	Point left(int((rand() % (max - min)) + min), segments_count, 0);
		//	Point right(int(left.x + rand() % (max - min) + min + 1), segments_count, 1);
		//	line[i] = left;
		//	line[i + 1] = right;
		//	segments_count++;
		//}
		////Calculate other number of segments:
		//for (int i = 1; i < RankSize; i++)
		//{
		//	SegmentsCount = (size) / RankSize + (mod > 0 ? 1 : 0);
		//	mod--;
		//	MPI_Send(&SegmentsCount, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		//	for (int j = 0; j < SegmentsCount; j++)
		//	{
		//		Point left(int((rand() % (max - min)) + min), segments_count, 0);
		//		Point right(int(left.x + rand() % (max - min) + min + 1), segments_count, 1);
		//		MPI_Send(&left, 3, MPI_INT, i, 0, MPI_COMM_WORLD);
		//		MPI_Send(&right, 3, MPI_INT, i, 0, MPI_COMM_WORLD);
		//		segments_count++;
		//	}
		//}
		//SegmentsCount = root_database;
	}
	else
	{
		MPI_Recv(&SegmentsCount, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		line = new Point[SegmentsCount*2];
		for (int i = 0; i < SegmentsCount*2; i++)
		{
			MPI_Recv(&line[i], 3, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			//MPI_Recv(&line[i + 1], 3, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}
	qsort(line, SegmentsCount * 2, sizeof(Point), comparePoint);
	vector<int> solution = Cover(line, SegmentsCount);
	int solution_size = solution.size();
	int recv_buff;
	MPI_Reduce(&solution_size, &recv_buff, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if (rank == 0)
	{
		for (int i = 0; i < global_line.size(); i++)
		{
			global_line[i].Print(true);
		}
		cout << "Total points before deleting: " << recv_buff << endl;
		gather_solution = new int[recv_buff];
		for (int i = 0; i < solution.size(); i++)
			gather_solution[i] = solution[i];
		for (int i = solution.size(); i < recv_buff; i++)
			MPI_Recv(&gather_solution[i], 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	else
	{
		for (int i = 0; i < solution.size(); i++)
			MPI_Send(&solution[i], 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}

	if (rank == 0)
	{
		for (int i = 0; i < recv_buff; i++)
		{
			cout << gather_solution[i] << " ";
		}
		cout << endl;

		vector<int> segments_solutions(size);
		for (int i = 0; i < size; i++)
			segments_solutions[i] = 0;

		qsort(&global_line[0], global_line.size(), sizeof(Point), comparePoint);
		qsort(gather_solution, recv_buff, sizeof(int), compareInt);
		for (int i = 0; i < recv_buff; i++)
		{
			bool WorthIt = false;
			int j = FindByValue(gather_solution[i], global_line);
			while (global_line[j].x > 0)
			{
				segments_solutions[global_line[j].segment]++;
				if (segments_solutions[global_line[j].segment] == 1)
					WorthIt = true;
				j--;
			}
			if (WorthIt)
			{
				Final_Solution.push_back(gather_solution[i]);
			}
		}
		cout << "After deliting: " << Final_Solution.size() << endl;
		for (int i = 0; i < Final_Solution.size(); i++)
		{
			cout << Final_Solution[i] << " ";
		}
		cout << endl;
	}
	
	
	delete[] line;
	//Segment* line_raw;
	//Point* line_sorted;

	//int* Displ;
	//int* Sendcounts;
	//int DataSize;
	//if (rank == 0)
	//{
	//	int size = 3;
	//	line_raw = Generate(size, 0, 10);
	//	DataSize = DataDistr(Displ, Sendcounts, size * 2, RankSize);
	//}
	//else
	//{
	//	MPI_Recv(&DataSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	//	line_raw = new Segment[DataSize];
	//	for (int i = 0; i < DataSize; i++)
	//	{
	//		MPI_Recv(&line_raw[i].points.first, 3, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	//		MPI_Recv(&line_raw[i].points.second, 3, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	//	}
	//}





	/*Point a;
	if (rank == 0)
	{
		a.x = 1; a.segment = 2; a.side = 0;
		MPI_Send(&a, 3, MPI_INT, 1, 0, MPI_COMM_WORLD);
	}
	else
	{
		MPI_Recv(&a, 3, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		cout << a.x << " " << a.segment << " " << a.side << endl;
	}*/

	/*delete[] line;*/
	MPI_Finalize();
}