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

int FindByValue(int x, vector<int> &table)
{
	for (int i = 0; i < table.size(); i++)
		if (table[i] == x) return i;
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
				while(unmarked.size() != 0)
				{
					marked.push_back(unmarked.top());//And mark the segment
					unmarked.pop();
				}
			}
		}
	}

	return result;
}
//For qsort
int comparePoint(const void* a, const void* b)
{
	if (((Point*)(a))->x <	((Point*)(b))->x) return -1;
	if (((Point*)(a))->x >	((Point*)(b))->x) return 1;
	if (((Point*)(a))->x == ((Point*)(b))->x)
	{
		if (((Point*)(a))->side < ((Point*)(b))->side) return -1;
		if (((Point*)(a))->side > ((Point*)(b))->side) return  1;
		if (((Point*)(a))->side == ((Point*)(b))->side)
		{
			if (((Point*)(a))->segment < ((Point*)(b))->segment) return -1;
			if (((Point*)(a))->segment > ((Point*)(b))->segment) return  1;
			if (((Point*)(a))->segment == ((Point*)(b))->segment) return 0;
		}
	}
}
int compareInt(const void* a, const void* b)
{
	if (*(int*)a	<	*(int*)b) return -1;
	if (*(int*)a	>	*(int*)b) return 1;
	if (*(int*)a	==	*(int*)b) return 0;
}
int comparePointBySegments(const void* a, const void* b)
{
	if (((Point*)(a))->segment <	((Point*)(b))->segment) return -1;
	if (((Point*)(a))->segment >	((Point*)(b))->segment) return 1;
	if (((Point*)(a))->segment ==	((Point*)(b))->segment) return 0;
}
//size - how many segments will generate
int DataDistr(Point* &line, vector<Point>& global_line, int &size, int &RankSize, int &max, int &min)
{
	int root_SegmentsCount;
	int SegmentsCount;
	int segment_id = 0;
	int mod = (size) % RankSize;

	//Calculate root's number of segmetns:
	SegmentsCount = (size) / RankSize + (mod > 0 ? 1 : 0);
	mod--;
	root_SegmentsCount = SegmentsCount;//buffer
	line = new Point[SegmentsCount * 2];
	for (int i = 0; i < SegmentsCount * 2 - 1; i = i + 2)
	{
		line[i] = global_line[segment_id];
		line[i + 1] = global_line[segment_id+1];
		segment_id += 2;
	}

	//Calculate other number of segments:
	for (int i = 1; i < RankSize; i++)
	{
		SegmentsCount = (size) / RankSize + (mod > 0 ? 1 : 0);
		mod--;
		MPI_Send(&SegmentsCount, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		for (int j = 0; j < SegmentsCount; j++)
		{
			MPI_Send(&global_line[segment_id], 3, MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Send(&global_line[segment_id+1], 3, MPI_INT, i, 0, MPI_COMM_WORLD);
			segment_id += 2;
		}
	}
	return root_SegmentsCount;
}

vector<Point> GenerateGlobalLine(int size, int max, int min)
{
	srand(time(NULL));
	vector<Point> global_line;
	int segment_id = 0;
	for (int j = 0; j < size; j++)
	{
		Point left(int((rand() % (max - min)) + min), segment_id, 0);			 //for easy reading
		Point right(int(left.x + rand() % (max - min) + min + 1), segment_id, 1);//for easy reading
		global_line.push_back(left);
		global_line.push_back(right);
		segment_id++;
	}
	return global_line;
}

void Print(vector<Point> &line)
{
	for (int i = 0; i < line.size(); i++)
	{
		line[i].Print(true);
	}
}
void Print(int* line, int size)
{
	for (int i = 0; i < size; i++)
	{
		cout << line[i] << " ";
	}
	cout << endl;
}
void Print(vector<int>& solution)
{
	for (int i = 0; i < solution.size(); i++)
	{
		cout << solution[i] << " ";
	}
	cout << endl;
}

void main(int argc, char** argv)
{
	double TimeStart, TimeEnd;
	int rank;
	int RankSize;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &RankSize);
	Point* line;//part of global_line. Unique for each process.
	vector<Point> global_line;//for sync and duplicate deleting; contains all segments
	int SegmentsCount;//# segments for each process
	int size;//total # of segments
	int* gather_solution;//serves for MPI_Gather
	vector<int> Final_Solution;//final minimal solution
	if (rank == 0)
	{
		int max, min;
		max = 10;
		min = 0;
		size = atoi(argv[1]);

		//[ Generating data ]
		global_line = GenerateGlobalLine(size, max, min);
		
		//Print(global_line);
		//[ Sequential ]
		TimeStart = MPI_Wtime();
		line = new Point[size*2];
		for (int i = 0; i < size * 2; i++)
		{
			line[i] = global_line[i];
		}
		cout << endl;
		qsort(line, size*2, sizeof(Point), comparePoint);
		vector<int> solution = Cover(line, size);
		TimeEnd = MPI_Wtime();
		Print(solution);
		cout << "time = " << TimeEnd - TimeStart << endl;
		delete[] line;

		//[ MPI ]
		TimeStart = MPI_Wtime();
		SegmentsCount = DataDistr(line, global_line, size, RankSize, max, min);
	}
	else
	{
		//Synced with DataDistr()::start
		MPI_Recv(&SegmentsCount, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		line = new Point[SegmentsCount*2];
		for (int i = 0; i < SegmentsCount*2; i++)
		{
			MPI_Recv(&line[i], 3, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		//Synced with DataDistr()::end
	}

	//Shared section
	qsort(line, SegmentsCount * 2, sizeof(Point), comparePoint);//local sorting
	vector<int> solution = Cover(line, SegmentsCount);//local covering
	int solution_size = solution.size();
	int recv_buff;//Total solution size
	MPI_Reduce(&solution_size, &recv_buff, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	//Gather solutions in one
	if (rank == 0)
	{
		gather_solution = new int[recv_buff];
		for (int i = 0; i < solution.size(); i++)//fill gather_solution with root solution
			gather_solution[i] = solution[i];
		for (int i = solution.size(); i < recv_buff; i++)//fill gather_solution with other solutions
			MPI_Recv(&gather_solution[i], 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	else
	{
		for (int i = 0; i < solution.size(); i++)//fill gather_solution with other solutions
			MPI_Send(&solution[i], 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}

	//Delete extra points
	//Алгоритм: По global_line и gather_solution строю картину происходящего: какие отрезки покрыти и сколькими точками.
	//Если при удалении точки один из отрезков остается без точек покрытий, то удалять нельзя!
	//Иначе, добавляем в ответ - Final_Solution
	if (rank == 0)
	{
		qsort(gather_solution, recv_buff, sizeof(int), compareInt);
		vector<vector <int> > segments_list;//i - №segment in global_line
		segments_list.resize(size);
		
		//Fill segments_list
		for (int i = 0; i < recv_buff; i++)
		{
			for (int j = 0; j < global_line.size()-1; j = j+2)
			{
				if ((gather_solution[i] >= global_line[j].x) && (gather_solution[i] <= global_line[j + 1].x))
				{
					segments_list[global_line[j].segment].push_back(gather_solution[i]);
				}
			}
		}
		//Delete extra
		for (int i = 0; i < recv_buff; i++)//For all solutions
		{
			bool CanDelete = true;
			for (int j = 0; j < segments_list.size(); j++)//For all segments
			{
				for (int k = 0; k < segments_list[j].size(); k++)//For all segment's points
				{
					if (segments_list[j][k] == gather_solution[i])//If our solution contains in the segment
						if (segments_list[j].size() == 1)//But if segment has only one point (this solution)
						{
							CanDelete = false;
							goto NextSolution;//Skip and add to Final_Solution
						}
				}
			}

		NextSolution:
			if (!CanDelete)
			{
					Final_Solution.push_back(gather_solution[i]);
			}
			else
			{
				for (int j = 0; j < segments_list.size(); j++)//For all segments
				{
					for (int k = 0; k < segments_list[j].size(); k++)//For all segment's points
					{
						if (segments_list[j][k] == gather_solution[i])//If our solution contains in the segment
								segments_list[j].erase(segments_list[j].begin() + k);
					}
				}
			}
		
		}

		Final_Solution.erase(unique(Final_Solution.begin(), Final_Solution.end()), Final_Solution.end());
		TimeEnd = MPI_Wtime();
		Print(Final_Solution);
		cout << "time = " << TimeEnd - TimeStart << endl;

		delete[] gather_solution;
	}
	delete[] line;
	MPI_Finalize();
}