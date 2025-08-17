/*
 * IMS_v1.0.h
 *
 *  Created on: 2025年3月6日
 *      Author: LinLin
 */

#ifndef IMS_V1_0_H_
#define IMS_V1_0_H_

#include <algorithm>
#include <assert.h>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <libgen.h>
#include <limits.h>
#include <list>
#include <map>
#include <math.h>
#include <set>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <unistd.h>
#include <utility>
#include <vector>
//#include <fcntl.h>
//#include <sys/file.h>

#include "tools.h"

#define MAX_VALUE 99999999
string root_path;

class Node{
public:
	int idx;
	vector<pair<int, int>> edges;

	Node(int val):idx(val){}

	void add_relation(int target, int weight) {
        edges.emplace_back(target, weight);
    }
};

class Graph{
public:
	int k;
	int nnode;
	int nedge;
	vector<Node> nodes;

	Graph() : k(-1), nnode(-1), nedge(-1){}

	// 构造函数（通过读图构造）
	Graph(string filename, int nk)
	{
		ifstream fin(filename);		// 打开文件
		// 判断是否成功打开文件
		if (!fin.is_open())
		{
			cerr << "Can not open the file!" << filename << endl;
			exit(-2);
		}

		// 检查文件流是否处于错误状态
		if (fin.fail())
		{
			cerr << "Error occurred during file operation." << filename << endl;
			exit(-3);
		}

		// 判断文件是否为空
		if (fin.eof())
		{
			cerr << "Empty file" << filename << endl;
			exit(-4);
		}

		// 开始读图
		string first_line;
		getline(fin, first_line);
		istringstream first_ss(first_line);

		int elem1, elem2;
		first_ss >> elem1 >> elem2;
		k = nk;
		nnode = elem1;
		nedge = elem2;

		// 初始化节点列表
		nodes.reserve(nnode);
		for (int i = 0; i < nnode; ++i) {
			nodes.emplace_back(i); // 假设节点编号从0开始
		}

		// 逐行读取邻接关系
		string line;
		int source_index = -1;
		while (getline(fin, line)) {
			if (line.empty()) continue; // 跳过空行

			istringstream ss(line);
			vector<string> parts;
			string part;

			// 分割行内容
			while (ss >> part) {
				parts.push_back(part);
			}

			// 跳过行首冗余信息
			source_index += 1;
			if (source_index >= nnode) {
				cerr << "错误: 超出顶点数量限制" << endl;
				break;
			}

			// 解析邻接关系（从parts[1]开始）
			for (size_t i = 1; i < parts.size(); i += 2) { // i=1跳过行首编号
				if (i + 1 >= parts.size()) {
					cerr << "格式错误: 第 " << source_index + 1 << " 行数据不完整" << endl;
					break;
				}

				// 转换目标节点编号（假设输入文件从1开始编号）
				int target_index = stoi(parts[i]) - 1;
				int weight = stoi(parts[i + 1]);

				// 添加邻接关系
				if (target_index >= 0 && target_index < nnode) {
					nodes[source_index].add_relation(target_index, weight);
				} else {
					cerr << "无效的目标节点索引: " << target_index << endl;
				}
			}
		}

		fin.close();
	}


	void print_graph()
	{
		for (int i = 0; i < nnode; i++) { // 遍历所有节点
			cout << nodes[i].idx; // 输出当前节点编号（不带换行）

			// 遍历当前节点的所有邻接关系
			int nnbh = int(nodes[i].edges.size());
			for (int j = 0; j < nnbh; j++) {
				int adjv = nodes[i].edges[j].first;   // 邻接点索引
				int adjw = nodes[i].edges[j].second;  // 边权值

				cout << " " << adjv << " " << adjw;
			}

			// 当前节点处理完毕后换行
			cout << "\n";
		}
	}
};

class Solution{
public:
	int cost;
	double btime;
	vector<int> ptn;		// partioning
	vector<int> sc;			// size_cluster
//	vector<vector<int>>

	// 计算
	int cal_pcost(const Graph &graph, const int &v1)
	{
		int pcost = 0;
		int size = int(graph.nodes[v1].edges.size());

		// 遍历所有邻接点
		for(int i = 0; i < size; i++)
		{
			int v2 = graph.nodes[v1].edges[i].first;
			int w = graph.nodes[v1].edges[i].second;

			if(ptn[v1] == ptn[v2] && w < 0) pcost += abs(w);
			if(ptn[v1] != ptn[v2] && w > 0) pcost += w;
		}

		return pcost;
	}

	// 计算目标函数值
	int cal_cost(const Graph &graph)
	{
		int ccost = 0;
		for(int i = 0; i < graph.nnode; i++)
		{
			int v1 = graph.nodes[i].idx;
			int size = int(graph.nodes[v1].edges.size());
			for(int j = 0; j < size; j++)
			{
				int v2 = graph.nodes[v1].edges[j].first;
				int w = graph.nodes[v1].edges[j].second;

				if(ptn[v1] == ptn[v2] && w < 0)
					ccost += abs(w);

				if(ptn[v1] != ptn[v2] && w > 0)
					ccost +=  w;
			}
		}

		ccost /= 2;
		return ccost;
	}

	// 构造一个空解
	Solution() : cost(MAX_VALUE), btime(0.0), ptn(), sc() {}

	// 构造一个随机初始解，并更新每个分区中的点数和解对应的cost
	Solution(const Graph &graph, clock_t cstime)
		: btime(0.0), ptn(graph.nnode, 0), sc(graph.k, 0) {

#if(DEBUG)
		printf("debug--1.1.1\n");fflush(stdout);
#endif
		// 给每个分区随机分配一个点
		int *randlist = new int[graph.nnode];			// 生成随机序列
		Generate_Rand_List(randlist, graph.nnode);
#if(DEBUG)
		printf("debug--1.1.2\n");fflush(stdout);
#endif
		for(int pid = 0; pid < graph.k; pid++)
		{
//			printf("debug--1.1.2.1, pid=%d\n", pid);fflush(stdout);
			int v = randlist[pid];
//			printf("debug--1.1.2.2, v=%d\n", v);fflush(stdout);
			ptn[v] = pid;
//			printf("debug--1.1.2.3\n");fflush(stdout);
			sc[pid]++;
		}
#if(DEBUG)
		printf("debug--1.1.3\n");fflush(stdout);
#endif
		// 剩余点随机分配
		for(int i = graph.k; i < graph.nnode; i++)
		{
			int v = randlist[i];			// 选点
			int p = rand() % graph.k;				// 选分区
			ptn[v] = p;						// 赋值
			sc[p]++;						// 更新分区大小
		}

		delete[] randlist;

		// 计算cost
#if(DEBUG)
		printf("debug--1.1.3\n");fflush(stdout);
#endif
		cost = cal_cost(graph);
		btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
	}

	// 通过复制另一个解来构造新解
	Solution(const Solution &sol1)
		: cost(sol1.cost), btime(sol1.btime), ptn(sol1.ptn), sc(sol1.sc) {}

	// 复制另一个解
	void cpy(const Solution &sol1)
	{
		cost = sol1.cost;
		ptn = sol1.ptn;					// 深拷贝，改变该数组不会影响到sol1.ptn
		sc = sol1.sc;
		btime = sol1.btime;
	}

	bool verify(const Graph &graph)
	{
		bool is_true = false;

		int vcost = 0;
		for(int i = 0; i < graph.nnode; i++)
		{
			int v1 = graph.nodes[i].idx;
			int size = int(graph.nodes[v1].edges.size());
			for(int j = 0; j < size; j++)
			{
				int v2 = graph.nodes[v1].edges[j].first;
				int w = graph.nodes[v1].edges[j].second;

				if(ptn[v1] == ptn[v2] && w < 0)
					vcost += abs(w);

				if(ptn[v1] != ptn[v2] && w > 0)
					vcost +=  w;
			}
		}

		vcost /= 2;
		return vcost;

		if(vcost == cost)
			is_true = true;

		return is_true;
	}

	// 使用默认析构函数
	~Solution() {}



};

class Gain_node {
public:
	int elem1;						// 要移动的点
	int elem2;						// 移动到的分区
	int elem3;						// new pid
	int delta;
	int type;						// type=0表示是移动，type=1表示是交换，type=2表示是push

	Gain_node() : elem1(-1), elem2(-1), elem3(-1), delta(MAX_VALUE), type(-1){}
	Gain_node(int e1, int e2, int e3, int delta, int type) : elem1(e1), elem2(e2), elem3(e3), delta(delta), type(type){}

	void clear()
	{
		elem1 = -1;
		elem2 = -1;
		elem3 = -1;
		delta = MAX_VALUE;
		type = -1;
	}

//	void set(int vertex, int cluster, double move_gain, int type)
//	{
//		this->vertex = vertex;
//		this->cluster = cluster;
//		this->move_gain = move_gain;
//	}
	~Gain_node() {}
};

//Graph graph;
clock_t start, end;
double *each_run_rlt;
double *each_run_time;
double *each_hit_time;
double avg_cost, avg_time, avg_htime, std_dev;
double sum_avg_cost = 0, sum_avg_time = 0, sum_avg_htime = 0;

Solution glb_best_sol;									// global best sol
Solution ml_best_sol;
Solution ils_best_sol;
Solution ts_best_sol;									// tabu search best sol

double glb_best_cost;
double rd_best_cost;
double ils_best_cost;
double ts_best_cost;


//char filename[1000] = "./Benchmarks/R0_signed_graph/slashdot-zoo.graph";
char filename[1000] = "./instances/soc-sign-epinions.graph";

int seed = 14;
double time_limit = 200;
int K = 2;
int runs = 1;

int *asc_nodes1, *asc_nodes2, *asc_nodes3, *asc_nodes4;
int *desc_nodes1, *desc_nodes2, *desc_nodes3, *desc_nodes4;
int *desc_swap_nodes1;
vector<vector<int>> Descending_swap_v;
int **pos_gamma, **neg_gamma;

// 记录最好的结果和最好的时间
int rbcost;
double rbtime;


// IMS算法参数
int max_nipv = 10;				// IMS run length 50, 20
double pct_sp = 0.05;			// percentage of strong perturb
double pct_wp = 0.01;			// percentage of weak perturb
int sp, wp;

/*
 * 读参数
 */
void Read_Parameters(int argc, char **argv)
{
//	for(int i = 1; i < argc; i++)
//	{
//		cout << argv[i] << endl;
//	}

	for (int i = 1; i < argc; i += 2)							// [0]表示'-'，[1]表示类型，[2]表示数据
	{
		if (argv[i][0] != '-')
		{
			exit(0);
		}
		else if (argv[i][1] == 'i')  // The file name
		{
			strncpy(filename, argv[i + 1], 1000);
//		cerr << filename << endl;
		}
		else if (argv[i][1] == 's') // The maximum time
			seed = atoi(argv[i + 1]);
		else if (argv[i][1] == 'r') // The maximum time
			time_limit = atof(argv[i + 1]);			// ils 的迭代次数
		else if (argv[i][1] == 'k') //not used
			K = atoi(argv[i + 1]);			// ts的最大迭代次数比率
//		else if (argv[i][1] == 't')
//			tt_pct = atof(argv[i + 1]);				// 禁忌期比率
//		else if (argv[i][1] == 'd')
//			td_pct = atof(argv[i + 1]);				// 连续未改进次数比率
//		else if (argv[i][1] == 'n') //not used
//			np_pct = atof(argv[i + 1]);			// ts的最大迭代次数比率
//		else if (argv[i][1] == 'b')
//			param_beta = atof(argv[i + 1]);				// 连续未改进次数比率
//		else if (argv[i][1] == 'r')
//			param_brate = atof(argv[i + 1]);				// 禁忌期比率


	}
	// check parameters
//	if (strlen(filename) == 0)
//	{
//		cerr << "No input data" << endl;
//		exit(1);
//	}
}


/*
 * 输出标题
 */
void output_header(int run)
{
	char outputfile[1000];
	char *graph_name = basename(filename);
	sprintf(outputfile, "%soutput_dir/results_%s_%d.txt", root_path.c_str(), graph_name, K);
	FILE *opf = fopen(outputfile, "a");
	if (!opf) {
		perror("Failed to open output file");
		exit(-1);
	}
	fprintf(opf, "\n\n------------------------- IMS Round %d -------------------------\n", run);
	fclose(opf);
}


/*
 * 分配内存
 */
void allocate_memory(const Graph &graph)
{
	Descending_swap_v.resize(graph.k);
	asc_nodes1 = new int[graph.nnode];
	desc_nodes1 = new int[graph.nnode];
	asc_nodes2 = new int[graph.nnode];
	desc_nodes2 = new int[graph.nnode];
	asc_nodes3 = new int[graph.nnode];
	desc_nodes3 = new int[graph.nnode];
	asc_nodes4 = new int[graph.nnode];
	desc_nodes4 = new int[graph.nnode];
	desc_swap_nodes1 = new int[graph.nnode];

	pos_gamma = new int*[graph.k];
	neg_gamma = new int*[graph.k];
	for(int i = 0; i < graph.k; i++)
	{
		pos_gamma[i] = new int[graph.nnode];
	}
	for(int i = 0; i < graph.k; i++)
	{
		neg_gamma[i] = new int[graph.nnode];
	}

}


/*
 * 释放内存
 */
void free_memory(const Graph &graph)
{
	delete[] asc_nodes1;
	delete[] asc_nodes2;
	delete[] asc_nodes3;
	delete[] asc_nodes4;
	delete[] desc_nodes1;
	delete[] desc_nodes2;
	delete[] desc_nodes3;
	delete[] desc_nodes4;
	delete[] desc_swap_nodes1;

	for(int i = 0;  i < graph.k; i++)
		delete[] pos_gamma[i];
	delete[] pos_gamma;

	for(int i = 0; i < graph.k; i++)
		delete[] neg_gamma[i];
	delete[] neg_gamma;
}


/*
 * RH算法
 */
void multistart_relocation_heuristic(const Graph &graph, Solution &bsol, double tlimit)
{
	// 开始重启迭代
	int riter = 0;								// restart_iter
	clock_t cstime = clock();
	while ((clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC) < tlimit)
	{
#if(DEBUG)
		printf("debug--1.1\n");fflush(stdout);
#endif
		// 生成初始解
		Solution csol = Solution(graph, cstime);		// cur_sol
		Solution tsol = Solution(csol);			// temp_sol

#if(DEBUG)
		printf("debug--1.2\n");fflush(stdout);
#endif
		bool improved = true;
		while(improved)
		{
			improved = false;
#if(DEBUG)
			printf("debug--1.3\n");fflush(stdout);
#endif
			for(int i = 0; i < graph.nnode; i++)
			{
				int vtx1 = i;
#if(DEBUG)
				printf("debug--1.5\n");fflush(stdout);
#endif
				for(int clst2 = 0; clst2 < graph.k; clst2++)
				{
					int clst1 = csol.ptn[vtx1];
					if(csol.sc[clst1] == 1 || clst2 == clst1)
						continue;

					tsol.ptn[vtx1] = clst2;
					if(tsol.cal_pcost(graph, vtx1) < csol.cal_pcost(graph, vtx1))
					{
						csol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
						csol.sc[clst1] -= 1;
						csol.ptn[vtx1] = clst2;
						csol.sc[clst2] += 1;
						improved = true;
					}
					else
					{
						tsol.cpy(csol);
					}
				}
			}
		}

		// 计算cost
#if(DEBUG)
		printf("debug--1.6\n");fflush(stdout);
#endif
		csol.cost = csol.cal_cost(graph);
		if(csol.cost < rbcost)
		{
			rbcost = csol.cost;
			rbtime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
		}
		if(csol.cost < bsol.cost)
		{
//			printf("RH round %d, time=%.4f, ccost=%d, improve=%d\n", riter, (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC), csol.cost, bsol.cost - csol.cost);fflush(stdout);
			bsol.cpy(csol);
		}
		riter++;
#if(DPTN)
		printf("RH round %d, time=%.4f, ccost=%d, rbcost=%d\n", riter, (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC), csol.cost, rbcost);fflush(stdout);
#endif
	}

//	printf("=========================== RH END ===========================\n");fflush(stdout);

	// 写入文件
#if(DEBUG)
	printf("debug--1.7\n");fflush(stdout);
#endif
	char outputfile[1000];
	char *graph_name = basename(filename);
	sprintf(outputfile, "%soutput_dir/results_%s_%d.txt", root_path.c_str(), graph_name, K);
	FILE *opf = fopen(outputfile, "a");
	if (!opf) {
		perror("Failed to open RH output file");
		exit(-1);
	}
	fprintf(opf, "RH restarts:%d, RH best cost=%d, RH best time=%.4f\n", riter, bsol.cost, bsol.btime);
	fclose(opf);
}


/*
 * 从文件中读进来初始解
 */
void read_RH_sol(const Graph &graph, Solution &csol, char *instancefile, const int &fno)
{
	// 读图
	char *graph_name = basename(instancefile);  // 获取当前图的名称
	char rltfile[1000];			// 用于存储文件名
	if(tuning)
		sprintf(rltfile, "../init_sols/%s_%d_%d", graph_name, graph.k, fno);	// 生成新的前缀zzzzmutex&graph_name&tid（每个时间块文件各自拥有一个互斥文件）
	else
		sprintf(rltfile, "%sinit_sols/%s_%d_%d", root_path.c_str(), graph_name, graph.k, fno);
	// 打开解文件
	ifstream sol_file(rltfile);
	if (!sol_file.is_open()) {
		cerr << "无法打开解文件: " << rltfile << endl;
	}

	csol.btime = 0.0;
	// 分配内存
	csol.ptn.resize(graph.nnode, 0);
	csol.sc.resize(graph.k, 0);

	// 读取cost
	sol_file >> csol.cost;

	// 读取solution
	for(int i = 0; i < graph.nnode; i++)
	{
		sol_file >> csol.ptn[i];
		csol.sc[csol.ptn[i]]++;
	}

	// 更新rbtime
	if(csol.cost < rbcost)
	{
		rbcost = csol.cost;
		rbtime = (clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
	}

	sol_file.close();
}

/*
 * 计算常用指标并输出：best cost, average cost, average time, hit, dev
 */
void cal_indicators(int *each_run_rlt, int *each_run_time, const int &runs, int &bcost, double &avg_cost, double &avg_time)
{
	double sum_cost = 0, sum_time = 0;
	std_dev = 0;

	bcost = MAX_VALUE;
	for (int i = 0; i < runs; i++)
	{
		sum_cost += each_run_rlt[i];
		sum_time += each_run_time[i];

		if(each_run_rlt[i] < bcost)
			bcost = each_run_rlt[i];
	}

	avg_cost = sum_cost / runs;
	avg_time = sum_time / runs;

	sum_avg_cost += avg_cost;
	sum_avg_time += avg_time;

	for (int i = 0; i < runs; i++)
		std_dev += pow(each_run_rlt[i] - avg_cost, 2) / runs;

	std_dev = sqrt(std_dev);

	// 输出
//	cout << "best cost: " << glb_best_cost << "\nhit:" << hit << "\navg_cost: " << avg_cost << "\navg_time: " << avg_time << "\nstd_dev: " << std_dev << endl;
}


/*
 * 按 结点的影响力=∑|正边权|+|负边权| 确定节点的遍历顺序
 */
void calculate_v1(const Graph &graph)
{
	// 重置
	for(int i = 0; i < graph.nnode; i++)
	{
		asc_nodes1[i] = i;
		desc_nodes1[i] = i;
	}

	// 计算结点影响力
	int *obj_value1 = new int[graph.nnode];
	int *obj_value2 = new int[graph.nnode];
	memset(obj_value1, 0, sizeof(int)*graph.nnode);

	for(int i = 0; i < graph.nnode; i++)
	{
		int v1 = graph.nodes[i].idx;
		int size = int(graph.nodes[i].edges.size());
		for(int j = 0; j < size; j++)
		{
			obj_value1[v1] += abs(graph.nodes[i].edges[j].second);
		}
	}
	memcpy(obj_value2, obj_value1, sizeof(int)*graph.nnode);

	// 排序
	Quick_Sort_asc(asc_nodes1, obj_value1, 0, graph.nnode-1);
	Quick_Sort_desc(desc_nodes1, obj_value2, 0, graph.nnode-1);

	// 释放内存
	delete[] obj_value1;
	delete[] obj_value2;
}


/*
 * 按 结点的影响力=∑|正边权| 确定节点的遍历顺序
 */
void calculate_v2(const Graph &graph)
{
	// 重置
	for(int i = 0; i < graph.nnode; i++)
	{
		asc_nodes2[i] = i;
		desc_nodes2[i] = i;
	}

	// 计算结点影响力
	int *obj_value1 = new int[graph.nnode];
	int *obj_value2 = new int[graph.nnode];
	memset(obj_value1, 0, sizeof(int)*graph.nnode);

	for(int i = 0; i < graph.nnode; i++)
	{
		int v1 = graph.nodes[i].idx;
		int size = int(graph.nodes[i].edges.size());
		for(int j = 0; j < size; j++)
		{
			if(graph.nodes[i].edges[j].second > 0)
				obj_value1[v1] += graph.nodes[i].edges[j].second;
		}
	}
	memcpy(obj_value2, obj_value1, sizeof(int)*graph.nnode);

	// 排序
	Quick_Sort_asc(asc_nodes2, obj_value1, 0, graph.nnode-1);
	Quick_Sort_desc(desc_nodes2, obj_value2, 0, graph.nnode-1);

	// 释放内存
	delete[] obj_value1;
	delete[] obj_value2;
}


/*
 * 按 结点的影响力=∑|负边权| 确定节点的遍历顺序
 */
void calculate_v3(const Graph &graph)
{
	// 重置
	for(int i = 0; i < graph.nnode; i++)
	{
		asc_nodes3[i] = i;
		desc_nodes3[i] = i;
	}

	// 计算结点影响力
	int *obj_value1 = new int[graph.nnode];
	int *obj_value2 = new int[graph.nnode];
	memset(obj_value1, 0, sizeof(int)*graph.nnode);

	for(int i = 0; i < graph.nnode; i++)
	{
		int v1 = graph.nodes[i].idx;
		int size = int(graph.nodes[i].edges.size());
		for(int j = 0; j < size; j++)
		{
			if(graph.nodes[i].edges[j].second < 0)
				obj_value1[v1] += abs(graph.nodes[i].edges[j].second);
		}
	}
	memcpy(obj_value2, obj_value1, sizeof(int)*graph.nnode);

	// 排序
	Quick_Sort_asc(asc_nodes3, obj_value1, 0, graph.nnode-1);
	Quick_Sort_desc(desc_nodes3, obj_value2, 0, graph.nnode-1);

	// 释放内存
	delete[] obj_value1;
	delete[] obj_value2;
}


/*
 * 根据当前解，按节点的分区间的|正边权|+分区里的|负边权|之和确定节点的遍历顺序
 */
void calculate_v4(const Graph &graph, const Solution &csol)
{
	// 重置
	for(int i = 0; i < graph.nnode; i++)
	{
		asc_nodes4[i] = i;
		desc_nodes4[i] = i;
	}

	// 计算结点影响力
	int *obj_value1 = new int[graph.nnode];
	int *obj_value2 = new int[graph.nnode];
	memset(obj_value1, 0, sizeof(int)*graph.nnode);

	for(int i = 0; i < graph.nnode; i++)
	{
		int v1 = graph.nodes[i].idx;
		int size = int(graph.nodes[i].edges.size());
		for(int j = 0; j < size; j++)
		{
			if(graph.nodes[i].edges[j].second > 0 && csol.ptn[i] != csol.ptn[j])
				obj_value1[v1] += graph.nodes[i].edges[j].second;
			if(graph.nodes[i].edges[j].second < 0  && csol.ptn[i] == csol.ptn[j])
				obj_value1[v1] += abs(graph.nodes[i].edges[j].second);
		}
	}
	memcpy(obj_value2, obj_value1, sizeof(int)*graph.nnode);

	// 排序
	Quick_Sort_asc(asc_nodes4, obj_value1, 0, graph.nnode-1);
	Quick_Sort_desc(desc_nodes4, obj_value2, 0, graph.nnode-1);

	// 释放内存
	delete[] obj_value1;
	delete[] obj_value2;
}

/*
 * 根据csol初始化Gamma表
 */
void init_Gamma(const Graph &graph, const Solution &csol)
{
	for(int i = 0; i < graph.k; i++)
	{
		memset(neg_gamma[i], 0, sizeof(int)*graph.nnode);
		memset(pos_gamma[i], 0, sizeof(int)*graph.nnode);
	}

	for(int v1 = 0; v1 < graph.nnode; v1++)
	{
		int size = int(graph.nodes[v1].edges.size());
		for(int j = 0; j < size; j++)
		{
			int v2 = graph.nodes[v1].edges[j].first;
			int w = graph.nodes[v1].edges[j].second;
			int p = csol.ptn[v2];

			if(w < 0)
				neg_gamma[p][v1] += abs(w);
			else
				pos_gamma[p][v1] += w;
		}
	}
}

/*
 * 更新Gamma表，并移动
 */
void update_Gamma(const Graph &graph, Solution &csol, const Gain_node &gnode)
{
	if(gnode.type == 0)
	{
		int v1 = gnode.elem1;
		int p1 = csol.ptn[v1];
		int p2 = gnode.elem2;
		int size = graph.nodes[v1].edges.size();
		for(int j = 0; j < size; j++)
		{
			int v2 = graph.nodes[v1].edges[j].first;
			int w = graph.nodes[v1].edges[j].second;

			if(w < 0)
			{
				neg_gamma[p1][v2] -= abs(w);
				neg_gamma[p2][v2] += abs(w);
			}
			else
			{
				pos_gamma[p1][v2] -= w;
				pos_gamma[p2][v2] += w;
			}
		}

		// 移动
		csol.sc[p1]--;
		csol.ptn[v1] = p2;
		csol.sc[p2]++;
		csol.cost += gnode.delta;
	}
	else if(gnode.type == 1)
	{
		int v1 = gnode.elem1;
		int v2 = gnode.elem2;
		int p1 = csol.ptn[v1];
		int p2 = csol.ptn[v2];

		// v1 → p2
		int size = graph.nodes[v1].edges.size();
		for(int j = 0; j < size; j++)
		{
			int v2 = graph.nodes[v1].edges[j].first;
			int w = graph.nodes[v1].edges[j].second;

			if(w < 0)
			{
				neg_gamma[p1][v2] -= abs(w);
				neg_gamma[p2][v2] += abs(w);
			}
			else
			{
				pos_gamma[p1][v2] -= w;
				pos_gamma[p2][v2] += w;
			}
		}

		// v2 → p1
		size = graph.nodes[v2].edges.size();
		for(int j = 0; j < size; j++)
		{
			int v1 = graph.nodes[v2].edges[j].first;
			int w = graph.nodes[v2].edges[j].second;

			if(w < 0)
			{
				neg_gamma[p2][v1] -= abs(w);
				neg_gamma[p1][v1] += abs(w);
			}
			else
			{
				pos_gamma[p2][v1] -= w;
				pos_gamma[p1][v1] += w;
			}
		}

		// 移动
		csol.ptn[v1] = p2;
		csol.ptn[v2] = p1;
		csol.cost += gnode.delta;
	}

	// 验证
	if(csol.cost != csol.cal_cost(graph))
	{
		printf("ccost=%d, vcost=%d", csol.cost, csol.cal_cost(graph)); fflush(stdout);
		exit(-2);
//		assert(csol.cost == csol.cal_cost(graph));
	}



}


/*
 * 使用gamma表进行快速更新的
 */
Solution local_search(const Graph &graph, const Solution &bsol, clock_t cstime)
{
//	printf("========================== LS BEGIN ==========================\n");fflush(stdout);

	Solution csol = Solution(bsol);
	init_Gamma(graph, csol);
#if(DEBUG)
	printf("debug--1.2\n");fflush(stdout);
#endif
	bool improved = true;
	int iter = 0;
	while(improved)
	{
		improved = false;
#if(DEBUG)
		printf("debug--1.3\n");fflush(stdout);
#endif
		for(int i = 0; i < graph.nnode; i++)
		{
			int vtx1 = i;
			int clst1 = csol.ptn[vtx1];
			int size = int(graph.nodes[i].edges.size());
			if(csol.sc[clst1] == 1 || size == 0)
				continue;
#if(DEBUG)
			printf("debug--1.5\n");fflush(stdout);
#endif
			for(int clst2 = 0; clst2 < graph.k; clst2++)
			{
				clst1 = csol.ptn[vtx1];
				if(csol.sc[clst1] == 1 || clst2 == clst1)
					continue;

				int delta = neg_gamma[clst2][vtx1] - neg_gamma[clst1][vtx1] - pos_gamma[clst2][vtx1] + pos_gamma[clst1][vtx1];

				if(delta < 0)
				{
					Gain_node gn = Gain_node(vtx1, clst2, -1, delta, 0);
					update_Gamma(graph, csol, gn);
					if(csol.cost < rbcost)
					{
						rbcost = csol.cost;
						rbtime = (clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
					}
					csol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
//					printf("LS1 round %d, time=%.4f, ccost=%d, rbcost=%d\n", iter, (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC), csol.cost, rbcost);fflush(stdout);
					improved = true;
					assert(csol.cost == csol.cal_cost(graph));
				}
			}
		}
		iter++;
#if(DPTN)
		printf("LS1 round %d, time=%.4f, ccost=%d, rbcost=%d\n", iter, (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC), csol.cost, rbcost);fflush(stdout);
#endif
	}

//	printf("=========================== LS END ===========================\n");fflush(stdout);

	return csol;
}


/*
 * 使用gamma表进行快速更新，使用遍历顺序0
 */
Solution biased_local_search1(const Graph &graph, const Solution &bsol, clock_t cstime)
{
//	printf("========================== BLS1 BEGIN ==========================\n");fflush(stdout);

	Solution csol = Solution(bsol);
	init_Gamma(graph, csol);
#if(DEBUG)
	printf("debug--1.2\n");fflush(stdout);
#endif
	bool improved = true;
	int iter = 0;
	while(improved)
	{
		improved = false;
#if(DEBUG)
		printf("debug--1.3\n");fflush(stdout);
#endif
		for(int i = 0; i < graph.nnode; i++)
		{
			int vtx1 = asc_nodes1[i];
			int clst1 = csol.ptn[vtx1];
			int size = int(graph.nodes[i].edges.size());
			if(csol.sc[clst1] == 1 || size == 0)
				continue;
#if(DEBUG)
			printf("debug--1.5\n");fflush(stdout);
#endif
			for(int clst2 = 0; clst2 < graph.k; clst2++)
			{
				clst1 = csol.ptn[vtx1];
				if(csol.sc[clst1] == 1 || clst2 == clst1)
					continue;

				int delta = neg_gamma[clst2][vtx1] - neg_gamma[clst1][vtx1] - pos_gamma[clst2][vtx1] + pos_gamma[clst1][vtx1];

				if(delta < 0)
				{
					Gain_node gn = Gain_node(vtx1, clst2, -1, delta, 0);
					update_Gamma(graph, csol, gn);
					if(csol.cost < rbcost)
					{
						rbcost = csol.cost;
						rbtime = (clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
					}
					csol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
					assert(csol.cost == csol.cal_cost(graph));
					improved = true;
				}
			}
		}
		iter++;
#if(DPTN)
		printf("BLS1 round %d, time=%.4f, ccost=%d, rbcost=%d\n", iter, (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC), csol.cost, rbcost);fflush(stdout);
#endif

	}

//	printf("=========================== BLS1 END ===========================\n");fflush(stdout);

	return csol;
}


/*
 * 使用gamma表进行快速更新，使用遍历顺序2
 */
Solution biased_local_search2(const Graph &graph, const Solution &bsol, clock_t cstime)
{
//	printf("========================== BLS2 BEGIN ==========================\n");fflush(stdout);

	Solution csol = Solution(bsol);
	init_Gamma(graph, csol);
#if(DEBUG)
	printf("debug--1.2\n");fflush(stdout);
#endif
	bool improved = true;
	int iter = 0;
	while(improved)
	{
		improved = false;
#if(DEBUG)
		printf("debug--1.3\n");fflush(stdout);
#endif
		for(int i = 0; i < graph.nnode; i++)
		{
			int vtx1 = asc_nodes2[i];
			int clst1 = csol.ptn[vtx1];
			int size = int(graph.nodes[i].edges.size());
			if(csol.sc[clst1] == 1 || size == 0)
				continue;
#if(DEBUG)
			printf("debug--1.5\n");fflush(stdout);
#endif
			for(int clst2 = 0; clst2 < graph.k; clst2++)
			{
				clst1 = csol.ptn[vtx1];
				if(csol.sc[clst1] == 1 || clst2 == clst1)
					continue;

				int delta = neg_gamma[clst2][vtx1] - neg_gamma[clst1][vtx1]
						  - pos_gamma[clst2][vtx1] + pos_gamma[clst1][vtx1];

				if(delta < 0)
				{
					Gain_node gn = Gain_node(vtx1, clst2, -1, delta, 0);
					update_Gamma(graph, csol, gn);
					if(csol.cost < rbcost)
					{
						rbcost = csol.cost;
						rbtime = (clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
					}
					csol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
					assert(csol.cost == csol.cal_cost(graph));
					improved = true;
				}
			}
		}
		iter++;
#if(DPTN)
		printf("BLS2 round %d, time=%.4f, ccost=%d, rbcost=%d\n", iter, (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC), csol.cost, rbcost);fflush(stdout);
#endif
	}

//	printf("=========================== BLS2 END ===========================\n");fflush(stdout);

	return csol;
}


/*
 * 使用gamma表进行快速更新，使用遍历顺序3
 */
Solution biased_local_search3(const Graph &graph, const Solution &bsol, clock_t cstime)
{
//	printf("========================== BLS3 BEGIN ==========================\n");fflush(stdout);

	Solution csol = Solution(bsol);
	init_Gamma(graph, csol);
#if(DEBUG)
	printf("debug--1.2\n");fflush(stdout);
#endif
	bool improved = true;
	int iter = 0;
	while(improved)
	{
		improved = false;
#if(DEBUG)
		printf("debug--1.3\n");fflush(stdout);
#endif
		for(int i = 0; i < graph.nnode; i++)
		{
			int vtx1 = asc_nodes3[i];
			int clst1 = csol.ptn[vtx1];
			int size = int(graph.nodes[i].edges.size());
			if(csol.sc[clst1] == 1 || size == 0)
				continue;
#if(DEBUG)
			printf("debug--1.5\n");fflush(stdout);
#endif
			for(int clst2 = 0; clst2 < graph.k; clst2++)
			{
				clst1 = csol.ptn[vtx1];
				if(csol.sc[clst1] == 1 || clst2 == clst1)
					continue;

				int delta = neg_gamma[clst2][vtx1] - neg_gamma[clst1][vtx1] - pos_gamma[clst2][vtx1] + pos_gamma[clst1][vtx1];

				if(delta < 0)
				{
					Gain_node gn = Gain_node(vtx1, clst2, -1, delta, 0);
					update_Gamma(graph, csol, gn);
					if(csol.cost < rbcost)
					{
						rbcost = csol.cost;
						rbtime = (clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
					}
					csol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
					assert(csol.cost == csol.cal_cost(graph));
					improved = true;
				}
			}
		}

		iter++;
#if(DPTN)
		printf("BLS3 round %d, time=%.4f, ccost=%d, rbcost=%d\n", iter, (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC), csol.cost, rbcost);fflush(stdout);
#endif
	}

	return csol;
}


/*
 * 找到点x和点y的权重
 */
int get_weight(const Graph &graph, int vtx1, int vtx2) {
    int weight = 0;
    if(graph.nodes[vtx1].edges.size() <= graph.nodes[vtx2].edges.size())
    {
    	int nynbh = int(graph.nodes[vtx2].edges.size());
    	for(int i = 0; i < nynbh; i++)
    	{
    		if(vtx1 == graph.nodes[vtx2].edges[i].first)
    			weight = graph.nodes[vtx2].edges[i].second;
    	}
    }
    else
    {
    	int nxnbh = int(graph.nodes[vtx1].edges.size());
		for(int i = 0; i < nxnbh; i++)
		{
			if(vtx2 == graph.nodes[vtx1].edges[i].first)
				weight = graph.nodes[vtx1].edges[i].second;
		}
    }

    return weight;
}


// Function to calculate swap values
void swap_v(const Graph &graph, const Solution& csol) {
    int swap_n = int(graph.nnode / graph.k);
	vector<vector<pair<int, int>>> swap_v(graph.k);
	for(int vtx1 = 0; vtx1 < graph.nnode; vtx1++)
    {
    	int v_cost = 0;
    	int pid1 = csol.ptn[vtx1];
    	for(int pid2 = 0; pid2 < graph.k; pid2++)
    	{
    		if(pid2 != pid1)
    		{
    			v_cost += pos_gamma[pid2][vtx1];
    			v_cost -= neg_gamma[pid2][vtx1];
    		}
    		else
    		{
    			v_cost -= pos_gamma[pid2][vtx1];
				v_cost += neg_gamma[pid2][vtx1];
    		}
    	}
    	swap_v[pid1].push_back(make_pair(v_cost, vtx1));
    }

	for(int pid = 0; pid < graph.k; pid++)
	{
		sort(swap_v[pid].begin(), swap_v[pid].end(), [](const pair<int, int> &a, const pair<int, double> &b) {
		        return a.first > b.first;
		    });
		int len_dsvk = int(Descending_swap_v[pid].size());  // length of descending swap [v] to [k]
		if(len_dsvk >= swap_n)
		{
			for(int i = 0; i < swap_n; i++)
				Descending_swap_v[pid].push_back(swap_v[pid][i].second);
		}
		else
		{
			for(int i = 0; i < len_dsvk; i++)
				Descending_swap_v[pid].push_back(swap_v[pid][i].second);
		}
	}
}

// Function to perform swap local search
Solution swap_local_search(const Graph &graph, const Solution& bsol, clock_t cstime) {
	Solution csol = Solution(bsol);
	init_Gamma(graph, csol);
#if(DEBUG)
	printf("debug--1.2\n");fflush(stdout);
#endif
	bool improved = true;
	int iter = 0;
	while(improved)
	{
		improved = false;
#if(DEBUG)
		printf("debug--1.3\n");fflush(stdout);
#endif
		for(int pid1 = 0; pid1 < graph.k; pid1++)
		{
			for(const auto& vtx1 : Descending_swap_v[pid1])
			{
				for(int pid2 = pid1+1; pid2 < graph.k;pid2++)
				{
					for(const auto& vtx2 : Descending_swap_v[pid2])
					{
						if(csol.ptn[vtx1] == csol.ptn[vtx2])
							continue;
						int weight = get_weight(graph, vtx1, vtx2);
						int delta = 2*weight
									+ neg_gamma[csol.ptn[vtx2]][vtx1]
									- pos_gamma[csol.ptn[vtx2]][vtx1]
									- neg_gamma[csol.ptn[vtx1]][vtx1]
									+ pos_gamma[csol.ptn[vtx1]][vtx1]
									+ neg_gamma[csol.ptn[vtx1]][vtx2]
									- pos_gamma[csol.ptn[vtx1]][vtx2]
									- neg_gamma[csol.ptn[vtx2]][vtx2]
									+ pos_gamma[csol.ptn[vtx2]][vtx2];
						if(delta < 0)
						{
							csol.cost += delta;
							Gain_node gn = Gain_node(vtx1, vtx2, -1, delta, 1);

							update_Gamma(graph, csol, gn);
							improved = true;
						}
					}
				}
			}
		}

		iter++;
#if(DPTN)
		printf("SLS round %d, time=%.4f, ccost=%d, rbcost=%d\n", iter, (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC), csol.cost, rbcost);fflush(stdout);
#endif
	}

	return csol;
}

/*
 * 扰动
 */
Solution shake(const Graph &graph, Solution &csol, const int np)
{
	for(int i = 0; i < np; i++)
	{
		// 选点
		int vtx1 = rand() % graph.nnode;
		int vtx2 = rand() % graph.nnode;
		while(vtx1 == vtx2)
			vtx2 = rand() % graph.nnode;
		int clst1 = csol.ptn[vtx1];
		int clst2 = csol.ptn[vtx2];

		// 扰动
		int type = rand() % 2;
		if(type == 0 && csol.sc[clst1] > 1)				// move
		{
			csol.sc[clst1]--;
			csol.ptn[vtx1] = clst2;
			csol.sc[clst2]++;
		}
		else if(type == 1)			// swap
		{
			csol.ptn[vtx1] = clst2;
			csol.ptn[vtx2] = clst1;
		}
	}

	csol.cost = csol.cal_cost(graph);
	if(csol.cost < rbcost)
	{
		rbcost = csol.cost;
		rbtime = (clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
	}
	return csol;
}


/*
 * 扰动（写错的版本）
 */
Solution shake1_0(const Graph &graph, Solution &csol, const int np)
{
	calculate_v4(graph, csol);
	for(int i = 0; i < graph.nnode; i++)		// TODO：这里也测一下
	{
		// 选点
		// TODO: 测一下如果改回来会变好吗
		int vtx1 = rand() % graph.nnode;
		int vtx2 = rand() % graph.nnode;
		while(vtx1 == vtx2)
			vtx2 = rand() % graph.nnode;
		int clst1 = csol.ptn[vtx1];
		int clst2 = csol.ptn[vtx2];

		// 扰动
		int type = rand() % 2;
		if(type == 0)				// move
		{
			if(csol.sc[clst1] > 1)
			{
				csol.sc[clst1]--;
				csol.ptn[vtx1] = clst2;
				csol.sc[clst2]++;
			}
		}
		else if(type == 1)			// swap
		{
			csol.ptn[vtx1] = clst2;
			csol.ptn[vtx2] = clst1;
		}

	}

	csol.cost = csol.cal_cost(graph);
	if(csol.cost < rbcost)
	{
		rbcost = csol.cost;
		rbtime = (clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
	}

	return csol;
}


/*
 * 扰动（有偏扰动）TODO
 */
Solution shake1(const Graph &graph, Solution &csol, const int np)
{
	int n = 0;
	calculate_v4(graph, csol);
	for(int i = 0; i < graph.nnode; i++)		// TODO：这里也测一下
	{
		// 选点
		// TODO: 测一下如果改回来会变好吗
		int vtx1 = desc_nodes4[i];
		int clst1 = csol.ptn[vtx1];

		if(csol.sc[clst1] == 1)
			continue;

		int clst2 = rand() % (graph.k-1);
		while(clst2 == clst1)
			clst2 = rand() % (graph.k-1);

		// 扰动
		csol.sc[clst1]--;
		csol.ptn[vtx1] = clst2;
		csol.sc[clst2]++;

		n++;
		if(n >= np)
			break;
	}

	csol.cost = csol.cal_cost(graph);
	if(csol.cost < rbcost)
	{
		rbcost = csol.cost;
		rbtime = (clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
	}

	return csol;
}


/*
 * MS算法
 */
Solution Maxima_search(const Graph &graph, Solution &csol, clock_t cstime, double tlimit)
{
	Solution bsol = Solution(csol);			// best_sol

	assert(csol.cost > 0);
	Solution bsol0 = Solution(local_search(graph, csol, cstime));
	Solution bsol1 = Solution(biased_local_search1(graph, csol, cstime));
	Solution bsol2 = Solution(biased_local_search2(graph, csol, cstime));
	Solution bsol3 = Solution(biased_local_search3(graph, csol, cstime));

	if(bsol0.cost < bsol.cost)
	{
		bsol.cpy(bsol0);
		bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
	}
	if(bsol1.cost < bsol.cost)
	{
		bsol.cpy(bsol1);
		bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
	}
	if(bsol2.cost < bsol.cost)
	{
		bsol.cpy(bsol2);
		bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
	}
	if(bsol3.cost < bsol.cost)
	{
		bsol.cpy(bsol3);
		bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
	}

	// 开始重启迭代
	int non_improve = 0;										// restart_iter
	while (non_improve < max_nipv && (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC) < tlimit)
	{
#if(DEBUG)
		printf("debug--1.1\n");fflush(stdout);
#endif
		shake(graph, csol, wp);

		// 生成初始解
		Solution nsol0 = Solution(local_search(graph, csol, cstime));
		Solution nsol1 = Solution(biased_local_search1(graph, csol, cstime));
		Solution nsol2 = Solution(biased_local_search2(graph, csol, cstime));
		Solution nsol3 = Solution(biased_local_search3(graph, csol, cstime));

		if(nsol0.cost < csol.cost) csol.cpy(nsol0);
		if(nsol1.cost < csol.cost) csol.cpy(nsol1);
		if(nsol2.cost < csol.cost) csol.cpy(nsol2);
		if(nsol3.cost < csol.cost) csol.cpy(nsol3);

		if(csol.cost < bsol.cost)
		{
			if(csol.cost < rbcost)
			{
				rbcost = csol.cost;
				rbtime = (clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
			}
			bsol.cpy(csol);
			bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
			non_improve = 0;
		}
		else
			non_improve++;

//		printf("non improve=%d, time=%.4f, ccost=%d, rbcost=%d\n",
//				non_improve,
//				(clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC),
//				csol.cost,
//				rbcost);
		fflush(stdout);
	}

	return bsol;
}


/*
 * MS算法（使用了swap的版本）
 */
Solution rel_Maxima_search(const Graph &graph, Solution &csol, clock_t cstime, double tlimit)
{
	Solution bsol = Solution(csol);			// best_sol

	assert(csol.cost > 0);
	Solution bsol0 = Solution(local_search(graph, csol, cstime));
	Solution bsol1 = Solution(biased_local_search1(graph, csol, cstime));
	Solution bsol2 = Solution(biased_local_search2(graph, csol, cstime));
	Solution bsol3 = Solution(biased_local_search3(graph, csol, cstime));

	if(bsol0.cost < bsol.cost)
	{
		bsol.cpy(bsol0);
		bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
	}
	if(bsol1.cost < bsol.cost)
	{
		bsol.cpy(bsol1);
		bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
	}
	if(bsol2.cost < bsol.cost)
	{
		bsol.cpy(bsol2);
		bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
	}
	if(bsol3.cost < bsol.cost)
	{
		bsol.cpy(bsol3);
		bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
	}

	swap_v(graph, bsol);
	bsol.cpy(swap_local_search(graph, bsol, cstime));


	// 开始重启迭代
	int non_improve = 0;										// restart_iter
	while (non_improve < max_nipv && (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC) < tlimit)
	{
#if(DEBUG)
		printf("debug--1.1\n");fflush(stdout);
#endif
		shake(graph, csol, wp);

		// 生成初始解
		Solution nsol0 = Solution(local_search(graph, csol, cstime));
		Solution nsol1 = Solution(biased_local_search1(graph, csol, cstime));
		Solution nsol2 = Solution(biased_local_search2(graph, csol, cstime));
		Solution nsol3 = Solution(biased_local_search3(graph, csol, cstime));

		if(nsol0.cost < csol.cost) csol.cpy(nsol0);
		if(nsol1.cost < csol.cost) csol.cpy(nsol1);
		if(nsol2.cost < csol.cost) csol.cpy(nsol2);
		if(nsol3.cost < csol.cost) csol.cpy(nsol3);

		swap_v(graph, csol);
		Solution nsol4 = Solution(swap_local_search(graph, csol, cstime));
		if(nsol4.cost < csol.cost) csol.cpy(nsol3);


		if(csol.cost < bsol.cost)
		{
			if(csol.cost < rbcost)
			{
				rbcost = csol.cost;
				rbtime = (clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
			}
			bsol.cpy(csol);
			bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
			non_improve = 0;
		}
		else
			non_improve++;

//		printf("non improve=%d, time=%.4f, ccost=%d, rbcost=%d\n",
//				non_improve,
//				(clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC),
//				csol.cost,
//				rbcost);
		fflush(stdout);
	}

	return bsol;
}



/*
 * 验证分区和cost
 */
void verify(const Graph &graph, Solution &bsol)
{
	// 1.验证分区数
	for(int i = 0; i < graph.k; i++)
	{
		if(bsol.sc[i] <= 0)
			cerr << "分区" << i << "的点数小于0" << endl;
	}

	// 2.验证cost
	int vcost = bsol.cal_cost(graph);
	if(vcost != bsol.cost)
		cerr << "cost验证未通过" << endl;

	// 3.

}


/*
 * 将解写入文件
 */
void write_IMSS_sol(const Graph &graph, const Solution &csol, char *instancefile, const int &fno) {
    // 获取文件名
    char *graph_name = basename(instancefile);

    // 拼接文件路径
    char rltfile_detail[1000];
    sprintf(rltfile_detail, "%sresults/IMSS_%s_%d_%d", root_path.c_str(), graph_name, graph.k, fno);

    // 用 ofstream 打开文件
    ofstream fout(rltfile_detail);
    if (!fout.is_open()) {
        cerr << "无法打开文件: " << rltfile_detail << endl;
        exit(EXIT_FAILURE);
    }

    // 写入 cost
    fout << "cost=" << csol.cost << endl;
    fout << "cost=" << csol.cost << endl;
    fout << "best_cost=" << rbcost << endl;
    fout << "best_time=" << rbtime << endl;

    // 写入 ptn 数组
    fout << "ptn=";
    for (int i = 0; i < graph.nnode; i++) {
        fout << csol.ptn[i] << " ";
    }
    fout << endl;

    // 写入 sc 数组
    fout << "sc[]=";
    for (int i = 0; i < graph.k; i++) {
        fout << csol.sc[i] << " ";
    }
    fout << endl;

    fout.close();
}


/*
 * IMS
 */
void iterated_maxima_search(const Graph &graph, Solution &csol, double tlimit)
{
#if(DEBUG)
	printf("debug--2.1\n");fflush(stdout);
#endif
	Solution bsol = Solution(csol);						// best_sol
	init_Gamma(graph, csol);
	sp = int(pct_sp * graph.nnode);
	wp = int(pct_wp * graph.nnode);

#if(DEBUG)
	printf("debug--2.2\n");fflush(stdout);
#endif
	// 开始重启迭代
	int iter = 0;										// iter
	clock_t cstime = clock();
#if(DEBUG)
	printf("debug--2.3\n");fflush(stdout);
#endif
	while ((clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC) < tlimit)
	{
//		Solution nsol1 = Solution(Maxima_search(graph, csol, cstime, tlimit));
		Solution nsol1 = Solution(rel_Maxima_search(graph, csol, cstime, tlimit));

#if(DEBUG)
		printf("debug--2.4\n");fflush(stdout);
#endif

		if(nsol1.cost < bsol.cost)
		{
//			printf("IMS strong restarts:%d, time=%f, best cost=%d, improve=%d\n",
//					iter,
//					(clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC),
//					nsol1.cost,
//					bsol.cost - nsol1.cost);
			fflush(stdout);
			bsol.cpy(nsol1);
			bsol.btime = nsol1.btime;
		}

		iter++;
		csol.cpy(shake(graph, nsol1, sp));
	}

//	printf("IMS FINAL restarts:%d, btime=%f, best cost=%d\n",
//			iter,
//			bsol.btime,
//			bsol.cost);
	fflush(stdout);

	// 写入文件
#if(DEBUG)
	printf("debug--2.9\n");fflush(stdout);
#endif
	char outputfile[1000];
	char *graph_name = basename(filename);
	sprintf(outputfile, "%soutput_dir/results_%s_%d.txt", root_path.c_str(), graph_name, K);
	FILE *opf = fopen(outputfile, "a");
	if (!opf) {
		perror("Failed to open output file");
		exit(-1);
	}
	fprintf(opf, "IMS restarts:%d, IMS best cost=%d, IMS best time=%.4f\n", iter, rbcost, rbtime);
	fclose(opf);

	// 释放内存

	csol.cpy(bsol);
}



/*
 * 执行算例instancefile
 */
void IMS_run(char *instancefile, double timelimit)
{
	// 读图
	Graph graph = Graph(instancefile, K);
	allocate_memory(graph);

	calculate_v1(graph);
	calculate_v2(graph);
	calculate_v3(graph);

	// run 10 times algorithms
	int crun = 0;
	while(crun < runs)
	{
		output_header(crun + 1);
		start = clock();
		rbcost = MAX_VALUE;
		rbtime = MAX_VALUE;
		start = clock();

		// 1.RH
#if(DEBUG)
		printf("debug--1\n");fflush(stdout);
#endif
		Solution bsol = Solution();
		read_RH_sol(graph, bsol, instancefile, crun);

		// 2.IMSS
#if(DEBUG)
		printf("debug--2\n");fflush(stdout);
#endif
		iterated_maxima_search(graph, bsol, timelimit);
		verify(graph, bsol);
		write_IMSS_sol(graph, bsol, instancefile, seed);

		printf("IMSS Round %d: best cost=%d\n", crun, rbcost);

		crun++;
	}

	// 释放内存
	free_memory(graph);
}



#endif /* IMS_V1_0_H_ */
