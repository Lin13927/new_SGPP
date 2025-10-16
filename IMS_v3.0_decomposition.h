/*
 * IMS_v3.0.h
 *
 *  Created on: 2025年10月10日
 *      Author: LinLin Zhi
 */

#ifndef IMS_V3_0_H_
#define IMS_V3_0_H_

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

#include <unordered_set>

#include "tools.h"

#define MAX_VALUE 99999999

#define is_Verify 0

#define rMS_mode1 0       // Maxima search模式一：[概率]控制NBL和BL
#define rMS_mode2 0       // Maxima search模式二：先NBL再BL
#define rMS_mode3 0       // Maxima search模式三：先BL再NBL
#define rMS_mode4 1       // Maxima search模式四：只有BL，作为调试base

#define NDBLS_mode1 0  	  // ND-based_LS模式一：Add和Swap两种邻域都从初始解出发
#define NDBLS_mode2 0  	  // ND-based_LS模式二：先Add，再Swap
#define NDBLS_mode3 0  	  // ND-based_LS模式三：先Swap，再Add
#define NDBLS_mode4 1  	  // ND-based_LS模式四：只有Add，+以概率选择Swap

#define is_probSwap 0     // 用[概率]控制是否要Swap
#define is_trueGamma 1    // 使用gamma总表来获取delta

string root_path;

class Node
{
public:
	int idx;
	vector<pair<int, int> > edges;
	// 在C++中，pair是标准模板库STL提供的一个模板类，用于将两个不同类型或相同类型的对象组合成一个整体。


	Node(int val): idx(val)
	{
	}


	void add_relation(int target, int weight)
	{
		edges.emplace_back(target, weight);
	}
};

class Graph
{
public:
	int k;
	int nnode;
	int nedge;
	vector<Node> nodes;


	Graph() : k(-1), nnode(-1), nedge(-1)
	{
	}


	// 构造函数（通过读图构造）
	// TODO 考虑了节点的权重，若不考虑的话，与2019算法结果有差距吗？
	Graph(string filename, int nk)
	{
		cout << "Reading graph from " << filename << "..." << endl;
		ifstream fin(filename); // 打开文件
		// 判断是否成功打开文件
		if (!fin.is_open())
		{
			cerr << "Can not open the file! " << filename << endl;
			exit(-2);
		}

		// 检查文件流是否处于错误状态
		if (fin.fail())
		{
			cerr << "Error occurred during file operation. " << filename << endl;
			exit(-3);
		}

		// 判断文件是否为空
		if (fin.eof())
		{
			cerr << "Empty file " << filename << endl;
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

		cout << "K=" << k << ", V=" << nnode << ", E="<< nedge << endl;

		// 初始化节点列表
		/* C++的vector对象可以通过reserve方法来设置vector对象的容量，通过resize方法来改变vector对象的大小。
		 * reserve所设置的容量指的是vector容器中可容纳元素的预留空间个数，但并不真正创建元素对象。
		 * resize则是直接改变vector容器中元素的个数，并且创建对象或者销毁空间。
		 */
		nodes.reserve(nnode);
		for (int i = 0; i < nnode; ++i)
		{
			/* C++11新增，功能与push_back相同，向vector容器尾部添加一个元素
			 */
			nodes.emplace_back(i); // 假设节点编号从0开始
		}

		// 逐行读取邻接关系
		string line;
		int source_index = -1;
		while (getline(fin, line))
		{
			if (line.empty()) continue; // 跳过空行

			istringstream ss(line);
			vector<string> parts;
			string part;

			// 分割行内容
			while (ss >> part)
			{
				parts.push_back(part);
			}

			// 跳过行首冗余信息
			source_index += 1;
			if (source_index >= nnode)
			{
				cerr << "错误: 超出顶点数量限制" << endl;
				break;
			}

			// 解析邻接关系（从parts[1]开始）
			for (size_t i = 1; i < parts.size(); i += 2)
			{
				// i = 1跳过行首编号
				if (i + 1 >= parts.size())
				{
					cerr << "格式错误: 第 " << source_index + 1 << " 行数据不完整" << endl;
					break;
				}

				// 转换目标节点编号（假设输入文件从1开始编号）
				int target_index = stoi(parts[i]) - 1;
				int weight = stoi(parts[i + 1]);

				// 添加邻接关系
				if (target_index >= 0 && target_index < nnode)
				{
					nodes[source_index].add_relation(target_index, weight);
				}
				else
				{
					cerr << "无效的目标节点索引: " << target_index << endl;
				}
			}
		}

		fin.close();
		cout << "Graph loaded successfully." << endl << endl;
	}


	void print_graph()
	{
		for (int i = 0; i < nnode; i++)
		{
			// 遍历所有节点
			cout << nodes[i].idx; // 输出当前节点编号（不带换行）

			// 遍历当前节点的所有邻接关系
			int nnbh = int(nodes[i].edges.size());
			for (int j = 0; j < nnbh; j++)
			{
				int adjv = nodes[i].edges[j].first; // 邻接点索引
				int adjw = nodes[i].edges[j].second; // 边权值

				cout << " " << adjv << " " << adjw;
			}

			// 当前节点处理完毕后换行
			cout << "\n";
		}
	}
};

class Solution
{
public:
	int cost;
	double btime;
	vector<int> ptn; // partitioning
	vector<int> sc; // size_cluster
	//	vector<vector<int>>


	// 计算某个点v1的边权和 = 正边 + 负边
	int cal_pcost(const Graph &graph, const int &v1)
	{
		int pcost = 0;
		int size = int(graph.nodes[v1].edges.size());

		// 遍历所有邻接点
		for (int i = 0; i < size; i++)
		{
			int v2 = graph.nodes[v1].edges[i].first;
			int w = graph.nodes[v1].edges[i].second;

			if (ptn[v1] == ptn[v2] && w < 0) pcost += abs(w);
			if (ptn[v1] != ptn[v2] && w > 0) pcost += w;
		}

		return pcost;
	}


	// 暴力计算目标函数值，一般用于验证
	int cal_cost(const Graph &graph)
	{
		int ccost = 0;
		for (int i = 0; i < graph.nnode; i++)
		{
			int v1 = graph.nodes[i].idx;
			int size = int(graph.nodes[v1].edges.size());
			for (int j = 0; j < size; j++)
			{
				int v2 = graph.nodes[v1].edges[j].first;
				int w = graph.nodes[v1].edges[j].second;

				if (ptn[v1] == ptn[v2] && w < 0)
					ccost += abs(w);

				if (ptn[v1] != ptn[v2] && w > 0)
					ccost += w;
			}
		}

		ccost /= 2;
		return ccost;
	}


	// 构造一个空解
	Solution() : cost(MAX_VALUE), btime(0.0), ptn(), sc()
	{
	}


	// 构造一个随机初始解，并更新每个分区中的点数和解对应的cost
	Solution(const Graph &graph, clock_t cstime)
		: btime(0.0), ptn(graph.nnode, 0), sc(graph.k, 0)
	{
#if(DEBUG)
		printf("debug--1.1.1\n");fflush(stdout);
#endif
		// 给每个分区随机分配一个点
		int *randlist = new int[graph.nnode];
		Generate_Rand_List(randlist, graph.nnode); // 生成随机序列
#if(DEBUG)
		printf("debug--1.1.2\n");fflush(stdout);
#endif
		for (int pid = 0; pid < graph.k; pid++)
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
		for (int i = graph.k; i < graph.nnode; i++)
		{
			int v = randlist[i]; // 选点
			int p = rand() % graph.k; // 选分区
			ptn[v] = p; // 赋值
			sc[p]++; // 更新分区大小
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
		: cost(sol1.cost), btime(sol1.btime), ptn(sol1.ptn), sc(sol1.sc)
	{
	}


	// 复制另一个解
	void cpy(const Solution &sol1)
	{
		cost = sol1.cost;
		ptn = sol1.ptn; // 深拷贝，改变该数组不会影响到sol1.ptn
		sc = sol1.sc;
		btime = sol1.btime;
	}


	// 只验证cost
	bool verify(const Graph &graph)
	{
		bool is_true = false;

		int vcost = 0;
		for (int i = 0; i < graph.nnode; i++)
		{
			int v1 = graph.nodes[i].idx;
			int size = int(graph.nodes[v1].edges.size());
			for (int j = 0; j < size; j++)
			{
				int v2 = graph.nodes[v1].edges[j].first;
				int w = graph.nodes[v1].edges[j].second;

				if (ptn[v1] == ptn[v2] && w < 0)
					vcost += abs(w);

				if (ptn[v1] != ptn[v2] && w > 0)
					vcost += w;
			}
		}

		vcost /= 2;
		return vcost;

		if (vcost == cost)
			is_true = true;

		return is_true;
	}


	// 使用默认析构函数
	~Solution()
	{
	}
};

class Gain_node
{
public:
	int elem1; // 要移动的点
	int elem2; // 移动到的分区
	int elem3; // new pid
	int delta;
	int type;  // type = 0表示是移动，type = 1表示是交换

	Gain_node() : elem1(-1), elem2(-1), elem3(-1), delta(MAX_VALUE), type(-1)
	{
	}


	Gain_node(int e1, int e2, int e3, int delta, int type) : elem1(e1), elem2(e2), elem3(e3), delta(delta), type(type)
	{
	}


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


	~Gain_node()
	{
	}
};

//Graph graph;
clock_t start, end;
double *each_run_rlt;
double *each_run_time;
double *each_hit_time;
double avg_cost, avg_time, avg_htime, std_dev;
double sum_avg_cost = 0, sum_avg_time = 0, sum_avg_htime = 0;

int *asc_nodes1, *asc_nodes2, *asc_nodes3, *asc_nodes4;
int *desc_nodes1, *desc_nodes2, *desc_nodes3, *desc_nodes4;
vector<vector<int> > Descending_swap_v;
int **pos_gamma, **neg_gamma, **true_gamma;
int **active_matrix1, **active_matrix2;

// 记录最好的结果和最好的时间
int rbcost;
double rbtime;

//
char filename[1001] = "./instances/slashdot-zoo.graph"; //
int seed = 0;
double time_limit = 120;
int K = 64;
int runs = 1;
//int counter = 0; // Swap counters

// IMS算法参数
int max_nipv = 20; // IMS run length 10 (candidate values 50, 20)
double pct_sp = 0.2; // percentage of strong perturb  0.1  【感觉这个值不能太小】
double pct_wp = 0.01; // percentage of weak perturb   0.01
int sp, wp;
double param_q = 0.01; // 控制Swap的参数，biased_LS里面启动Swap的概率很小【= 1%】
double param_b = 0.8; // 控制ND-based_local_search和biased_local_search的参数【= 0.8效果还可以】


// ====================================================
// ==================== 函数声明========================
// ====================================================
// 主要算法函数
void Read_Parameters(int argc, char **argv);
void output_header(int run);
void allocate_memory(const Graph &graph);
void free_memory(const Graph &graph);

// 读入RH算法结果作为初始解
//void multistart_relocation_heuristic(const Graph &graph, Solution &bsol, double tlimit);
void read_RH_sol(const Graph &graph, Solution &csol, char *instancefile, const int &fno);

// 节点排序计算函数
void calculate_v1(const Graph &graph);
void calculate_v2(const Graph &graph);
void calculate_v3(const Graph &graph);
void calculate_v4(const Graph &graph, const Solution &csol);

// Gamma表相关函数
void init_Gamma(const Graph &graph, const Solution &csol);
void update_Gamma(const Graph &graph, Solution &csol, const Gain_node &gnode);

// Neighborhood decomposition matrix相关函数
void init_active_matrix(const Graph &graph);
void update_active_matrix(const Graph &graph, int pid1, int pid2);

// Swap操作相关函数
void cal_swap_value(const Graph &graph, const Solution &csol);
int get_weight(const Graph &graph, int vtx1, int vtx2);

// Local search算法
Solution local_search(const Graph &graph, const Solution &bsol, clock_t cstime);  // 最初始的版本
Solution biased_local_search1(const Graph &graph, const Solution &bsol, clock_t cstime);
Solution biased_local_search2(const Graph &graph, const Solution &bsol, clock_t cstime);
Solution biased_local_search3(const Graph &graph, const Solution &bsol, clock_t cstime);
Solution biased_local_search4(const Graph &graph, const Solution &bsol, clock_t cstime);
Solution swap_local_search(const Graph &graph, const Solution &bsol, clock_t cstime);
Solution local_search_decomposition(const Graph &graph, const Solution &bsol, clock_t cstime);
Solution swap_local_search_decomposition(const Graph &graph, const Solution &bsol, clock_t cstime);
Solution ND_based_LS(const Graph &graph, Solution &csol, clock_t cstime, double tlimit);
Solution Biased_LS(const Graph &graph, Solution &csol, clock_t cstime, double tlimit);

// Perturbation操作
Solution shake(const Graph &graph, Solution &csol, const int np);
Solution shake1(const Graph &graph, Solution &csol, const int np);

// 主算法
//Solution Maxima_search(const Graph &graph, Solution &csol, clock_t cstime, double tlimit);  // 最初始的版本
Solution rel_Maxima_search(const Graph &graph, Solution &csol, clock_t cstime, double tlimit);
Solution rel_Maxima_search1(const Graph &graph, Solution &csol, clock_t cstime, double tlimit);
Solution rel_Maxima_search2(const Graph &graph, Solution &csol, clock_t cstime, double tlimit);
Solution rel_Maxima_search3(const Graph &graph, Solution &csol, clock_t cstime, double tlimit);
void iterated_maxima_search(const Graph &graph, Solution &csol, double tlimit);

// 统计、验证和输出函数
void cal_indicators(int *each_run_rlt, int *each_run_time, const int &runs, int &bcost, double &avg_cost, double &avg_time);
void verify(const Graph &graph, Solution &bsol);
void write_IMSSDN_sol(const Graph &graph, const Solution &csol, char *instancefile, const int &fno);
void IMS_run(char *instancefile, double timelimit);
// ====================================================
// ====================================================
// ====================================================


/*
 * 读参数
 */
void Read_Parameters(int argc, char **argv)
{
	for(int i = 1; i < argc; i++)
	{
		cout << argv[i] << " ";
	}
	cout << endl;

	for (int i = 1; i < argc; i += 2) // [0]表示'-'，[1]表示类型，[2]表示数据
	{
		if (argv[i][0] != '-')
		{
			exit(0);
		}
		else if (argv[i][1] == 'i') // The file name
		{
			strncpy(filename, argv[i + 1], 1000);
//			cerr << filename << endl;
		}
		else if (argv[i][1] == 's') // seed
			seed = atoi(argv[i + 1]);
		else if (argv[i][1] == 'r') // The maximum time
			time_limit = atof(argv[i + 1]);
		else if (argv[i][1] == 'k') // K
			K = atoi(argv[i + 1]);
	}

//	check parameters
	if (strlen(filename) == 0)
	{
		cerr << "No input data" << endl;
		exit(1);
	}
	cout << "Parameter loaded successfully." << endl << endl;
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
	if (!opf)
	{
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

	pos_gamma = new int *[graph.k];
	neg_gamma = new int *[graph.k];
	true_gamma = new int *[graph.k];
	active_matrix1 = new int *[graph.k];
	active_matrix2 = new int *[graph.k];
	for (int i = 0; i < graph.k; i++)
	{
		pos_gamma[i] = new int[graph.nnode];
		neg_gamma[i] = new int[graph.nnode];
		true_gamma[i] = new int[graph.nnode];
		active_matrix1[i] = new int[graph.k];
		active_matrix2[i] = new int[graph.k];
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

	for (int i = 0; i < graph.k; i++)
	{
		delete[] pos_gamma[i];
		delete[] neg_gamma[i];
		delete[] true_gamma[i];
		delete[] active_matrix1[i];
		delete[] active_matrix2[i];
	}
	delete[] pos_gamma;
	delete[] neg_gamma;
	delete[] true_gamma;
	delete[] active_matrix1;
	delete[] active_matrix2;
}


///*
// * RH算法
// */
//void multistart_relocation_heuristic(const Graph &graph, Solution &bsol, double tlimit)
//{
//	// 开始重启迭代
//	int riter = 0; // restart_iter
//	clock_t cstime = clock();
//	while ((clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC) < tlimit)
//	{
//#if(DEBUG)
//		printf("debug--1.1\n");fflush(stdout);
//#endif
//		// 生成初始解
//		Solution csol = Solution(graph, cstime); // cur_sol
//		Solution tsol = Solution(csol); // temp_sol
//
//#if(DEBUG)
//		printf("debug--1.2\n");fflush(stdout);
//#endif
//		bool improved = true;
//		while (improved)
//		{
//			improved = false;
//#if(DEBUG)
//			printf("debug--1.3\n");fflush(stdout);
//#endif
//			for (int vtx1 = 0; vtx1 < graph.nnode; vtx1++)
//			{
////				int vtx1 = i;
//#if(DEBUG)
//				printf("debug--1.5\n");fflush(stdout);
//#endif
//				for (int clst2 = 0; clst2 < graph.k; clst2++)
//				{
//					int clst1 = csol.ptn[vtx1];
//					if (csol.sc[clst1] == 1 || clst2 == clst1)
//						continue;
//
//					tsol.ptn[vtx1] = clst2;
//					if (tsol.cal_pcost(graph, vtx1) < csol.cal_pcost(graph, vtx1))
//					{
//						csol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
//						csol.sc[clst1] -= 1;
//						csol.ptn[vtx1] = clst2;
//						csol.sc[clst2] += 1;
//						improved = true;
//					}
//					else
//					{
//						tsol.cpy(csol); // 则退回到当前解，并从当前解出发继续搜索
//					}
//				}
//			}
//		}
//
//		// 计算cost
//#if(DEBUG)
//		printf("debug--1.6\n");fflush(stdout);
//#endif
//		csol.cost = csol.cal_cost(graph);
//		if (csol.cost < rbcost)
//		{
//			rbcost = csol.cost;
//			rbtime = csol.btime;
//		}
//		if (csol.cost < bsol.cost)
//		{
////			printf("RH round %d, time=%.4f, ccost=%d, improve=%d\n", riter, (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC), csol.cost, bsol.cost - csol.cost);fflush(stdout);
//			bsol.cpy(csol);
//		}
//		riter++;
//#if(DPTN)
//		printf("RH round %d, time=%.4f, ccost=%d, rbcost=%d\n", riter, (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC), csol.cost, rbcost);fflush(stdout);
//#endif
//	}
//
////	printf("=========================== RH END ===========================\n");fflush(stdout);
//
//	// 写入文件
//#if(DEBUG)
//	printf("debug--1.7\n");fflush(stdout);
//#endif
//	char outputfile[1000];
//	char *graph_name = basename(filename);
//	sprintf(outputfile, "%soutput_dir/results_%s_%d.txt", root_path.c_str(), graph_name, K);
//	FILE *opf = fopen(outputfile, "a");
//	if (!opf)
//	{
//		perror("Failed to open RH output file");
//		exit(-1);
//	}
//	fprintf(opf, "RH restarts:%d, RH best cost=%d, RH best time=%.4f\n", riter, bsol.cost, bsol.btime);
//	fclose(opf);
//}


/*
 * 从文件中读进来初始解
 * 在C++中，将函数形参声明为const主要有以下作用：防止参数被修改并增强代码安全性。
 */
void read_RH_sol(const Graph &graph, Solution &csol, char *instancefile, const int &fno)
{
	// 读图
	char *graph_name = basename(instancefile); // 获取当前图的名称
	char rltfile[1000]; // 用于存储文件名
	if (tuning)
		sprintf(rltfile, "./init_sols/%s_%d_%d", graph_name, graph.k, fno);  // TODO ..和.来控制路径对错
	else
		sprintf(rltfile, "%sinit_sols/%s_%d_%d", root_path.c_str(), graph_name, graph.k, fno);
	// 打开解文件
	ifstream sol_file(rltfile);
	if (!sol_file.is_open())
	{
		cerr << "无法打开解文件: " << rltfile << endl;
//		exit(-20);
	}

	// 这里其实很凑巧的构造了一个初始不可行解
	csol.btime = 0.0;
	// 分配内存
	csol.ptn.resize(graph.nnode, 0);
	csol.sc.resize(graph.k, 0);

	// 读取cost
	sol_file >> csol.cost;

	// 读取solution
	// TODO 注意：解里面也保存了size，只是这里没读，可以作为测试
	// 这里有一个非常大的问题，假设文件没有正确打开，这里也可以构造一个不可行解
	// 其只有一个分区，点都在这个分区里面，但是从这个解开始能够跑到质量很好的解，这一点很神奇
	//【从侧面印证这个算法对初始解非常不敏感，这里测试一下】
	for (int i = 0; i < graph.nnode; i++)
	{
		sol_file >> csol.ptn[i];
		csol.sc[csol.ptn[i]]++;
	}

	// !更新rbtime
	if (csol.cost < rbcost)
	{
		rbcost = csol.cost;
		rbtime = (clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
	}

	sol_file.close();
}


/*
 * 按 节点的影响力=∑(|正边权|+|负边权|) 确定节点的遍历顺序
 */
void calculate_v1(const Graph &graph)
{
	// 重置
	for (int i = 0; i < graph.nnode; i++)
	{
		asc_nodes1[i] = i;
		desc_nodes1[i] = i;
	}

	// 计算节点影响力
	int *obj_value1 = new int[graph.nnode];
	int *obj_value2 = new int[graph.nnode];
	memset(obj_value1, 0, sizeof(int) * graph.nnode);

	for (int i = 0; i < graph.nnode; i++)
	{
		int v1 = graph.nodes[i].idx;
		int size = int(graph.nodes[i].edges.size());
		for (int j = 0; j < size; j++)
		{
			obj_value1[v1] += abs(graph.nodes[i].edges[j].second);
		}
	}
	memcpy(obj_value2, obj_value1, sizeof(int) * graph.nnode);

	// 排序
	Quick_Sort_asc(asc_nodes1, obj_value1, 0, graph.nnode - 1);
	Quick_Sort_desc(desc_nodes1, obj_value2, 0, graph.nnode - 1);

	// 释放内存
	delete[] obj_value1;
	delete[] obj_value2;
}


/*
 * 按 节点的影响力=∑|正边权| 确定节点的遍历顺序
 */
void calculate_v2(const Graph &graph)
{
	// 重置
	for (int i = 0; i < graph.nnode; i++)
	{
		asc_nodes2[i] = i;
		desc_nodes2[i] = i;
	}

	// 计算节点影响力
	int *obj_value1 = new int[graph.nnode];
	int *obj_value2 = new int[graph.nnode];
	memset(obj_value1, 0, sizeof(int) * graph.nnode);

	for (int i = 0; i < graph.nnode; i++)
	{
		int v1 = graph.nodes[i].idx;
		int size = int(graph.nodes[i].edges.size());
		for (int j = 0; j < size; j++)
		{
			if (graph.nodes[i].edges[j].second > 0)
				obj_value1[v1] += graph.nodes[i].edges[j].second;
		}
	}
	memcpy(obj_value2, obj_value1, sizeof(int) * graph.nnode);

	// 排序
	Quick_Sort_asc(asc_nodes2, obj_value1, 0, graph.nnode - 1);
	Quick_Sort_desc(desc_nodes2, obj_value2, 0, graph.nnode - 1);

	// 释放内存
	delete[] obj_value1;
	delete[] obj_value2;
}


/*
 * 按 节点的影响力=∑|负边权| 确定节点的遍历顺序
 */
void calculate_v3(const Graph &graph)
{
	// 重置
	for (int i = 0; i < graph.nnode; i++)
	{
		asc_nodes3[i] = i;
		desc_nodes3[i] = i;
	}

	// 计算节点影响力
	int *obj_value1 = new int[graph.nnode];
	int *obj_value2 = new int[graph.nnode];
	memset(obj_value1, 0, sizeof(int) * graph.nnode);

	for (int i = 0; i < graph.nnode; i++)
	{
		int v1 = graph.nodes[i].idx;
		assert(v1==i);
		int size = int(graph.nodes[i].edges.size());
		for (int j = 0; j < size; j++)
		{
			if (graph.nodes[i].edges[j].second < 0)
				obj_value1[v1] += abs(graph.nodes[i].edges[j].second); // 这里不要忘记abs()
		}
	}
	memcpy(obj_value2, obj_value1, sizeof(int) * graph.nnode);

	// 排序
	Quick_Sort_asc(asc_nodes3, obj_value1, 0, graph.nnode - 1);
	Quick_Sort_desc(desc_nodes3, obj_value2, 0, graph.nnode - 1);

	// 释放内存
	delete[] obj_value1;
	delete[] obj_value2;
}


/*
 * 根据当前解，按【节点的】 分区间的|正边权| + 分区里的|负边权| 之和确定节点的遍历顺序
 * 其实这里计算的就是每个节点对整体不平衡的贡献程度
 * 用于扰动的测试
 */
void calculate_v4(const Graph &graph, const Solution &csol)
{
	// 重置
	for (int i = 0; i < graph.nnode; i++)
	{
		asc_nodes4[i] = i;
		desc_nodes4[i] = i;
	}

	// 计算节点影响力
	int *obj_value1 = new int[graph.nnode];
	int *obj_value2 = new int[graph.nnode];
	memset(obj_value1, 0, sizeof(int) * graph.nnode);

	for (int i = 0; i < graph.nnode; i++)
	{
		int v1 = graph.nodes[i].idx;
		int size = int(graph.nodes[i].edges.size());
		for (int j = 0; j < size; j++)
		{
			if (graph.nodes[i].edges[j].second > 0 && csol.ptn[i] != csol.ptn[j])
				obj_value1[v1] += graph.nodes[i].edges[j].second;
			if (graph.nodes[i].edges[j].second < 0 && csol.ptn[i] == csol.ptn[j])
				obj_value1[v1] += abs(graph.nodes[i].edges[j].second);
		}
	}
	memcpy(obj_value2, obj_value1, sizeof(int) * graph.nnode);

	// 排序
	Quick_Sort_asc(asc_nodes4, obj_value1, 0, graph.nnode - 1);
	Quick_Sort_desc(desc_nodes4, obj_value2, 0, graph.nnode - 1);

	// 释放内存
	delete[] obj_value1;
	delete[] obj_value2;
}


/*
 * 根据csol初始化Gamma表
 */
void init_Gamma(const Graph &graph, const Solution &csol)
{
	for (int i = 0; i < graph.k; i++)
	{
		memset(neg_gamma[i], 0, sizeof(int) * graph.nnode);
		memset(pos_gamma[i], 0, sizeof(int) * graph.nnode);
		memset(true_gamma[i], 0, sizeof(int) * graph.nnode);
	}

	// 更新分开的gamma表
	for (int v1 = 0; v1 < graph.nnode; v1++)
	{
		int size = int(graph.nodes[v1].edges.size());
		for (int j = 0; j < size; j++)
		{
			int v2 = graph.nodes[v1].edges[j].first;
			int w = graph.nodes[v1].edges[j].second;
			int p = csol.ptn[v2];

			if (w < 0)
				neg_gamma[p][v1] += abs(w);
			else
				pos_gamma[p][v1] += w;
		}
	}

#if(is_trueGamma)
	// 更新总的gamma表，进一步加速
	for (int vtx1 = 0; vtx1 < graph.nnode; vtx1++)
	{
		int clst1 = csol.ptn[vtx1];
		for (int clst2 = 0; clst2 < graph.k; clst2++)
		{
			if (clst2 == csol.ptn[vtx1])
			{
				true_gamma[clst2][vtx1] = MAX_VALUE;
			}
			else
			{
				true_gamma[clst2][vtx1] = neg_gamma[clst2][vtx1] - neg_gamma[clst1][vtx1]
										+ pos_gamma[clst1][vtx1] - pos_gamma[clst2][vtx1];
			}
		}
	}
#endif
}


/*
 * 更新 Gamma 表，并移动
 */
void update_Gamma(const Graph &graph, Solution &csol, const Gain_node &gnode)
{
	if (gnode.type == 0)  // Add
	{
		int v1 = gnode.elem1;
		int p1 = csol.ptn[v1];
		int p2 = gnode.elem2;
		int size = graph.nodes[v1].edges.size();
		for (int j = 0; j < size; j++)
		{
			int v2 = graph.nodes[v1].edges[j].first;
			int w = graph.nodes[v1].edges[j].second;

			if (w < 0)
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

#if(is_trueGamma)
		/*
		 * 更新相关的true_gamma（涉及移动点和它的相邻点）
		 * 这里相当于是解已经更新好了，再来维护的true_gamma表，
		 * 所以是在移动后的解上直接计算和更新
		 */
		// 移动点的true_gamma
		int clst1 = csol.ptn[v1];
		for (int clst2 = 0; clst2 < graph.k; clst2++)
		{
			if (clst2 == csol.ptn[v1])
			{
				true_gamma[clst2][v1] = MAX_VALUE; // 把0换掉
			}
			else
			{
				true_gamma[clst2][v1] = neg_gamma[clst2][v1] - neg_gamma[clst1][v1]
									  + pos_gamma[clst1][v1] - pos_gamma[clst2][v1];
			}
		}

		// 移动点相邻点vn的true_gamma
		for (int i = 0; i < size; i++)
		{
			int vn = graph.nodes[v1].edges[i].first;
			int clst1 = csol.ptn[vn];
			for (int clst2 = 0; clst2 < graph.k; clst2++)
			{
				if (clst2 == csol.ptn[vn])
				{
					true_gamma[clst2][vn] = MAX_VALUE;
				}
				else
				{
					true_gamma[clst2][vn] = neg_gamma[clst2][vn] - neg_gamma[clst1][vn]
										  + pos_gamma[clst1][vn] - pos_gamma[clst2][vn];
				}
			}
		}
#endif
	}
	else if (gnode.type == 1)  // Swap
	{
		int v1 = gnode.elem1;
		int v2 = gnode.elem2;
		int p1 = csol.ptn[v1];
		int p2 = csol.ptn[v2];

		// v1 → p2
		int size = graph.nodes[v1].edges.size();
		for (int j = 0; j < size; j++)
		{
			int v2 = graph.nodes[v1].edges[j].first;
			int w = graph.nodes[v1].edges[j].second;

			if (w < 0)
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
		csol.ptn[v1] = p2;

#if(is_trueGamma)
		// 移动点v1的true_gamma
		int clst1 = csol.ptn[v1];
		for (int clst2 = 0; clst2 < graph.k; clst2++)
		{
			if (clst2 == csol.ptn[v1])
			{
				true_gamma[clst2][v1] = MAX_VALUE; // 把0换掉
			}
			else
			{
				true_gamma[clst2][v1] = neg_gamma[clst2][v1] - neg_gamma[clst1][v1]
									  + pos_gamma[clst1][v1] - pos_gamma[clst2][v1];
			}
		}

		// 移动点v1相邻点vn的true_gamma
		for (int i = 0; i < size; i++)
		{
			int vn = graph.nodes[v1].edges[i].first;
			int clst1 = csol.ptn[vn];
			for (int clst2 = 0; clst2 < graph.k; clst2++)
			{
				if (clst2 == csol.ptn[vn])
				{
					true_gamma[clst2][vn] = MAX_VALUE;
				}
				else
				{
					true_gamma[clst2][vn] = neg_gamma[clst2][vn] - neg_gamma[clst1][vn]
										  + pos_gamma[clst1][vn] - pos_gamma[clst2][vn];
				}
			}
		}
#endif

		// v2 → p1
		size = graph.nodes[v2].edges.size();
		for (int j = 0; j < size; j++)
		{
			int v1 = graph.nodes[v2].edges[j].first;
			int w = graph.nodes[v2].edges[j].second;

			if (w < 0)
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
		csol.ptn[v2] = p1;

#if(is_trueGamma)
		// 移动点v2的true_gamma
		clst1 = csol.ptn[v2];
		for (int clst2 = 0; clst2 < graph.k; clst2++)
		{
			if (clst2 == csol.ptn[v2])
			{
				true_gamma[clst2][v2] = MAX_VALUE; // 把0换掉
			}
			else
			{
				true_gamma[clst2][v2] = neg_gamma[clst2][v2] - neg_gamma[clst1][v2]
									  + pos_gamma[clst1][v2] - pos_gamma[clst2][v2];
			}
		}

		// 移动点v2相邻点vn的true_gamma
		for (int i = 0; i < size; i++)
		{
			int vn = graph.nodes[v2].edges[i].first;
			int clst1 = csol.ptn[vn];
			for (int clst2 = 0; clst2 < graph.k; clst2++)
			{
				if (clst2 == csol.ptn[vn])
				{
					true_gamma[clst2][vn] = MAX_VALUE;
				}
				else
				{
					true_gamma[clst2][vn] = neg_gamma[clst2][vn] - neg_gamma[clst1][vn]
										  + pos_gamma[clst1][vn] - pos_gamma[clst2][vn];
				}
			}
		}
#endif
		// 移动
		csol.cost += gnode.delta;
	}

#if(is_Verify)
	// 验证
	if (csol.cost != csol.cal_cost(graph))
	{
		printf("ccost=%d, vcost=%d", csol.cost, csol.cal_cost(graph));fflush(stdout);
		exit(-666);
		assert(csol.cost == csol.cal_cost(graph));
	}
#endif
}


/*
 * 初始化active_matrix全为1
 */
void init_active_matrix(const Graph &graph)
{
	for (int i = 0; i < graph.k; i++) {
		memset(active_matrix1[i], 1, sizeof(int) * graph.k);
		memset(active_matrix2[i], 1, sizeof(int) * graph.k);
		active_matrix1[i][i] = 0;
		active_matrix2[i][i] = 0;
	}
}


/*
 * 更新由移动引起的active_matrix的变化
 * 如果Add和Swap这两个邻域都用了的话，需要同时更新！
 */
void update_active_matrix(const Graph &graph, int pid1, int pid2)
{
	// 更新行
	memset(active_matrix1[pid1], 1, sizeof(int) * graph.k);
	memset(active_matrix1[pid2], 1, sizeof(int) * graph.k);
	memset(active_matrix2[pid1], 1, sizeof(int) * graph.k);
	memset(active_matrix2[pid2], 1, sizeof(int) * graph.k);

	// 更新列
	for (int i = 0; i < graph.k; i++)
	{
		active_matrix1[i][pid1] = 1;
		active_matrix1[i][pid2] = 1;
		active_matrix2[i][pid1] = 1;
		active_matrix2[i][pid2] = 1;
	}
}


/*
 * 找到点x和点y的权重，因为是稀疏图，时间复杂度 << O(n)
 */
int get_weight(const Graph &graph, int vtx1, int vtx2)
{
	int weight = 0;
	// 因为是无向图，所以遍历邻居节点多的那个点，找到它们之间边的权重
	if (graph.nodes[vtx1].edges.size() <= graph.nodes[vtx2].edges.size())
	{
		int nynbh = int(graph.nodes[vtx2].edges.size());
		for (int i = 0; i < nynbh; i++)
		{
			if (vtx1 == graph.nodes[vtx2].edges[i].first)
				weight = graph.nodes[vtx2].edges[i].second;
		}
	}
	else
	{
		int nxnbh = int(graph.nodes[vtx1].edges.size());
		for (int i = 0; i < nxnbh; i++)
		{
			if (vtx2 == graph.nodes[vtx1].edges[i].first)
				weight = graph.nodes[vtx1].edges[i].second;
		}
	}

	return weight;
}


///*
// * 计算交换的权重，并更新Descending_swap_v中的值
// * @param graph
// * @param csol
// */
//void cal_swap_value(const Graph &graph, const Solution &csol)
//{
//	int swap_n = int(graph.nnode / graph.k);
//	vector<vector<pair<int, int> > > swap_v(graph.k);
//	for (int vtx1 = 0; vtx1 < graph.nnode; vtx1++)
//	{
//		int v_cost = 0;
//		int pid1 = csol.ptn[vtx1];
//		for (int pid2 = 0; pid2 < graph.k; pid2++)
//		{
//			if (pid2 != pid1) // 如果分区不同
//			{
//				v_cost += pos_gamma[pid2][vtx1]; // 加上正的边权
//				v_cost -= neg_gamma[pid2][vtx1]; // 减去负的边权
//			}
//			else // 如果分区相同
//			{
//				v_cost -= pos_gamma[pid2][vtx1]; // 减去正的边权
//				v_cost += neg_gamma[pid2][vtx1]; // 加上负的边权
//			}
//		}
//		swap_v[pid1].push_back(make_pair(v_cost, vtx1));
//	}

//	for (int pid = 0; pid < graph.k; pid++)
//	{
//		// 每个分区中的点按照v_cost的【降序排序】
//		sort(swap_v[pid].begin(), swap_v[pid].end(), [](const pair<int, int> &a, const pair<int, double> &b)
//		{
//			return a.first > b.first;
//		});
////		int len_dsvk = int(Descending_swap_v[pid].size()); // pid分区中点的个数，length of descending swap [v] to [k]
//		int len_dsvk = int(swap_v[pid].size()); // pid分区中点的个数，length of descending swap [v] to [k]
//		if (len_dsvk >= swap_n)
//		{
//			for (int i = 0; i < swap_n; i++)
//				Descending_swap_v[pid].push_back(swap_v[pid][i].second);
//		}
//		else
//		{
//			for (int i = 0; i < len_dsvk; i++)
//				Descending_swap_v[pid].push_back(swap_v[pid][i].second);
//		}
//	}
//}

/**
 * 计算交换的权重，并更新Descending_swap_v中的值
 * @param graph
 * @param csol
 */
void cal_swap_value(const Graph &graph, const Solution &csol)
{
	// 重置数组，这一步很重要
	for (vector<int> &vec: Descending_swap_v)
	{
		vec.clear(); // 清空每个分区里的元素
	}
	int swap_n = int(graph.nnode / graph.k); // 平均每个分区中的点数
	vector<vector<pair<int, int> > > swap_v(graph.k); // pid: 分区，<权重和，该分区中的点>
	for (int vtx1 = 0; vtx1 < graph.nnode; vtx1++)
	{
		int v_cost = 0;
		int pid1 = csol.ptn[vtx1];
		for (int pid2 = 0; pid2 < graph.k; pid2++)
		{
			if (pid2 != pid1) // 如果分区不同
			{
				v_cost += pos_gamma[pid2][vtx1]; // 加正的边权
				v_cost -= neg_gamma[pid2][vtx1]; // 减去负的边权
			}
			else // 如果分区相同
			{
				v_cost -= pos_gamma[pid2][vtx1]; // 减去正的边权
				v_cost += neg_gamma[pid2][vtx1]; // 加上负的边权
			}
		}
		swap_v[pid1].push_back(make_pair(v_cost, vtx1));
	}

	/*
	 * swap_n是平均每个分区中应该分到的点数，这样存点是为了防止某个分区中的点太多，对这个分区的评估开销更多。
	 * 所以存储待交换的节点时，只取了前面几个比较有价值的交换点，也节省一部分时间。
	 * 这种做法也有点类似于之前说的【重点节点】和【从属节点】
	 */
	for (int pid = 0; pid < graph.k; pid++)
	{
		// 对swap_v中每个分区的点按照权重值进行【降序排序】
		sort(swap_v[pid].begin(), swap_v[pid].end(), [](const pair<int, int> &a, const pair<int, double> &b)
		{
			return a.first > b.first;
		});
		int len_dsvk = int(swap_v[pid].size()); // length of descending swap [v] to [k]
		// 这里存比较点数比较【少】的部分
		if (len_dsvk >= swap_n)
		{
			for (int i = 0; i < swap_n; i++)
				Descending_swap_v[pid].push_back(swap_v[pid][i].second);
		}
		else // len_dsvk < swap_n
		{
			for (int i = 0; i < len_dsvk; i++)
				Descending_swap_v[pid].push_back(swap_v[pid][i].second);
		}
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
	while (improved)
	{
		improved = false;
#if(DEBUG)
		printf("debug--1.3\n");fflush(stdout);
#endif
		for (int i = 0; i < graph.nnode; i++)
		{
			int vtx1 = i;
			int clst1 = csol.ptn[vtx1];
			int size = int(graph.nodes[i].edges.size());
			if (csol.sc[clst1] == 1 || size == 0)
				continue;
#if(DEBUG)
			printf("debug--1.5\n");fflush(stdout);
#endif
			for (int clst2 = 0; clst2 < graph.k; clst2++)
			{
				if (clst2 == clst1) continue;
				if (csol.sc[clst1] == 1) break;
#if(is_trueGamma)
				int delta = true_gamma[clst2][vtx1];
#else
				int delta = neg_gamma[clst2][vtx1] - neg_gamma[clst1][vtx1]
						  + pos_gamma[clst1][vtx1] - pos_gamma[clst2][vtx1];
#endif
				if (delta < 0)
				{
					Gain_node gn = Gain_node(vtx1, clst2, -1, delta, 0);
					update_Gamma(graph, csol, gn);
					improved = true;
#if(is_Verify)
					assert(csol.cost == csol.cal_cost(graph));
#endif
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
	while (improved)
	{
		improved = false;
#if(DEBUG)
		printf("debug--1.3\n");fflush(stdout);
#endif
		for (int i = 0; i < graph.nnode; i++)
		{
			int vtx1 = asc_nodes1[i];
			int clst1 = csol.ptn[vtx1];
			int size = int(graph.nodes[i].edges.size());
			if (csol.sc[clst1] == 1 || size == 0)
				continue;
#if(DEBUG)
			printf("debug--1.5\n");fflush(stdout);
#endif
			for (int clst2 = 0; clst2 < graph.k; clst2++)
			{
				if (clst2 == clst1) continue;
				if (csol.sc[clst1] == 1) break;
#if(is_trueGamma)
				int delta = true_gamma[clst2][vtx1];
#else
				int delta = neg_gamma[clst2][vtx1] - neg_gamma[clst1][vtx1]
						  + pos_gamma[clst1][vtx1] - pos_gamma[clst2][vtx1];
#endif
				if (delta < 0)
				{
					Gain_node gn = Gain_node(vtx1, clst2, -1, delta, 0);
					update_Gamma(graph, csol, gn);
#if(is_Verify)
					assert(csol.cost == csol.cal_cost(graph));
#endif
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
	while (improved)
	{
		improved = false;
#if(DEBUG)
		printf("debug--1.3\n");fflush(stdout);
#endif
		for (int i = 0; i < graph.nnode; i++)
		{
			int vtx1 = asc_nodes2[i];
			int clst1 = csol.ptn[vtx1];
			int size = int(graph.nodes[i].edges.size());
			if (csol.sc[clst1] == 1 || size == 0)
				continue;
#if(DEBUG)
			printf("debug--1.5\n");fflush(stdout);
#endif
			for (int clst2 = 0; clst2 < graph.k; clst2++)
			{
				if (clst2 == clst1) continue;
				if (csol.sc[clst1] == 1) break;
#if(is_trueGamma)
				int delta = true_gamma[clst2][vtx1];
#else
				int delta = neg_gamma[clst2][vtx1] - neg_gamma[clst1][vtx1]
						  + pos_gamma[clst1][vtx1] - pos_gamma[clst2][vtx1];
#endif
				if (delta < 0)
				{
					Gain_node gn = Gain_node(vtx1, clst2, -1, delta, 0);
					update_Gamma(graph, csol, gn);
#if(is_Verify)
					assert(csol.cost == csol.cal_cost(graph));
#endif
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
	while (improved)
	{
		improved = false;
#if(DEBUG)
		printf("debug--1.3\n");fflush(stdout);
#endif
		for (int i = 0; i < graph.nnode; i++)
		{
			int vtx1 = asc_nodes3[i];
			int clst1 = csol.ptn[vtx1];
			int size = int(graph.nodes[i].edges.size());
			if (csol.sc[clst1] == 1 || size == 0)
				continue;
#if(DEBUG)
			printf("debug--1.5\n");fflush(stdout);
#endif
			for (int clst2 = 0; clst2 < graph.k; clst2++)
			{
				if (clst2 == clst1) continue;
				if (csol.sc[clst1] == 1) break;
#if(is_trueGamma)
				int delta = true_gamma[clst2][vtx1];
#else
				int delta = neg_gamma[clst2][vtx1] - neg_gamma[clst1][vtx1]
				          + pos_gamma[clst1][vtx1] - pos_gamma[clst2][vtx1];
#endif

				if (delta < 0)
				{
					Gain_node gn = Gain_node(vtx1, clst2, -1, delta, 0);
					update_Gamma(graph, csol, gn);
#if(is_Verify)
					assert(csol.cost == csol.cal_cost(graph));
#endif
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
 * 使用gamma表进行快速更新，使用遍历顺序4（动态更新）
 * 这种有个问题：是在搜索之前更新节点权重（目前是这样），还是在每次更新之后更新节点权重？
 */
Solution biased_local_search4(const Graph &graph, const Solution &bsol, clock_t cstime)
{
//	printf("========================== BLS4 BEGIN ==========================\n");fflush(stdout);
	calculate_v4(graph, bsol);
	Solution csol = Solution(bsol);
	init_Gamma(graph, csol);
#if(DEBUG)
	printf("debug--1.2\n");fflush(stdout);
#endif
	bool improved = true;
	int iter = 0;
	while (improved)
	{
		improved = false;
#if(DEBUG)
		printf("debug--1.3\n");fflush(stdout);
#endif
		for (int i = 0; i < graph.nnode; i++)
		{
			int vtx1 = asc_nodes4[i];
			int clst1 = csol.ptn[vtx1];
			int size = int(graph.nodes[i].edges.size());
			if (csol.sc[clst1] == 1 || size == 0)
				continue;
#if(DEBUG)
			printf("debug--1.5\n");fflush(stdout);
#endif
			for (int clst2 = 0; clst2 < graph.k; clst2++)
			{
				if (clst2 == clst1) continue;
				if (csol.sc[clst1] == 1) break;
#if(is_trueGamma)
				int delta = true_gamma[clst2][vtx1];
#else
				int delta = neg_gamma[clst2][vtx1] - neg_gamma[clst1][vtx1]
					      + pos_gamma[clst1][vtx1] - pos_gamma[clst2][vtx1];
#endif
				if (delta < 0)
				{
					Gain_node gn = Gain_node(vtx1, clst2, -1, delta, 0);
					update_Gamma(graph, csol, gn);
#if(is_Verify)
					assert(csol.cost == csol.cal_cost(graph));
#endif
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


///*
// * 使用gamma表进行快速更新的普通local search
// */
//Solution local_search_best(const Graph &graph, const Solution &bsol, clock_t cstime)
//{
////	printf("========================== LS BEGIN ==========================\n");fflush(stdout);
//
//	Solution csol = Solution(bsol);
//	init_Gamma(graph, csol);
//#if(DEBUG)
//	printf("debug--1.2\n");fflush(stdout);
//#endif
//	bool improved = true;
//	int iter = 0;
//	while (improved)
//	{
//		improved = false;
//		int best_delta = INT_MAX; // 最好的delta
//		Gain_node best_gn; // 最好的邻域交换结构
//#if(DEBUG)
//		printf("debug--1.3\n");fflush(stdout);
//#endif
//		for (int i = 0; i < graph.nnode; i++)
//		{
//			int vtx1 = i;
//			int clst1 = csol.ptn[vtx1];
////			int size = int(graph.nodes[i].edges.size());
////			if (csol.sc[clst1] == 1 || size == 0) // size = 0表示该点为孤立点
////				continue;
//#if(DEBUG)
//			printf("debug--1.5\n");fflush(stdout);
//#endif
//			for (int clst2 = 0; clst2 < graph.k; clst2++)
//			{
//				if (csol.sc[clst1] == 1 || clst2 == clst1)
//					continue;
//
//#if(is_trueGamma)
//				int delta = true_gamma[clst2][vtx1];
//#else
//				int delta = neg_gamma[clst2][vtx1] - neg_gamma[clst1][vtx1]
//						  - pos_gamma[clst2][vtx1] + pos_gamma[clst1][vtx1];
//#endif
//
//#if(is_bestImprove)
//				if (delta < best_delta)
//				{
//					best_delta = delta;
//					best_gn = Gain_node(vtx1, clst2, -1, delta, 0);
//				}
//#else
//				if (delta < 0)
//				{
//					best_delta = delta;
//					best_gn = Gain_node(vtx1, clst2, -1, delta, 0);
//					goto update;
//				}
//#endif
//			}
//		}
//
//update:
//		if (best_delta < 0)
//		{
////			Gain_node gn = Gain_node(vtx1, clst2, -1, delta, 0);
//			update_Gamma(graph, csol, best_gn);
////
////			if (csol.cost < rbcost)
////			{
////				rbcost = csol.cost;
////				rbtime = (clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
////			}
////			csol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
////			assert(csol.cost == csol.cal_cost(graph));
//			improved = true;
//		}
//		iter++;
//#if(DPTN)
//		printf("LS1 round %d, time=%.4f, ccost=%d, rbcost=%d\n", iter, (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC), csol.cost, rbcost);fflush(stdout);
//#endif
//	}
//
////	printf("=========================== LS END ===========================\n");fflush(stdout);
//
//	return csol;
//}
//
//
///*
// * 使用gamma表进行快速更新，使用遍历顺序1：节点的影响力=∑(|正边权|+|负边权|) 的【升序】遍历节点
// */
//Solution biased_local_search1_best(const Graph &graph, const Solution &bsol, clock_t cstime)
//{
////	printf("========================== BLS1 BEGIN ==========================\n");fflush(stdout);
//
//	Solution csol = Solution(bsol);
//	init_Gamma(graph, csol);
//#if(DEBUG)
//	printf("debug--1.2\n");fflush(stdout);
//#endif
//	bool improved = true;
//	int iter = 0;
//	while (improved)
//	{
//		improved = false;
//		int best_delta = INT_MAX; // 最好的delta
//		Gain_node best_gn; // 最好的邻域交换结构
//#if(DEBUG)
//		printf("debug--1.3\n");fflush(stdout);
//#endif
//		for (int i = 0; i < graph.nnode; i++)
//		{
//			int vtx1 = asc_nodes1[i];
//			int clst1 = csol.ptn[vtx1];
////			int size = int(graph.nodes[i].edges.size());
////			if (csol.sc[clst1] == 1 || size == 0) // size = 0表明是孤立点
////				continue;
//#if(DEBUG)
//			printf("debug--1.5\n");fflush(stdout);
//#endif
//			for (int clst2 = 0; clst2 < graph.k; clst2++)
//			{
//				clst1 = csol.ptn[vtx1];
//				if (csol.sc[clst1] == 1 || clst2 == clst1)
//					continue;
//
//#if(is_trueGamma)
//				int delta = true_gamma[clst2][vtx1];
//#else
//				int delta = neg_gamma[clst2][vtx1] - neg_gamma[clst1][vtx1]
//						  - pos_gamma[clst2][vtx1] + pos_gamma[clst1][vtx1];
//#endif
//
//#if(is_bestImprove)
//				if (delta < best_delta)
//				{
//					best_delta = delta;
//					best_gn = Gain_node(vtx1, clst2, -1, delta, 0);
//				}
//#else
//				if (delta < 0)
//				{
//					best_delta = delta;
//					best_gn = Gain_node(vtx1, clst2, -1, delta, 0);
//					goto update;
//				}
//#endif
//			}
//		}
//
//update:
//		if (best_delta < 0)
//		{
////			Gain_node gn = Gain_node(vtx1, clst2, -1, delta, 0);
//			update_Gamma(graph, csol, best_gn);
//
////			if (csol.cost < rbcost)
////			{
////				rbcost = csol.cost;
////				rbtime = (clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
////			}
////			csol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
////			assert(csol.cost == csol.cal_cost(graph));
//			improved = true;
//
//		}
//		iter++;
//#if(DPTN)
//		printf("BLS1 round %d, time=%.4f, ccost=%d, rbcost=%d\n", iter, (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC), csol.cost, rbcost);fflush(stdout);
//#endif
//	}
//
////	printf("=========================== BLS1 END ===========================\n");fflush(stdout);
//
//	return csol;
//}
//
//
///*
// * 使用gamma表进行快速更新，使用遍历顺序2：节点的影响力=∑|正边权| 的【升序】遍历节点
// */
//Solution biased_local_search2_best(const Graph &graph, const Solution &bsol, clock_t cstime)
//{
////	printf("========================== BLS2 BEGIN ==========================\n");fflush(stdout);
//
//	Solution csol = Solution(bsol);
//	init_Gamma(graph, csol);
//#if(DEBUG)
//	printf("debug--1.2\n");fflush(stdout);
//#endif
//	bool improved = true;
//	int iter = 0;
//	while (improved)
//	{
//		improved = false;
//		int best_delta = INT_MAX; // 最好的delta
//		Gain_node best_gn; // 最好的邻域交换结构
//#if(DEBUG)
//		printf("debug--1.3\n");fflush(stdout);
//#endif
//		for (int i = 0; i < graph.nnode; i++)
//		{
//			int vtx1 = asc_nodes2[i];
//			int clst1 = csol.ptn[vtx1];
////			int size = int(graph.nodes[i].edges.size());
////			if (csol.sc[clst1] == 1 || size == 0)
////				continue;
//#if(DEBUG)
//			printf("debug--1.5\n");fflush(stdout);
//#endif
//			for (int clst2 = 0; clst2 < graph.k; clst2++)
//			{
//				clst1 = csol.ptn[vtx1];
//				if (csol.sc[clst1] == 1 || clst2 == clst1)
//					continue;
//
//#if(is_trueGamma)
//				int delta = true_gamma[clst2][vtx1];
//#else
//				int delta = neg_gamma[clst2][vtx1] - neg_gamma[clst1][vtx1]
//						  - pos_gamma[clst2][vtx1] + pos_gamma[clst1][vtx1];
//#endif
//
//#if(is_bestImprove)
//				if (delta < best_delta)
//				{
//					best_delta = delta;
//					best_gn = Gain_node(vtx1, clst2, -1, delta, 0);
//				}
//#else
//				if (delta < 0)
//				{
//					best_delta = delta;
//					best_gn = Gain_node(vtx1, clst2, -1, delta, 0);
//					goto update;
//				}
//#endif
//			}
//		}
//
//update:
//		if (best_delta < 0)
//		{
////			Gain_node gn = Gain_node(vtx1, clst2, -1, delta, 0);
//			update_Gamma(graph, csol, best_gn);
//
////			if (csol.cost < rbcost)
////			{
////				rbcost = csol.cost;
////				rbtime = (clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
////			}
////			csol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
////			assert(csol.cost == csol.cal_cost(graph));
//			improved = true;
//
//		}
//		iter++;
//#if(DPTN)
//		printf("BLS2 round %d, time=%.4f, ccost=%d, rbcost=%d\n", iter, (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC), csol.cost, rbcost);fflush(stdout);
//#endif
//	}
//
////	printf("=========================== BLS2 END ===========================\n");fflush(stdout);
//
//	return csol;
//}
//
//
///*
// * 使用gamma表进行快速更新，使用遍历顺序3：节点的影响力=∑|负边权| 的【升序】遍历节点
// */
//Solution biased_local_search3_best(const Graph &graph, const Solution &bsol, clock_t cstime)
//{
//	//	printf("========================== BLS3 BEGIN ==========================\n");fflush(stdout);
//
//	Solution csol = Solution(bsol);
//	init_Gamma(graph, csol);
//#if(DEBUG)
//	printf("debug--1.2\n");fflush(stdout);
//#endif
//	bool improved = true;
//	int iter = 0;
//	while (improved)
//	{
//		improved = false;
//		int best_delta = INT_MAX; // 最好的delta
//		Gain_node best_gn; // 最好的邻域交换结构
//#if(DEBUG)
//		printf("debug--1.3\n");fflush(stdout);
//#endif
//
//		for (int i = 0; i < graph.nnode; i++)
//		{
//			int vtx1 = asc_nodes3[i];
//			int clst1 = csol.ptn[vtx1];
////			int size = int(graph.nodes[i].edges.size());
////			if (csol.sc[clst1] == 1 || size == 0)
////				continue;
//#if(DEBUG)
//			printf("debug--1.5\n");fflush(stdout);
//#endif
//			for (int clst2 = 0; clst2 < graph.k; clst2++)
//			{
//				clst1 = csol.ptn[vtx1];
//				if (csol.sc[clst1] == 1 || clst2 == clst1)
//					continue;
//
//#if(is_trueGamma)
//				int delta = true_gamma[clst2][vtx1];
//#else
//				int delta = neg_gamma[clst2][vtx1] - neg_gamma[clst1][vtx1]
//						  - pos_gamma[clst2][vtx1] + pos_gamma[clst1][vtx1];
//#endif
//
//#if(is_bestImprove)
//				if (delta < best_delta)
//				{
//					best_delta = delta;
//					best_gn = Gain_node(vtx1, clst2, -1, delta, 0);
//				}
//#else
//				if (delta < 0)
//				{
//					best_delta = delta;
//					best_gn = Gain_node(vtx1, clst2, -1, delta, 0);
//					goto update;
//				}
//#endif
//			}
//		}
//
//update:
//		if (best_delta < 0)
//		{
////			Gain_node gn = Gain_node(vtx1, clst2, -1, delta, 0);
//			update_Gamma(graph, csol, best_gn);
//
////			if (csol.cost < rbcost)
////			{
////				rbcost = csol.cost;
////				rbtime = (clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
////			}
////			csol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
////			assert(csol.cost == csol.cal_cost(graph));
//			improved = true;
//
//		}
//		iter++;
//#if(DPTN)
//		printf("BLS3 round %d, time=%.4f, ccost=%d, rbcost=%d\n", iter, (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC), csol.cost, rbcost);fflush(stdout);
//#endif
//	}
//
//	return csol;
//}
//
//
///*
// * 使用gamma表进行快速更新，使用遍历顺序4（动态更新）
// * 这种有个问题：是在搜索之前更新节点权重（目前是这样），还是在每次更新之后更新节点权重？
// */
//Solution biased_local_search4_best(const Graph &graph, const Solution &bsol, clock_t cstime)
//{
////	printf("========================== BLS3 BEGIN ==========================\n");fflush(stdout);
//	calculate_v4(graph, bsol);
//	Solution csol = Solution(bsol);
//	init_Gamma(graph, csol);
//#if(DEBUG)
//	printf("debug--1.2\n");fflush(stdout);
//#endif
//	bool improved = true;
//	int iter = 0;
//	while (improved)
//	{
//		improved = false;
//		int best_delta = INT_MAX; // 最好的delta
//		Gain_node best_gn; // 最好的邻域交换结构
//#if(DEBUG)
//		printf("debug--1.3\n");fflush(stdout);
//#endif
//
//		for (int i = 0; i < graph.nnode; i++)
//		{
//			int vtx1 = asc_nodes4[i];
//			int clst1 = csol.ptn[vtx1];
////			int size = int(graph.nodes[i].edges.size());
////			if (csol.sc[clst1] == 1 || size == 0)
////				continue;
//#if(DEBUG)
//			printf("debug--1.5\n");fflush(stdout);
//#endif
//			for (int clst2 = 0; clst2 < graph.k; clst2++)
//			{
//				clst1 = csol.ptn[vtx1];
//				if (csol.sc[clst1] == 1 || clst2 == clst1)
//					continue;
//
//#if(is_trueGamma)
//				int delta = true_gamma[clst2][vtx1];
//#else
//				int delta = neg_gamma[clst2][vtx1] - neg_gamma[clst1][vtx1]
//						  - pos_gamma[clst2][vtx1] + pos_gamma[clst1][vtx1];
//#endif
//
//#if(is_bestImprove)
//				if (delta < best_delta)
//				{
//					best_delta = delta;
//					best_gn = Gain_node(vtx1, clst2, -1, delta, 0);
//				}
//#else
//				if (delta < 0)
//				{
//					best_delta = delta;
//					best_gn = Gain_node(vtx1, clst2, -1, delta, 0);
//					goto update;
//				}
//#endif
//			}
//		}
//
//update:
//		if (best_delta < 0)
//		{
////			Gain_node gn = Gain_node(vtx1, clst2, -1, delta, 0);
//			update_Gamma(graph, csol, best_gn);
//
////			if (csol.cost < rbcost)
////			{
////				rbcost = csol.cost;
////				rbtime = (clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
////			}
////			csol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
////			assert(csol.cost == csol.cal_cost(graph));
//			improved = true;
//
//		}
//		iter++;
//#if(DPTN)
//		printf("BLS4 round %d, time=%.4f, ccost=%d, rbcost=%d\n", iter, (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC), csol.cost, rbcost);fflush(stdout);
//#endif
//	}
//
//	return csol;
//}


// Function to perform swap local search
Solution swap_local_search(const Graph &graph, const Solution &bsol, clock_t cstime)
{
	Solution csol = Solution(bsol);
	init_Gamma(graph, csol);
#if(DEBUG)
	printf("debug--1.2\n");fflush(stdout);
#endif
	bool improved = true;
	int iter = 0;
	while (improved)
	{
		improved = false;
#if(DEBUG)
		printf("debug--1.3\n");fflush(stdout);
#endif
		for (int pid1 = 0; pid1 < graph.k; pid1++)
		{
			for (int pid2 = pid1 + 1; pid2 < graph.k; pid2++)
			{
				/* const auto是C++11引入的语法特性，用于在声明变量时自动推导类型并设置为常量。
				 * 其核心作用是将变量声明为不可修改的常量，同时自动推导其类型。
				 */
				for (const auto &vtx1: Descending_swap_v[pid1])
				{
					for (const auto &vtx2: Descending_swap_v[pid2])
					{
						if (csol.ptn[vtx1] == csol.ptn[vtx2]) continue;
						int weight = get_weight(graph, vtx1, vtx2);
#if(is_trueGamma)
						int delta = 2 * weight
								  + true_gamma[pid2][vtx1] + true_gamma[pid1][vtx2];
#else
						int delta = 2 * weight
								  + neg_gamma[pid2][vtx1] - neg_gamma[pid1][vtx1]
								  + pos_gamma[pid1][vtx1] - pos_gamma[pid2][vtx1]
								  + neg_gamma[pid1][vtx2] - neg_gamma[pid2][vtx2]
								  + pos_gamma[pid2][vtx2] - pos_gamma[pid1][vtx2];
#endif
						if (delta < 0)
						{
							Gain_node gn = Gain_node(vtx1, vtx2, -1, delta, 1);
							update_Gamma(graph, csol, gn);
#if(is_Verify)
							assert(csol.cost == csol.cal_cost(graph));
#endif
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


/**
 * 有邻域分解的基于Add的局部搜索
 * @param graph
 * @param bsol
 * @param cstime
 * @return
 */
Solution local_search_decomposition(const Graph &graph, const Solution &bsol, clock_t cstime)
{
//	printf("========================== LSD BEGIN ==========================\n");fflush(stdout);
	Solution csol = Solution(bsol);
	init_Gamma(graph, csol);

	// 获取每个分区中的元素序列
	/* unordered_set是C++标准库STL中的一种无序集合容器，底层采用哈希表实现，支持快速查找、插入和删除操作，平均时间复杂度为O(1)。
	 * 它仅存储唯一元素，不保证元素顺序，适用于需要快速判断元素是否存在但不关心排序的场景。*/
	unordered_set<int> clst_elems[graph.k]; // 创建 k 个 unordered_set<int>，保存每个分区中的点
	for (int i = 0; i < graph.nnode; i++) // 根据划分结果，将节点编号添加到对应的 unordered_set<int> 中
	{
		clst_elems[csol.ptn[i]].insert(i);
	}
#if(DEBUG)
	printf("debug--1.2\n");fflush(stdout);
#endif
	bool improved = true, flag;
	int iter = 0;
 	while (improved)
	{
		improved = false;
#if(DEBUG)
		printf("debug--1.3\n");fflush(stdout);
#endif
		// 遍历分区1
		for (int clst1 = 0; clst1 < graph.k; clst1++)
		{
			// 检查是否可以移除点
			if (csol.sc[clst1] <= 1) continue;

			// 遍历分区2，评估分区1中的每个点到分区2的delta值
			for (int clst2 = 0; clst2 < graph.k; clst2++)
			{
				if (csol.sc[clst1] <= 1) break;
				if (clst2 == clst1) continue;
				if (active_matrix1[clst1][clst2] == 0) continue;
				active_matrix1[clst1][clst2] = 0;
				flag = false;

				// 将待遍历的点存储起来！！
				// 因为clst_elems一直在变动，然后它又是无序的，然后我们又一直在遍历其中的点，所以在服务器上会报断错误
				// 但为什么会报段错误我就不明白了
				vector<int> verticesArray;
				for (int vtx1: clst_elems[clst1])
				{
					verticesArray.push_back(vtx1);
				}

				// 遍历分区1中的点
//				for (int vtx1: clst_elems[clst1])
				for (int vtx1: verticesArray)
				{
					if (csol.sc[clst1] <= 1) break;

#if(is_trueGamma)
					int delta = true_gamma[clst2][vtx1];
//					printf("clst1=%d, clst2=%d, delta=%d\n", clst1, clst2, delta);
#else
					int delta = neg_gamma[clst2][vtx1] - neg_gamma[clst1][vtx1]
							  + pos_gamma[clst1][vtx1] - pos_gamma[clst2][vtx1];
#endif
					if (delta < 0)
					{
						Gain_node gn = Gain_node(vtx1, clst2, -1, delta, 0);
						int bvtx1 = gn.elem1;
						int bclst1 = csol.ptn[bvtx1];
						int bclst2 = gn.elem2;

						update_Gamma(graph, csol, gn);
						// 更新分区中的点队列
						clst_elems[bclst1].erase(bvtx1);
						clst_elems[bclst2].insert(bvtx1);

						improved = true;
						flag = true;
#if(is_Verify)
						assert(csol.sc[clst1] > 0 && csol.sc[clst2] > 0);
						assert(csol.cost == csol.cal_cost(graph));
#endif
					}
				}

				// 更新 Gamma 表
				if (flag == true)
				{
					update_active_matrix(graph, clst1, clst2);
				}
			}
			iter++;
#if(DPTN)
			printf("LS1 round %d, time=%.4f, ccost=%d, rbcost=%d\n",
					iter, (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC), csol.cost, rbcost);fflush(stdout);
#endif
		}
	}
//	printf("=========================== LSD END ===========================\n");fflush(stdout);

	return csol;
}


/**
 * 有邻域分解的基于Swap的局部搜索
 * @param graph
 * @param bsol
 * @param cstime
 * @return
 */
Solution swap_local_search_decomposition(const Graph &graph, const Solution &bsol, clock_t cstime)
{
	Solution csol = Solution(bsol);
	init_Gamma(graph, csol);

	// 获取每个分区中的元素序列
	unordered_set<int> clst_elems[graph.k]; // 创建 k 个 unordered_set<int>，保存每个分区中的点
	for (int i = 0; i < graph.nnode; i++) // 根据划分结果，将节点编号添加到对应的 unordered_set<int> 中
	{
		clst_elems[csol.ptn[i]].insert(i);
	}

#if(DEBUG)
	printf("debug--1.2\n");fflush(stdout);
#endif
	bool improved = true, flag;
	int iter = 0;
	while (improved)
	{
		improved = false;
#if(DEBUG)
		printf("debug--1.3\n");fflush(stdout);
#endif
		for (int pid1 = 0; pid1 < graph.k; pid1++)
		{
			for (int pid2 = pid1 + 1; pid2 < graph.k; pid2++)
			{
				if (active_matrix2[pid1][pid2] == 0) continue;
				active_matrix2[pid1][pid2] = active_matrix2[pid2][pid1] = 0;
				flag = false;

				// 将待遍历的点存储起来
				vector<int> vertices1Array, vertices2Array;
				for (int vtx1: clst_elems[pid1])
				{
					vertices1Array.push_back(vtx1);
				}
				for (int vtx2: clst_elems[pid2])
				{
					vertices2Array.push_back(vtx2);
				}

//				for (const int &vtx1: clst_elems[pid1])
//				{
//					for (const int &vtx2: clst_elems[pid2])
//					{
				for (const int &vtx1: vertices1Array)
				{
					for (const int &vtx2: vertices2Array)
					{
						if (csol.ptn[vtx1] == csol.ptn[vtx2]) continue;

						int weight = get_weight(graph, vtx1, vtx2); // 获取两个点之间的边权
#if(is_trueGamma)
						int delta = 2 * weight
								  + true_gamma[pid2][vtx1] + true_gamma[pid1][vtx2];
#else
						int delta = 2 * weight
						          + neg_gamma[pid2][vtx1] - neg_gamma[pid1][vtx1]
						          + pos_gamma[pid1][vtx1] - pos_gamma[pid2][vtx1]
						          + neg_gamma[pid1][vtx2] - neg_gamma[pid2][vtx2]
						          + pos_gamma[pid2][vtx2] - pos_gamma[pid1][vtx2];
#endif
						if (delta < 0)
						{
							Gain_node gn = Gain_node(vtx1, vtx2, -1, delta, 1);
							int bvtx1 = gn.elem1;
							int bvtx2 = gn.elem2;
							int bclst1 = csol.ptn[bvtx1];
							int bclst2 = csol.ptn[bvtx2];

							update_Gamma(graph, csol, gn);
							// 更新分区中的点队列
							clst_elems[bclst1].erase(bvtx1);
							clst_elems[bclst2].insert(bvtx1);
							clst_elems[bclst2].erase(bvtx2);
							clst_elems[bclst1].insert(bvtx2);
							improved = true;
							flag = true;
#if(is_Verify)
							assert(csol.cost == csol.cal_cost(graph));
#endif
						}
					}
				}

				// 更新 Gamma 表
				if (flag == true)
				{
					update_active_matrix(graph, pid1, pid2);
				}
			}

			iter++;
#if(DPTN)
			printf("SLS round %d, time=%.4f, ccost=%d, rbcost=%d\n", iter, (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC), csol.cost, rbcost);fflush(stdout);
#endif
		}
	} // end of while

	return csol;
}


/**
 * 基于邻域分解的局域搜索，一共测试三种模式：1.两种邻域都从初始解开始跑；2.先N1-add再N2-swap；3.先N2-swap再N1-add。
 * 写文章的时候，把N1和N2打包成一个整体，所以状态矩阵 active_matrix 只需要在最前面初始化一次
 * @param graph
 * @param csol
 * @param cstime
 * @param tlimit
 * @return
 */
Solution ND_based_LS(const Graph &graph, Solution &csol, clock_t cstime, double tlimit)
{
	Solution bsol = Solution(csol); // best_sol

#if(is_Verify)
	assert(csol.cost > 0);
#endif

#if(NDBLS_mode1) // 1.两种邻域都从初始解开始跑
	init_active_matrix(graph);
	Solution bsol0 = Solution(local_search_decomposition(graph, csol, cstime));  // Add
	Solution bsol1 = Solution(swap_local_search_decomposition(graph, csol, cstime));  // Swap

	if (bsol0.cost < bsol.cost)
	{
		bsol.cpy(bsol0);
		bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
	}
	if (bsol1.cost < bsol.cost)
	{
		bsol.cpy(bsol1);
		bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
	}
#elif(NDBLS_mode2) // 2.先N1-add再N2-swap
	init_active_matrix(graph);
	Solution bsol0 = Solution(local_search_decomposition(graph, csol, cstime));
	Solution bsol1 = Solution(swap_local_search_decomposition(graph, bsol0, cstime));

	if (bsol0.cost < bsol.cost)
	{
		bsol.cpy(bsol0);
		bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
	}
	if (bsol1.cost < bsol.cost)
	{
		bsol.cpy(bsol1);
		bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
	}
#elif(NDBLS_mode3) // 3.先N2-swap再N1-add
	init_active_matrix(graph);
	Solution bsol0 = Solution(swap_local_search_decomposition(graph, csol, cstime));
	Solution bsol1 = Solution(local_search_decomposition(graph, bsol0, cstime));

	if (bsol0.cost < bsol.cost)
	{
		bsol.cpy(bsol0);
		bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
	}
	if (bsol1.cost < bsol.cost)
	{
		bsol.cpy(bsol1);
		bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
	}
#elif(NDBLS_mode4) // 4.只有N1-add or 以概率选择N2-swap
	init_active_matrix(graph);
	Solution bsol0 = Solution(local_search_decomposition(graph, csol, cstime));
	if (bsol0.cost < bsol.cost)
	{
		bsol.cpy(bsol0);
		bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
	}
	#if(is_probSwap)
	double rnd = ((double) rand() / RAND_MAX);
	if (rnd < param_q /* K / graph.nnode*/)   // 依概率选择Swap，跟NDL/BL选择的概率设置一样
	{
		Solution bsol1 = Solution(swap_local_search_decomposition(graph, csol, cstime));
		if (bsol1.cost < bsol.cost)
		{
			bsol.cpy(bsol1);
			bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
		}
	}
	#endif
#endif

	return bsol;
}


/**
 * 有偏的局部搜索：【都是从同一个解出发】，先是三种Add LS，然后是以概率选择Swap LS
 * @param graph
 * @param csol
 * @param cstime
 * @param tlimit
 * @return
 */
Solution Biased_LS(const Graph &graph, Solution &csol, clock_t cstime, double tlimit)
{
	Solution bsol = Solution(csol); // best_sol

	Solution bsol0 = Solution(local_search(graph, csol, cstime));
	Solution bsol1 = Solution(biased_local_search1(graph, csol, cstime));
	Solution bsol2 = Solution(biased_local_search2(graph, csol, cstime));
	Solution bsol3 = Solution(biased_local_search3(graph, csol, cstime));
//	Solution bsol4 = Solution(biased_local_search4(graph, csol, cstime));

	if (bsol0.cost < bsol.cost)
	{
		bsol.cpy(bsol0);
		bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
	}
	if (bsol1.cost < bsol.cost)
	{
		bsol.cpy(bsol1);
		bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
	}
	if (bsol2.cost < bsol.cost)
	{
		bsol.cpy(bsol2);
		bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
	}
	if (bsol3.cost < bsol.cost)
	{
		bsol.cpy(bsol3);
		bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
	}
//	if (bsol4.cost < bsol.cost)
//	{
//		bsol.cpy(bsol4);
//		bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
//	}

#if(is_probSwap) // 这里Swap的可能性很小，才1%
	double rnd = ((double) rand() / RAND_MAX);
	if (rnd < param_q /* K / graph.nnode*/)   // 依概率选择Swap，跟NDL/BL选择的概率设置一样
	{
		cal_swap_value(graph, csol); // 更新权重、排序节点【因为一直要用Swap，每次使用前都要动态更新权重】
		Solution bsol5 = Solution(swap_local_search(graph, csol, cstime));
		if (bsol5.cost < bsol.cost) {
			bsol.cpy(bsol5);
			bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
		}
	}
#endif

	csol.cpy(bsol);  // rel_Maxima_search2()和rel_Maxima_search3()会用到
	return bsol;
}


/*
 * 纯随机扰动：随机选点，随机选择Add或Swap进行移动
 */
Solution shake(const Graph &graph, Solution &csol, const int np)
{
	for (int i = 0; i < np; i++)
	{
		// 选点
		int vtx1 = rand() % graph.nnode;
		int vtx2 = rand() % (graph.nnode - 1);
		if (vtx2 >= vtx1)
			vtx2 += 1;
		int clst1 = csol.ptn[vtx1];
		int clst2 = csol.ptn[vtx2];

		// 扰动
		if (clst1 != clst2)
		{
			int type = rand() % 2;
			if (type == 0 && csol.sc[clst1] > 1) // Add
			{
				csol.sc[clst1]--;
				csol.ptn[vtx1] = clst2;
				csol.sc[clst2]++;
			}
			else if (type == 1) // Swap
			{
				csol.ptn[vtx1] = clst2;
				csol.ptn[vtx2] = clst1;
			}
		}
	}

	csol.cost = csol.cal_cost(graph);
	if (csol.cost < rbcost)
	{
		rbcost = csol.cost;
		rbtime = (clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
	}
	return csol;
}


/*
 * 扰动（有偏扰动），效果差，放弃
 */
Solution shake1(const Graph &graph, Solution &csol, const int np)
{
	int nptr = 0;
	calculate_v4(graph, csol);
	for (int i = 0; i < graph.nnode; i++)
	{
		// 选点
		int vtx1 = desc_nodes4[i];  // 这里是降序排序，因为是把那些对解贡献差的点扰动掉
		int clst1 = csol.ptn[vtx1];

		if (csol.sc[clst1] == 1)
			continue;

		int clst2 = rand() % (graph.k - 1);
		while (clst2 == clst1)
			clst2 = rand() % (graph.k - 1);

		// 扰动
		csol.sc[clst1]--;
		csol.ptn[vtx1] = clst2;
		csol.sc[clst2]++;

		nptr++;
		if (nptr >= np)
			break;
	}

	csol.cost = csol.cal_cost(graph);
	if (csol.cost < rbcost)
	{
		rbcost = csol.cost;
		rbtime = (clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
	}

	return csol;
}


//====================================================================
//=========================== 四个版本 ================================
//====================================================================
///*
// * MS算法，最初始的版本
// */
//Solution Maxima_search(const Graph &graph, Solution &csol, clock_t cstime, double tlimit)
//{
//	Solution bsol = Solution(csol); // best_sol
//
//#if(is_Verify)
//	assert(csol.cost > 0);
//#endif
////	Solution bsol0 = Solution(local_search_decomposition(graph, csol, cstime));
//	Solution bsol1 = Solution(local_search(graph, csol, cstime));
//	Solution bsol2 = Solution(biased_local_search1(graph, csol, cstime));
//	Solution bsol3 = Solution(biased_local_search2(graph, csol, cstime));
//	Solution bsol4 = Solution(biased_local_search3(graph, csol, cstime));
//
////	if (bsol0.cost < bsol.cost)
////	{
////		bsol.cpy(bsol0);
////		bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
////	}
//	if (bsol1.cost < bsol.cost)
//	{
//		bsol.cpy(bsol1);
//		bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
//	}
//	if (bsol2.cost < bsol.cost)
//	{
//		bsol.cpy(bsol2);
//		bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
//	}
//	if (bsol3.cost < bsol.cost)
//	{
//		bsol.cpy(bsol3);
//		bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
//	}
//	if (bsol4.cost < bsol.cost)
//	{
//		bsol.cpy(bsol4);
//		bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
//	}
//	csol.cpy(bsol);
//
//	// 开始重启迭代
//	int non_improve = 0; // restart_iter
//	while (non_improve < max_nipv && (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC) < tlimit)
//	{
//#if(DEBUG)
//		printf("debug--1.1\n");fflush(stdout);
//#endif
//		shake(graph, csol, wp);
//
//		// 生成初始解
//		Solution nsol0 = Solution(local_search_decomposition(graph, csol, cstime));
//		Solution nsol1 = Solution(local_search(graph, csol, cstime));
//		Solution nsol2 = Solution(biased_local_search1(graph, csol, cstime));
//		Solution nsol3 = Solution(biased_local_search2(graph, csol, cstime));
//		Solution nsol4 = Solution(biased_local_search3(graph, csol, cstime));
//
//		if (nsol0.cost < csol.cost) csol.cpy(nsol0);
//		if (nsol1.cost < csol.cost) csol.cpy(nsol1);
//		if (nsol2.cost < csol.cost) csol.cpy(nsol2);
//		if (nsol3.cost < csol.cost) csol.cpy(nsol3);
//		if (nsol4.cost < csol.cost) csol.cpy(nsol4);
//
//		if (csol.cost < bsol.cost)
//		{
//			if (csol.cost < rbcost)
//			{
//				rbcost = csol.cost;
//				rbtime = (clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
//			}
//			bsol.cpy(csol);
//			bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
//			non_improve = 0;
//		}
//		else
//			non_improve++;
//
//#if(is_maxSearch)
//		csol.cpy(bsol); // maxima search，有点类似于intensification
//#endif
//
//		printf("IMS weak ni:%d, time=%.4f, best cost=%d\n",
//				non_improve,
//				(clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC),
//				bsol.cost);fflush(stdout);
//	}
//
//	return bsol;
//}


/*
 * MS算法
 */
Solution rel_Maxima_search(const Graph &graph, Solution &csol, clock_t cstime, double tlimit)
{
	Solution bsol = Solution(csol); // best_sol

#if(is_Verify)
	assert(csol.cost > 0);
#endif
	Solution bsol0 = Solution(local_search(graph, csol, cstime));
	Solution bsol1 = Solution(biased_local_search1(graph, csol, cstime));
	Solution bsol2 = Solution(biased_local_search2(graph, csol, cstime));
	Solution bsol3 = Solution(biased_local_search3(graph, csol, cstime));
//	Solution bsol4 = Solution(biased_local_search4(graph, csol, cstime));

	if (bsol0.cost < bsol.cost)
	{
		bsol.cpy(bsol0);
		bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
	}
	if (bsol1.cost < bsol.cost)
	{
		bsol.cpy(bsol1);
		bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
	}
	if (bsol2.cost < bsol.cost)
	{
		bsol.cpy(bsol2);
		bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
	}
	if (bsol3.cost < bsol.cost)
	{
		bsol.cpy(bsol3);
		bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
	}
//	if (bsol4.cost < bsol.cost)
//	{
//		bsol.cpy(bsol4);
//		bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
//	}

#if(is_probSwap)
	double rnd = ((double) rand() / RAND_MAX);
	if (rnd < param_q /** K / graph.nnode*/)
	{
		cal_swap_value(graph, csol);
		Solution bsol5 = Solution(swap_local_search(graph, csol, cstime));
		if (bsol5.cost < bsol.cost) {
			bsol.cpy(bsol5);
			bsol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
		}
//		counter++;
//		cout << "***SWAP works!***" << endl;
	}
#endif

	// 开始重启迭代
	int non_improve = 0; // restart_iter
	while (non_improve < max_nipv && (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC) < tlimit)
	{
#if(DEBUG)
		printf("debug--1.1\n");fflush(stdout);
#endif
		csol.cpy(bsol); // maxima search，有点类似于intensification
		csol = shake(graph, csol, wp);

		// 生成初始解
		Solution nsol0 = Solution(local_search(graph, csol, cstime));
		Solution nsol1 = Solution(biased_local_search1(graph, csol, cstime));
		Solution nsol2 = Solution(biased_local_search2(graph, csol, cstime));
		Solution nsol3 = Solution(biased_local_search3(graph, csol, cstime));
//		Solution nsol4 = Solution(biased_local_search4(graph, csol, cstime));

		if (nsol0.cost < csol.cost)
		{
			csol.cpy(nsol0);
			csol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
		}
		if (nsol1.cost < csol.cost)
		{
			csol.cpy(nsol1);
			csol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
		}
		if (nsol2.cost < csol.cost)
		{
			csol.cpy(nsol2);
			csol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
		}
		if (nsol3.cost < csol.cost)
		{
			csol.cpy(nsol3);
			csol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
		}
//		if (nsol4.cost < csol.cost)
//		{
//			csol.cpy(nsol4);
//			csol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
//		}

#if(is_probSwap)
		double rnd = ((double) rand() / RAND_MAX);
		if (rnd < param_q /** K / graph.nnode*/)
		{
			cal_swap_value(graph, csol);
			Solution nsol5 = Solution(swap_local_search(graph, csol, cstime));
			if (nsol5.cost < csol.cost) {
				csol.cpy(nsol5);
				csol.btime = (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC);
			}
		}
#endif

		if (csol.cost < bsol.cost)
		{
			if (csol.cost < rbcost)
			{
				rbcost = csol.cost;
				rbtime = csol.btime;
			}
			bsol.cpy(csol);

			printf("IMS weak ni:%d, time=%.4f, best cost=%d\n",
					non_improve,
					(clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC),
					bsol.cost);fflush(stdout);

			non_improve = 0;
		}
		else
			non_improve++;
	}

	return bsol;
}


/*
 * MS算法（以概率选择NBL和BL）
 */
Solution rel_Maxima_search1(const Graph &graph, Solution &csol, clock_t cstime, double tlimit)
{
	Solution bsol = Solution(csol); // best_sol
	Solution nsol = Solution(csol);

#if(is_Verify)
	assert(csol.cost > 0);
#endif

	double rnd = ((double) rand() / RAND_MAX);
	if (rnd < param_b/* * K / graph.nnode*/)  // param_b = Q
	{
		nsol = ND_based_LS(graph, csol, cstime, tlimit);
	}
	else
	{
		nsol = Biased_LS(graph, csol, cstime, tlimit);
	}
	if (nsol.cost < csol.cost) csol.cpy(nsol);
	if (csol.cost < bsol.cost)
	{
		if (csol.cost < rbcost)
		{
			rbcost = csol.cost;
			rbtime = csol.btime;
		}
		bsol.cpy(csol);
	}

	// 开始重启迭代
	int non_improve = 0; // restart_iter
	while (non_improve < max_nipv && (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC) < tlimit)
	{
#if(DEBUG)
		printf("debug--1.1\n");fflush(stdout);
#endif
		csol.cpy(bsol); // maxima search，有点类似于intensification
		shake(graph, csol, wp);

		double rnd = ((double) rand() / RAND_MAX);
		if (rnd < param_b /** K / graph.nnode*/)
		{
			nsol = ND_based_LS(graph, csol, cstime, tlimit);
		}
		else
		{
			nsol = Biased_LS(graph, csol, cstime, tlimit);
		}
		if (nsol.cost < csol.cost) csol.cpy(nsol);
		if (csol.cost < bsol.cost)
		{
			if (csol.cost < rbcost)
			{
				rbcost = csol.cost;
				rbtime = csol.btime;
			}
			bsol.cpy(csol);

			printf("IMS weak ni:%d, time=%.4f, best cost=%d\n",
					non_improve,
					(clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC),
					bsol.cost);fflush(stdout);

			non_improve = 0;
		}
		else
			non_improve++;
	}

	return bsol;
}


/*
 * MS算法（先NBL，再BL），顺序执行，不是从同一个解出发
 */
Solution rel_Maxima_search2(const Graph &graph, Solution &csol, clock_t cstime, double tlimit)
{
	Solution bsol = Solution(csol); // best_sol

#if(is_Verify)
	assert(csol.cost > 0);
#endif

	// 注意，这里是顺序执行，不是从同一个解出发
	Solution nsol1 = ND_based_LS(graph, csol, cstime, tlimit);
	Solution nsol2 = Biased_LS(graph, nsol1, cstime, tlimit);
	if (nsol2.cost < csol.cost) csol.cpy(nsol2);
	if (csol.cost < bsol.cost)
	{
		if (csol.cost < rbcost)
		{
			rbcost = csol.cost;
			rbtime = csol.btime;
		}
		bsol.cpy(csol);
	}

	// 开始重启迭代
	int non_improve = 0; // restart_iter
	while (non_improve < max_nipv && (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC) < tlimit)
	{
#if(DEBUG)
		printf("debug--1.1\n");fflush(stdout);
#endif
		csol.cpy(bsol); // maxima search，有点类似于intensification
		csol = shake(graph, csol, wp);

		Solution nsol1 = ND_based_LS(graph, csol, cstime, tlimit);
		Solution nsol2 = Biased_LS(graph, nsol1, cstime, tlimit);
		if (nsol2.cost < csol.cost) csol.cpy(nsol2);
		if (csol.cost < bsol.cost)
		{
			if (csol.cost < rbcost)
			{
				rbcost = csol.cost;
				rbtime = csol.btime;
			}
			bsol.cpy(csol);

			printf("IMS weak ni:%d, time=%.4f, best cost=%d\n",
					non_improve,
					(clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC),
					bsol.cost);fflush(stdout);

			non_improve = 0;
		}
		else
			non_improve++;
	}

	return bsol;
}


/*
 * MS算法（先BL，再NBL），顺序执行，不是从同一个解出发
 */
Solution rel_Maxima_search3(const Graph &graph, Solution &csol, clock_t cstime, double tlimit)
{
	Solution bsol = Solution(csol); // best_sol

#if(is_Verify)
	assert(csol.cost > 0);
#endif

	// 注意，这里是顺序执行，不是从同一个解出发
	Solution nsol1 = Biased_LS(graph, csol, cstime, tlimit);
	Solution nsol2 = ND_based_LS(graph, nsol1, cstime, tlimit);
	if (nsol2.cost < csol.cost) csol.cpy(nsol2);
	if (csol.cost < bsol.cost)
	{
		if (csol.cost < rbcost)
		{
			rbcost = csol.cost;
			rbtime = csol.btime;
		}
		bsol.cpy(csol);
	}

	// 开始重启迭代
	int non_improve = 0; // restart_iter
	while (non_improve < max_nipv && (clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC) < tlimit)
	{
#if(DEBUG)
		printf("debug--1.1\n");fflush(stdout);
#endif
		csol.cpy(bsol); // maxima search，有点类似于intensification
		csol = shake(graph, csol, wp);

		Solution nsol1 = Biased_LS(graph, csol, cstime, tlimit);
		Solution nsol2 = ND_based_LS(graph, nsol1, cstime, tlimit);
		if (nsol2.cost < csol.cost) csol.cpy(nsol2);
		if (csol.cost < bsol.cost)
		{
			if (csol.cost < rbcost)
			{
				rbcost = csol.cost;
				rbtime = csol.btime;
			}
			bsol.cpy(csol);

			printf("IMS weak ni:%d, time=%.4f, best cost=%d\n",
					non_improve,
					(clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC),
					bsol.cost);fflush(stdout);

			non_improve = 0;
		}
		else
			non_improve++;
	}

	return bsol;
}


/*
 * IMS
 */
void iterated_maxima_search(const Graph &graph, Solution &csol, double tlimit)
{
#if(DEBUG)
	printf("debug--2.1\n");fflush(stdout);
#endif
	Solution bsol = Solution(csol); // best_sol
	sp = int(pct_sp * graph.nnode); // strong perturbation strength
	wp = int(pct_wp * graph.nnode); // weak perturbation strength

#if(DEBUG)
	printf("debug--2.2\n");fflush(stdout);
#endif
	// 开始重启迭代
	int iter = 0;
	clock_t cstime = clock();
#if(DEBUG)
	printf("debug--2.3\n");fflush(stdout);
#endif
	while ((clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC) < tlimit)
	{
#if(rMS_mode1)
		Solution nsol1 = Solution(rel_Maxima_search1(graph, csol, cstime, tlimit));
#elif(rMS_mode2)
		Solution nsol1 = Solution(rel_Maxima_search2(graph, csol, cstime, tlimit));
#elif(rMS_mode3)
		Solution nsol1 = Solution(rel_Maxima_search3(graph, csol, cstime, tlimit));
#elif(rMS_mode4)
		Solution nsol1 = Solution(rel_Maxima_search(graph, csol, cstime, tlimit));
#endif

#if(DEBUG)
		printf("debug--2.4\n");fflush(stdout);
#endif
		if (nsol1.cost < bsol.cost)
		{
			bsol.cpy(nsol1);
		}
		printf("IMS strong restarts:%d, time=%.4f, best cost=%d\n\n",
				iter,
				(clock() - cstime) / static_cast<double>(CLOCKS_PER_SEC),
				bsol.cost);fflush(stdout);

		iter++;

		nsol1.cpy(bsol); // TODO maxima search，有点类似于intensification的搜索
		csol.cpy(shake(graph, nsol1, sp));
	}

	printf("IMS FINAL restarts:%d, btime=%f, best cost=%d\n",
			iter,
			bsol.btime,
			bsol.cost);fflush(stdout);

	// 写入文件
#if(DEBUG)
	printf("debug--2.9\n");fflush(stdout);
#endif
	char outputfile[1000];
	char *graph_name = basename(filename);
	sprintf(outputfile, "%soutput_dir/results_%s_%d.txt", root_path.c_str(), graph_name, K);  // root_path.c_str(): "./"
	FILE *opf = fopen(outputfile, "a");
	if (!opf)
	{
		perror("Failed to open output file");
		exit(-1);
	}
	fprintf(opf, "IMS restarts:%d, IMS best cost=%d, IMS best time=%.4f\n", iter, rbcost, rbtime);
	fclose(opf);

	// 返回IMS的最好解
	csol.cpy(bsol);
}


/*
 * 计算常用指标并输出：best cost, average cost, average time, hit, dev
 */
void cal_indicators(int *each_run_rlt, int *each_run_time, const int &runs, int &bcost, double &avg_cost,
                    double &avg_time)
{
	double sum_cost = 0, sum_time = 0;
	std_dev = 0;

	bcost = MAX_VALUE;
	for (int i = 0; i < runs; i++)
	{
		sum_cost += each_run_rlt[i];
		sum_time += each_run_time[i];

		if (each_run_rlt[i] < bcost)
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
 * 验证分区和cost
 */
void verify(const Graph &graph, Solution &bsol)
{
	// 1.验证分区数
	for (int i = 0; i < graph.k; i++)
	{
		if (bsol.sc[i] <= 0)
		{
			cerr << "分区" << i << "的点数小于0" << endl;
			exit(-999);
		}

	}

	// 2.验证cost
	int vcost = bsol.cal_cost(graph);
	if (vcost != bsol.cost)
		cerr << "cost验证未通过" << endl;

	// 3.
}


/*
 * 将解写入文件
 */
void write_IMSSDN_sol(const Graph &graph, const Solution &csol, char *instancefile, const int &fno)
{
	// 获取文件名
	char *graph_name = basename(instancefile);

	// 拼接文件路径
	char rltfile_detail[1000];
	sprintf(rltfile_detail, "%sresults/IMSS_%s_%d_%d.txt", root_path.c_str(), graph_name, graph.k, fno); // fno = seed

	// 用 ofstream 打开文件
	ofstream fout(rltfile_detail);
	if (!fout.is_open())
	{
		cerr << "无法打开文件: " << rltfile_detail << endl;
		exit(EXIT_FAILURE);
	}

	// 写入 cost
	fout << "cost=" << csol.cost << endl;
    fout << "cost=" << csol.cost << endl;
    fout << "best_cost=" << rbcost << endl;
    fout << "best_time=" << rbtime << endl;

	// 写入 ptn 数组
	//
	fout << "ptn=";
	for (int i = 0; i < graph.nnode; i++)
	{
		fout << csol.ptn[i] << " ";
	}
	fout << endl;

	// 写入 sc 数组
	// size of each cluster
	fout << "sc[]=";
	for (int i = 0; i < graph.k; i++)
	{
		fout << csol.sc[i] << " ";
	}
	fout << endl;

	fout.close();
}


/*
 * 执行算例instancefile
 */
void IMS_run(char *instancefile, double timelimit)
{
	// 读图
	Graph graph = Graph(instancefile, K);
	allocate_memory(graph);

	// 这里是从节点的角度，计算每个节点的影响力
	calculate_v1(graph); // 按 节点的影响力=∑(|正边权|+|负边权|) 确定节点的遍历顺序
	calculate_v2(graph); // 按 节点的影响力=∑|正边权| 确定节点的遍历顺序
	calculate_v3(graph); // 按 节点的影响力=∑|负边权| 确定节点的遍历顺序

	// run 10 times algorithms
	int crun = 0;
	while (crun < runs)
	{
		output_header(crun + 1);
		start = clock();
		rbcost = MAX_VALUE;
		rbtime = MAX_VALUE;

		// 1.RH 作为初始解，直接读进来（10个解，每个run读一个）
#if(DEBUG)
		printf("debug--1\n");fflush(stdout);
#endif
		Solution bsol = Solution();
		read_RH_sol(graph, bsol, instancefile, crun); // 当文件没有正确读进来时，这里会构造一个初始不可行解（边界错误）
		verify(graph, bsol); // XXX
//		bsol.verify(graph); // 只验证cost

		// 2.IMSS
#if(DEBUG)
		printf("debug--2\n");fflush(stdout);
#endif
		iterated_maxima_search(graph, bsol, timelimit);
		verify(graph, bsol); // XXX
		write_IMSSDN_sol(graph, bsol, instancefile, seed);

		printf("IMSS Round %d: best cost=%d, best time=%.4f\n", crun, rbcost, rbtime);

		crun++;
	}

	// 释放内存
	free_memory(graph);
}


#endif /* IMS_V3_0_H_ */
