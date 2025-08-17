/*
 * RH.h
 *
 *  Created on: 2025年4月20日
 *      Author: admin
 */

#ifndef RH_H_
#define RH_H_

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


//Graph graph;
clock_t start, end;

//char filename[1000] = "./Benchmarks/R0_signed_graph/slashdot-zoo.graph";
char filename[1000] = "./instances/soc-sign-epinions.graph";

int seed = 14;
double time_limit = 200;
int K = 2;
int runs = 1;

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
	}
}


/*
 * 分配内存
 */
void allocate_memory(const Graph &graph)
{
}


/*
 * 释放内存
 */
void free_memory(const Graph &graph)
{
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
		printf("init_cur_cost=%d\n", csol.cost);fflush(stdout);
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
						csol.sc[clst1]--;
						csol.ptn[vtx1] = clst2;
						csol.sc[clst2]++;
						improved = true;
						csol.cost = csol.cal_cost(graph);
						printf("cur_cost=%d\n", csol.cost);fflush(stdout);
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
	sprintf(outputfile, "./output_dir/results_%s_%d.txt", graph_name, K);
	FILE *opf = fopen(outputfile, "a");
	if (!opf) {
		perror("Failed to open RH output file");
		exit(-1);
	}
	fprintf(opf, "RH restart rounds:%d, RH best cost=%d, RH best time=%.4f\n", riter, bsol.cost, bsol.btime);
	fclose(opf);
}


/*
 * 将初始解写入文件
 */
void write_RH_sol(const Graph &graph, const Solution &csol, char *instancefile, const int &fno) {
    // 获取文件名
    char *graph_name = basename(instancefile);

    // 拼接文件路径
    char rltfile[1000];
    sprintf(rltfile, "./init_sols/%s_%d_%d", graph_name, graph.k, fno);

    // 用 ofstream 打开文件
    ofstream fout(rltfile);
    if (!fout.is_open()) {
        cerr << "无法打开文件: " << rltfile << endl;
        exit(EXIT_FAILURE);
    }

    // 写入 cost
    fout << csol.cost << endl;

    // 写入 ptn 数组
    for (int i = 0; i < graph.nnode; i++) {
        fout << csol.ptn[i] << " ";
    }
    fout << endl;

    // 写入 sc 数组
    for (int i = 0; i < graph.k; i++) {
        fout << csol.sc[i] << " ";
    }
    fout << endl;

    fout.close();
}


void RH_run(char *instancefile, double timelimit)
{
	// 读图
	Graph graph = Graph(instancefile, K);
	allocate_memory(graph);

	// run 10 times algorithms
	int crun = 0;
	while(crun < runs)
	{
		start = clock();

		// 1.RH
#if(DEBUG)
		printf("debug--1\n");fflush(stdout);
#endif
		Solution bsol = Solution();
		multistart_relocation_heuristic(graph, bsol, timelimit);
		write_RH_sol(graph, bsol, instancefile, seed);
		crun++;
//		exit(-1);
	}

	// 释放内存
	free_memory(graph);
}



#endif /* RH_H_ */
