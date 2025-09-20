//============================================================================
// Name        : 0.cpp
// Author      : Lin
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================
#ifdef _WIN32
#include <windows.h>
#endif
#include <iostream>
using namespace std;

#define isIMS 1
#define isTS 0
#define isVNS 0
#define isILS 0

#define DEBUG 0
#define DPTN 0

bool tuning = false;

#if(isIMS)
//#include "IMS_v1.0.h"
//#include "IMS_v2.0_swap.h"
// #include "statistic.h"
// #include "IMS_v2.0_swap_dc_purturb.h"
#include "IMS_v2.3_decomposition3.h"
#elif(isTS)
#include "TS.h"
#elif(isVNS)
#include "VNS.h"
#else
#include "RH.h"
#endif

int main(int argc, char **argv)
{
	cout << "!!! begin running!!!" << endl;
#ifdef _WIN32
	// 设置控制台代码页为 UTF-8
	SetConsoleOutputCP(CP_UTF8);
	// 或者使用中文代码页
	// SetConsoleOutputCP(936);
#endif

	// 输入参数
	if (tuning)
	{
		root_path = "./";
		Read_Parameters(argc, argv);
		srand(seed);
		string rltfile = "./results/result";
#if(isIMS)
		IMS_run(filename, time_limit);
#elif(isTS)
		TS_run(filename, rltfile, time_limit);
#elif(isVNS)
		VNS_run(filename, rltfile, time_limit);
#else
		RH_run(filename, time_limit);
#endif
	}
	else
	{
		root_path = "../";
		srand(seed);
		string dir = root_path + "instances/";
		string rltfile = root_path + "results/result";
		string listfile = "instances_list.txt";
		string path = dir + listfile;
		string *instance_list;
		double *tlist;
		int *klist;
		ifstream fin(path); // 打开文件

		// 判断是否成功打开文件
		if (!fin.is_open())
		{
			cerr << "Can not open the file!" << path << endl;
			exit(-14);
		}
		if (fin.fail()) // 检查文件流是否处于错误状态
		{
			cerr << "Error occurred during file operation." << path << endl;
			exit(-15);
		}
		if (fin.eof()) // 判断文件是否为空
		{
			cerr << "Empty file" << path << endl;
			exit(-16);
		}

		// 开始读内容
		int num_instances;
		fin >> num_instances;
		instance_list = new string[num_instances]; // 分配内存，并调整大小
		tlist = new double[num_instances];
		klist = new int[num_instances];

		int temp = 0;
		while (!fin.eof())
		{
			fin >> instance_list[temp];
			fin >> tlist[temp];
			fin >> klist[temp++];
			//			fin >> target_cost_list[temp++];
		}

		fin.close();

		// 遍历测试所有算例
		for (int i = 0; i < num_instances; i++)
		{
			strncpy(filename, instance_list[i].c_str(), 1000);

			time_limit = tlist[i];
			K = klist[i];
#if(isIMS)
			IMS_run(filename, time_limit);
#elif(isTS)
		TS_run(filename, rltfile, time_limit);
#elif(isVNS)
		VNS_run(filename, rltfile, time_limit);
#else
		RH_run(filename, time_limit);
#endif
		}

		// 释放内存
		delete[] klist;
		delete[] tlist;
		delete[] instance_list;
	}

	cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!

	return 0;
}
