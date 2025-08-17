/*
 * tool.h
 *
 *  Created on: 2025年3月9日
 *      Author: LinLin
 */

#ifndef TOOLS_H_
#define TOOLS_H_


void Generate_Rand_List(int *randlist, int len)
{
	assert(randlist != NULL);

	for (int i = 0; i < len; i++)
		randlist[i] = i;
//	printf("打乱前");
//	for(int i = 0; i < len; i++)
//	{
//		printf("%d ", randlist[i]);
//		fflush(stdout);
//	}
//	printf("\n");
//	fflush(stdout);

	// swap
	for (int i = 0; i < len; i++)
	{
		int randid = rand() % len;
		int tmp = randlist[i];
		randlist[i] = randlist[randid];
		randlist[randid] = tmp;
	}

//	printf("打乱后");
//	for(int i = 0; i < len; i++)
//	{
//		printf("%d ", randlist[i]);
//		fflush(stdout);
//	}
//	printf("\n");
//	fflush(stdout);
}


/* 升序排列 */
void Quick_Sort_asc(int *idlist, int *objlist, int l, int r)
{
	if (l < r)
	{
		int left = l, right = r;
		int pivot = objlist[l];							// 原本是double类型数据
		int id_pivot = idlist[l];
		while (left < right)
		{
			// Find first > pivot, from right to left
			while (left < right && objlist[right] >= pivot)
				right--;
			if (left < right)
			{
				objlist[left] = objlist[right];
				idlist[left++] = idlist[right];
			}

			// Find first <= pivot, from left to right
			while (left < right && objlist[left] < pivot)
				left++;
			if (left < right)
			{
				objlist[right] = objlist[left];
				idlist[right--] = idlist[left];
			}
		}
		objlist[left] = pivot;
		idlist[left] = id_pivot;
		Quick_Sort_asc(idlist, objlist, l, left - 1);
		Quick_Sort_asc(idlist, objlist, left + 1, r);
	}
}


/* 降序排列 */
void Quick_Sort_desc(int *idlist, int *objlist, int l, int r)
{
	if (l < r)
	{
		int left = l, right = r;
		int pivot = objlist[l];							// 原本是double类型数据
		int id_pivot = idlist[l];
		while (left < right)
		{
			// Find first > pivot, from right to left
			while (left < right && objlist[right] <= pivot)
				right--;
			if (left < right)
			{
				objlist[left] = objlist[right];
				idlist[left++] = idlist[right];
			}

			// Find first <= pivot, from left to right
			while (left < right && objlist[left] >= pivot)
				left++;
			if (left < right)
			{
				objlist[right] = objlist[left];
				idlist[right--] = idlist[left];
			}
		}
		objlist[left] = pivot;
		idlist[left] = id_pivot;
		Quick_Sort_desc(idlist, objlist, l, left - 1);
		Quick_Sort_desc(idlist, objlist, left + 1, r);
	}
}



#endif /* TOOLS_H_ */
