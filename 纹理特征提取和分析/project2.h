#include<stdio.h>
#include <iostream>
#include <opencv2\opencv.hpp>
using namespace cv;
using namespace std;

//灰度等级设为16
const int gray_level = 16;

//0度灰度共生矩阵
void getglcm_horison(Mat& input, Mat& dst)
{
	Mat src = input;
	CV_Assert(1 == src.channels());
	src.convertTo(src, CV_32S);
	int height = src.rows;
	int width = src.cols;
	int max_gray_level = 0;
	for (int j = 0; j < height; j++)//寻找像素灰度最大值
	{
		int* srcdata = src.ptr<int>(j);
		for (int i = 0; i < width; i++)
		{
			if (srcdata[i] > max_gray_level)
			{
				max_gray_level = srcdata[i];
			}
		}
	}
	max_gray_level++;//像素灰度最大值加1即为该矩阵所拥有的灰度级数
	if (max_gray_level > 16)//若灰度级数大于16，则将图像的灰度级缩小至16级，减小灰度共生矩阵的大小。
	{
		for (int i = 0; i < height; i++)
		{
			int* srcdata = src.ptr<int>(i);
			for (int j = 0; j < width; j++)
			{
				srcdata[j] = (int)srcdata[j] / gray_level;
			}
		}

		dst.create(gray_level, gray_level, CV_32SC1);
		dst = Scalar::all(0);
		for (int i = 0; i < height; i++)
		{
			int* srcdata = src.ptr<int>(i);
			for (int j = 0; j < width - 1; j++)
			{
				int rows = srcdata[j];
				int cols = srcdata[j + 1];
				dst.ptr<int>(rows)[cols]++;
			}
		}
	}
	else//若灰度级数小于16，则生成相应的灰度共生矩阵
	{
		dst.create(max_gray_level, max_gray_level, CV_32SC1);
		dst = Scalar::all(0);
		for (int i = 0; i < height; i++)
		{
			int* srcdata = src.ptr<int>(i);
			for (int j = 0; j < width - 1; j++)
			{
				int rows = srcdata[j];
				int cols = srcdata[j + 1];
				dst.ptr<int>(rows)[cols]++;
			}
		}
	}
}
//90度灰度共生矩阵
void getglcm_vertical(Mat& input, Mat& dst)
{
	Mat src = input;
	CV_Assert(1 == src.channels());
	src.convertTo(src, CV_32S);
	int height = src.rows;
	int width = src.cols;
	int max_gray_level = 0;
	for (int j = 0; j < height; j++)
	{
		int* srcdata = src.ptr<int>(j);
		for (int i = 0; i < width; i++)
		{
			if (srcdata[i] > max_gray_level)
			{
				max_gray_level = srcdata[i];
			}
		}
	}
	max_gray_level++;
	if (max_gray_level > 16)
	{
		for (int i = 0; i < height; i++)//将图像的灰度级缩小至16级，减小灰度共生矩阵的大小。
		{
			int* srcdata = src.ptr<int>(i);
			for (int j = 0; j < width; j++)
			{
				srcdata[j] = (int)srcdata[j] / gray_level;
			}
		}

		dst.create(gray_level, gray_level, CV_32SC1);
		dst = Scalar::all(0);
		for (int i = 0; i < height - 1; i++)
		{
			int* srcdata = src.ptr<int>(i);
			int* srcdata1 = src.ptr<int>(i + 1);
			for (int j = 0; j < width; j++)
			{
				int rows = srcdata[j];
				int cols = srcdata1[j];
				dst.ptr<int>(rows)[cols]++;
			}
		}
	}
	else
	{
		dst.create(max_gray_level, max_gray_level, CV_32SC1);
		dst = Scalar::all(0);
		for (int i = 0; i < height - 1; i++)
		{
			int* srcdata = src.ptr<int>(i);
			int* srcdata1 = src.ptr<int>(i + 1);
			for (int j = 0; j < width; j++)
			{
				int rows = srcdata[j];
				int cols = srcdata1[j];
				dst.ptr<int>(rows)[cols]++;
			}
		}
	}
}
//45度灰度共生矩阵
void getglcm_45(Mat& input, Mat& dst)
{
	Mat src = input;
	CV_Assert(1 == src.channels());
	src.convertTo(src, CV_32S);
	int height = src.rows;
	int width = src.cols;
	int max_gray_level = 0;
	for (int j = 0; j < height; j++)
	{
		int* srcdata = src.ptr<int>(j);
		for (int i = 0; i < width; i++)
		{
			if (srcdata[i] > max_gray_level)
			{
				max_gray_level = srcdata[i];
			}
		}
	}
	max_gray_level++;
	if (max_gray_level > 16)
	{
		for (int i = 0; i < height; i++)//将图像的灰度级缩小至16级，减小灰度共生矩阵的大小。
		{
			int* srcdata = src.ptr<int>(i);
			for (int j = 0; j < width; j++)
			{
				srcdata[j] = (int)srcdata[j] / gray_level;
			}
		}

		dst.create(gray_level, gray_level, CV_32SC1);
		dst = Scalar::all(0);
		for (int i = 0; i < height - 1; i++)
		{
			int* srcdata = src.ptr<int>(i);
			int* srcdata1 = src.ptr<int>(i + 1);
			for (int j = 0; j < width - 1; j++)
			{
				int rows = srcdata[j];
				int cols = srcdata1[j + 1];
				dst.ptr<int>(rows)[cols]++;
			}
		}
	}
	else
	{
		dst.create(max_gray_level, max_gray_level, CV_32SC1);
		dst = Scalar::all(0);
		for (int i = 0; i < height - 1; i++)
		{
			int* srcdata = src.ptr<int>(i);
			int* srcdata1 = src.ptr<int>(i + 1);
			for (int j = 0; j < width - 1; j++)
			{
				int rows = srcdata[j];
				int cols = srcdata1[j + 1];
				dst.ptr<int>(rows)[cols]++;
			}
		}
	}
}
//135度灰度共生矩阵
void getglcm_135(Mat& input, Mat& dst)
{
	Mat src = input;
	CV_Assert(1 == src.channels());
	src.convertTo(src, CV_32S);
	int height = src.rows;
	int width = src.cols;
	int max_gray_level = 0;
	for (int j = 0; j < height; j++)
	{
		int* srcdata = src.ptr<int>(j);
		for (int i = 0; i < width; i++)
		{
			if (srcdata[i] > max_gray_level)
			{
				max_gray_level = srcdata[i];
			}
		}
	}
	max_gray_level++;
	if (max_gray_level > 16)
	{
		for (int i = 0; i < height; i++)//将图像的灰度级缩小至16级，减小灰度共生矩阵的大小。
		{
			int* srcdata = src.ptr<int>(i);
			for (int j = 0; j < width; j++)
			{
				srcdata[j] = (int)srcdata[j] / gray_level;
			}
		}

		dst.create(gray_level, gray_level, CV_32SC1);
		dst = Scalar::all(0);
		for (int i = 0; i < height - 1; i++)
		{
			int* srcdata = src.ptr<int>(i);
			int* srcdata1 = src.ptr<int>(i + 1);
			for (int j = 1; j < width; j++)
			{
				int rows = srcdata[j];
				int cols = srcdata1[j - 1];
				dst.ptr<int>(rows)[cols]++;
			}
		}
	}
	else
	{
		dst.create(max_gray_level, max_gray_level, CV_32SC1);
		dst = Scalar::all(0);
		for (int i = 0; i < height - 1; i++)
		{
			int* srcdata = src.ptr<int>(i);
			int* srcdata1 = src.ptr<int>(i + 1);
			for (int j = 1; j < width; j++)
			{
				int rows = srcdata[j];
				int cols = srcdata1[j - 1];
				dst.ptr<int>(rows)[cols]++;
			}
		}
	}
}
//计算特征值
void feature_computer(Mat& src, double& Asm, double& Eng, double& Con, double& Idm, double& Rel)
{
	int height = src.rows;
	int width = src.cols;
	int total = 0;
	for (int i = 0; i < height; i++)
	{
		int* srcdata = src.ptr<int>(i);
		for (int j = 0; j < width; j++)
		{
			total += srcdata[j];//求图像所有像素的灰度值的和
		}
	}
	Mat copy;
	copy.create(height, width, CV_64FC1);
	for (int i = 0; i < height; i++)
	{
		int* srcdata = src.ptr<int>(i);
		double* copydata = copy.ptr<double>(i);
		for (int j = 0; j < width; j++)
		{
			copydata[j] = (double)srcdata[j] / (double)total;//图像每一个像素的的值除以像素总和
		}
	}
	for (int i = 0; i < height; i++)
	{
		double* srcdata = copy.ptr<double>(i);
		for (int j = 0; j < width; j++)
		{
			Asm += srcdata[j] * srcdata[j];//能量
			if (srcdata[j] > 0)
				Eng -= srcdata[j] * log(srcdata[j]);//熵             
			Con += (double)(i - j) * (double)(i - j) * srcdata[j];//对比度
			Idm += srcdata[j] / (1 + (double)(i - j) * (double)(i - j));//逆差矩
		}
	}
	double u1 = 0, u2 = 0, delta1 = 0, delta2 = 0, temp = 0;
	for (int i = 0; i < height; i++)
	{
		double* srcdata = copy.ptr<double>(i);
		temp = 0;
		for (int j = 0; j < width; j++)
		{
			temp += srcdata[j];
		}
		u1 += temp * i;
	}
	for (int j = 0; j < width; j++)
	{
		double* srcdata = copy.ptr<double>(j);
		temp = 0;
		for (int i = 0; i < height; i++)
		{
			temp += srcdata[i];
		}
		u2 += temp * j;
	}
	int n = 0;
	for (int i = 0; i < height; i++)
	{
		double* srcdata = copy.ptr<double>(i);
		temp = 0;
		for (int j = 0; j < width; j++)
		{
			temp += srcdata[j];
			n += i * j * srcdata[j];
		}
		delta1 += temp * (i - u1) * (i - u1);
	}
	for (int j = 0; j < width; j++)
	{
		double* srcdata = copy.ptr<double>(j);
		temp = 0;
		for (int i = 0; i < height; i++)
		{
			temp += srcdata[i];
		}
		delta2 += temp * (j - u2) * (j - u2);
	}
	Rel = (n - u1 * u2) / (delta1 * delta2);//相关
}

