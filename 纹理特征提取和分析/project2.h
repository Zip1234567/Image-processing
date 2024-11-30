#include<stdio.h>
#include <iostream>
#include <opencv2\opencv.hpp>
using namespace cv;
using namespace std;

//�Ҷȵȼ���Ϊ16
const int gray_level = 16;

//0�ȻҶȹ�������
void getglcm_horison(Mat& input, Mat& dst)
{
	Mat src = input;
	CV_Assert(1 == src.channels());
	src.convertTo(src, CV_32S);
	int height = src.rows;
	int width = src.cols;
	int max_gray_level = 0;
	for (int j = 0; j < height; j++)//Ѱ�����ػҶ����ֵ
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
	max_gray_level++;//���ػҶ����ֵ��1��Ϊ�þ�����ӵ�еĻҶȼ���
	if (max_gray_level > 16)//���Ҷȼ�������16����ͼ��ĻҶȼ���С��16������С�Ҷȹ�������Ĵ�С��
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
	else//���Ҷȼ���С��16����������Ӧ�ĻҶȹ�������
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
//90�ȻҶȹ�������
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
		for (int i = 0; i < height; i++)//��ͼ��ĻҶȼ���С��16������С�Ҷȹ�������Ĵ�С��
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
//45�ȻҶȹ�������
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
		for (int i = 0; i < height; i++)//��ͼ��ĻҶȼ���С��16������С�Ҷȹ�������Ĵ�С��
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
//135�ȻҶȹ�������
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
		for (int i = 0; i < height; i++)//��ͼ��ĻҶȼ���С��16������С�Ҷȹ�������Ĵ�С��
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
//��������ֵ
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
			total += srcdata[j];//��ͼ���������صĻҶ�ֵ�ĺ�
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
			copydata[j] = (double)srcdata[j] / (double)total;//ͼ��ÿһ�����صĵ�ֵ���������ܺ�
		}
	}
	for (int i = 0; i < height; i++)
	{
		double* srcdata = copy.ptr<double>(i);
		for (int j = 0; j < width; j++)
		{
			Asm += srcdata[j] * srcdata[j];//����
			if (srcdata[j] > 0)
				Eng -= srcdata[j] * log(srcdata[j]);//��             
			Con += (double)(i - j) * (double)(i - j) * srcdata[j];//�Աȶ�
			Idm += srcdata[j] / (1 + (double)(i - j) * (double)(i - j));//����
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
	Rel = (n - u1 * u2) / (delta1 * delta2);//���
}

