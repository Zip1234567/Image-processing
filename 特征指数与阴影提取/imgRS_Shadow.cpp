// imgRS_Shadow.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//cpp


#include<windows.h>
#include <iostream>
#include"opencv2/opencv.hpp"
using namespace std;
#pragma comment(lib,"opencv_world480d.lib")
using namespace cv;

Mat Shadow(Mat Mhsv) {
    float* pImData = (float*)Mhsv.data;
    int nHeight = Mhsv.rows;//获取图像的高
    int nWidth = Mhsv.cols;//获取图像的宽
    int nChannels = Mhsv.channels();//获取图像通道数目 BGR/单通道
    Mat M = Mat::zeros(nHeight, nWidth, CV_32FC1);
    float* pImDataI = (float*)M.data;
    for (int i = 0; i < nHeight; i++) {
        for (int j = 0; j < nWidth; j++) {
            float temph = *(pImData + (nWidth * i + j) * nChannels);
            float temps = *(pImData + (nWidth * i + j) * nChannels + 1);
            float tempv = *(pImData + (nWidth * i + j) * nChannels + 2);
            temph = temph * 255.0 / 179;
            temps = temps * 255;
            tempv = tempv * 255;
            float temp3 = (float)(temph - tempv) / (float)(temph + temps + tempv);
            *(pImDataI + i * nWidth + j) = temp3;
        }
    }
    return M;
}

Mat C1C2C3(Mat Img) {
    Mat HSVMat;
    int nHeight = Img.rows;//获取图像的高
    int nWidth = Img.cols;//获取图像的宽
    int nChannels = Img.channels();//获取图像通道数目 BGR/单通道
    double* pImData = (double*)Img.data;
    Mat Mnew = Mat::zeros(nHeight, nWidth, CV_8UC1);
    unsigned char* pImDataNew = Mnew.data;
    for (int i = 0; i < nHeight; i++) {
        for (int j = 0; j < nWidth; j++) {
            double B = *(pImData + (nWidth * i + j) * nChannels);
            double G = *(pImData + (nWidth * i + j) * nChannels + 1);
            double R = *(pImData + (nWidth * i + j) * nChannels + 2);
            double C1 = atan(R / max(G, B));
            double C2 = atan(G / max(R, B));
            double C3 = atan(B / max(G, R));
            if (C3 > 0.35 && B < 34) {
                *(pImDataNew + (nWidth * i + j)) = 255;
            }
            else {
                *(pImDataNew + (nWidth * i + j)) = 0;
            }
        }
    }
    return Mnew;
}

Mat LineTransform(Mat M, float a, float b) {
    double* pImData = (double*)M.data;
    int nHeight = M.rows;//获取图像的高
    int nWidth = M.cols;//获取图像的宽
    int nChannels = M.channels();//获取图像通道数目 BGR
    Mat Mnew = Mat::zeros(nHeight, nWidth, CV_64FC1);
    double* pImDataNew = (double*)Mnew.data;
    for (int k = 0; k < nChannels; k++) {//统计出原图像的最大最小灰度级
        double nMingray = 255, nMaxgray = 0;
        for (int i = 0; i < nHeight; i++) {
            for (int j = 0; j < nWidth; j++) {
                double temp = *(pImData + (nWidth * i + j) * nChannels + k);
                if (temp < nMingray) {
                    nMingray = temp;
                }
                if (temp > nMaxgray) {
                    nMaxgray = temp;
                }
            }
        }
        for (int i = 0; i < nHeight; i++) {
            for (int j = 0; j < nWidth; j++) {
                double temp = *(pImData + (nWidth * i + j) * nChannels + k);
                temp = (double)(a + (b - a) / (nMaxgray - nMingray) * (temp - nMingray));

                *(pImDataNew + (nWidth * i + j) * nChannels + k) = temp;
            }
        }
    }
    return Mnew;
}


int main() {
    string strFileName = "C:\\Users\\18440\\Desktop\\Opencv\\imgRS_Shadow\\Color.bmp";
    Mat Mcolor = imread(strFileName, IMREAD_ANYCOLOR);
    if (Mcolor.empty()) {
        cout << "error";
        return 0;
    }
    namedWindow("原图", WINDOW_NORMAL);
    imshow("原图", Mcolor);
    waitKey();
    
    //HSV变换 hsv cvtcolor只能输入8U16U32F
    //Mcolor.convertTo(Mcolor, CV_32F, 1. / 255);//unchar2CV_32F
    //Mat Mhsv = Mat::zeros(Mcolor.rows, Mcolor.cols, CV_32F);
    //cvtColor(Mcolor, Mhsv, COLOR_BGR2HSV);
    //namedWindow("hsv图", WINDOW_NORMAL);
    //imshow("hsv图", Mhsv);
    //waitKey();



    //计算指数
    Mat Mnew = Mat::zeros(Mcolor.rows, Mcolor.cols, CV_32FC1);
    //Mnew = Shadow(Mhsv);

    Mcolor.convertTo(Mcolor, CV_64F);
    Mnew = C1C2C3(Mcolor);
    
    namedWindow("计算C1C2C3指数图", WINDOW_NORMAL);
    imshow("计算C1C2C3指数图", Mnew);
    waitKey();

    Mnew.convertTo(Mnew, CV_64F);

    //线性拉伸
    Mat Mout = Mat::zeros(Mnew.rows, Mnew.cols, CV_64FC1);
    Mout = LineTransform(Mnew, 0, 255);
    Mout.convertTo(Mout, CV_8UC1);
    namedWindow("线性拉伸图", WINDOW_NORMAL);
    imshow("线性拉伸图", Mout);
    //imwrite("阴影二值化前.jpg", Mout);
    waitKey();

    //阈值分割 迭代 OTUS 最小误差
    //阈值化
    threshold(Mout, Mout, 195, 255, THRESH_BINARY);
    //OSTU法
    //threshold(Mout, Mout, 0, 255, THRESH_OTSU);
    //自适应不同亮暗区域
    //adaptiveThreshold(Mout, Mout, 255, ADAPTIVE_THRESH_MEAN_C, THRESH_BINARY, 11, 2);

    namedWindow("二值化图", WINDOW_NORMAL);
    imshow("二值化图", Mout);
    imwrite("C二值化.jpg", Mout);
    waitKey();



    //Mat element1 = cv::getStructuringElement(
    //    cv::MORPH_ELLIPSE, cv::Size(7, 7));
    //// 腐蚀膨胀操作
    //dilate(Mout, Mout, element1);


    // 腐蚀膨胀操作
    Mat element2 = cv::getStructuringElement(
        cv::MORPH_ELLIPSE, cv::Size(5, 5));
    morphologyEx(Mout, Mout, MORPH_CLOSE, element2);
    Mat element3 = cv::getStructuringElement(
        cv::MORPH_ELLIPSE, cv::Size(3, 3));
    morphologyEx(Mout, Mout, MORPH_OPEN, element3);
    //erode(Mout, Mout, element3);
    namedWindow("形态学处理图", WINDOW_NORMAL);
    imshow("形态学处理图", Mout);
    imwrite("C形态学.jpg", Mout);
    waitKey();

    return 0;
}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
