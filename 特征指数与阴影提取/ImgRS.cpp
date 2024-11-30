// ImgRS.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//项目 属性 目录 库目录 包含目录

#include<windows.h>
#include <iostream>
#include"opencv2/opencv.hpp"
using namespace std;
#pragma comment(lib,"opencv_world480d.lib")
using namespace cv;


//读入数据
//转换数据
//

//各个函数计算
Mat RVI(Mat ImgR, Mat ImgNIR) {
    double* pImDataR = (double*)ImgR.data;
    int nHeight = ImgR.rows;//获取图像的高
    int nWidth = ImgR.cols;//获取图像的宽
    int nChannels = ImgR.channels();//获取图像通道数目 BGR/单通道
    double* pImDataNir = (double*)ImgNIR.data; 
    Mat M = Mat::zeros(nHeight, nWidth, CV_64FC1);
    double* pImDataRVI = (double*)M.data;
    for (int i = 0; i < nHeight; i++) {
        for (int j = 0; j < nWidth; j++) {
            double temp1 = *(pImDataR + i * nWidth + j);
            double temp2 = *(pImDataNir + i * nWidth + j);
            double temp3 = (double)temp2 / temp1;
            *(pImDataRVI + i * nWidth + j) = temp3;
        }
    }
    return M;
}

Mat NDVI(Mat ImgR, Mat ImgNIR) {
    double* pImDataR = (double*)ImgR.data;
    int nHeight = ImgR.rows;//获取图像的高
    int nWidth = ImgR.cols;//获取图像的宽
    int nChannels = ImgR.channels();//获取图像通道数目 BGR/单通道
    double* pImDataNir = (double*)ImgNIR.data; //
    Mat M = Mat::zeros(nHeight, nWidth, CV_64FC1);
    double* pImDataRVI = (double*)M.data;
    for (int i = 0; i < nHeight; i++) {
        for (int j = 0; j < nWidth; j++) {
            double temp1 = *(pImDataR + i * nWidth + j);
            double temp2 = *(pImDataNir + i * nWidth + j);
            double temp3 = (double)(temp2 - temp1) / (double)(temp1 + temp2);
            *(pImDataRVI + i * nWidth + j) = temp3;
        }
    }
    return M;
}

Mat SAVI(Mat ImgR, Mat ImgNIR) {
    double* pImDataR = (double*)ImgR.data;
    int nHeight = ImgR.rows;//获取图像的高
    int nWidth = ImgR.cols;//获取图像的宽
    int nChannels = ImgR.channels();//获取图像通道数目 BGR/单通道
    double* pImDataNir = (double*)ImgNIR.data; 
    Mat M = Mat::zeros(nHeight, nWidth, CV_64FC1);
    double* pImDataRVI = (double*)M.data;
    for (int i = 0; i < nHeight; i++) {
        for (int j = 0; j < nWidth; j++) {
            double temp1 = *(pImDataR + i * nWidth + j);
            double temp2 = *(pImDataNir + i * nWidth + j);
            double temp3 = (double)(temp2 - temp1) / (double)(temp1 + temp2 + 0.5) * (1 + 0.5);
            *(pImDataRVI + i * nWidth + j) = temp3;
        }
    }
    return M;
}

Mat SAVI2(Mat ImgR, Mat ImgNIR) {
    double* pImDataR = (double*)ImgR.data;
    int nHeight = ImgR.rows;//获取图像的高
    int nWidth = ImgR.cols;//获取图像的宽
    int nChannels = ImgR.channels();//获取图像通道数目 BGR/单通道
    double* pImDataNir = (double*)ImgNIR.data;
    Mat M = Mat::zeros(nHeight, nWidth, CV_64FC1);
    double* pImDataRVI = (double*)M.data;
    for (int i = 0; i < nHeight; i++) {
        for (int j = 0; j < nWidth; j++) {
            double temp1 = *(pImDataR + i * nWidth + j);
            double temp2 = *(pImDataNir + i * nWidth + j);
            double temp3 = (double)(2*temp2+1-(sqrt(pow((2*temp2+1),2)-8*(double)(temp2-temp1)))) / 2;
            *(pImDataRVI + i * nWidth + j) = temp3;
        }
    }
    return M;
}

Mat NDWI(Mat ImgR, Mat ImgNIR) {
    double* pImDataR = (double*)ImgR.data;
    int nHeight = ImgR.rows;//获取图像的高
    int nWidth = ImgR.cols;//获取图像的宽
    int nChannels = ImgR.channels();//获取图像通道数目 BGR/单通道
    double* pImDataNir = (double*)ImgNIR.data;
    Mat M = Mat::zeros(nHeight, nWidth, CV_64FC1);
    double* pImDataRVI = (double*)M.data;
    for (int i = 0; i < nHeight; i++) {
        for (int j = 0; j < nWidth; j++) {
            double temp1 = *(pImDataR + i * nWidth + j);
            double temp2 = *(pImDataNir + i * nWidth + j);
            double temp3 = (double)(temp1 - temp2) / (double)(temp1 + temp2 ) ;
            *(pImDataRVI + i * nWidth + j) = temp3;
        }
    }
    return M;
}

Mat EWI(Mat ImgR, Mat ImgNIR, Mat ImgMIR) {
    double* pImDataR = (double*)ImgR.data;
    int nHeight = ImgR.rows;//获取图像的高
    int nWidth = ImgR.cols;//获取图像的宽
    int nChannels = ImgR.channels();//获取图像通道数目 BGR/单通道
    double* pImDataNir = (double*)ImgNIR.data;
    double* pImDataMIR = (double*)ImgMIR.data;
    Mat M = Mat::zeros(nHeight, nWidth, CV_64FC1);
    double* pImDataRVI = (double*)M.data;
    for (int i = 0; i < nHeight; i++) {
        for (int j = 0; j < nWidth; j++) {
            double temp1 = *(pImDataR + i * nWidth + j);
            double temp2 = *(pImDataNir + i * nWidth + j);
            double temp4 = *(pImDataMIR + i * nWidth + j);
            double temp3 = (double)(temp1 - temp2 - temp4) / (double)(temp1 + temp2 + temp4);
            *(pImDataRVI + i * nWidth + j) = temp3;
        }
    }
    return M;
}

Mat CWI(Mat ImgG, Mat ImgR, Mat ImgNIR, Mat ImgMIR) {
    double* pImDataG = (double*)ImgG.data;
    double* pImDataR = (double*)ImgR.data;
    int nHeight = ImgR.rows;//获取图像的高
    int nWidth = ImgR.cols;//获取图像的宽
    int nChannels = ImgR.channels();//获取图像通道数目 BGR/单通道
    double* pImDataNir = (double*)ImgNIR.data;
    double* pImDataMIR = (double*)ImgMIR.data;
    Mat M = Mat::zeros(nHeight, nWidth, CV_64FC1);
    double* pImDataRVI = (double*)M.data;
    for (int i = 0; i < nHeight; i++) {
        for (int j = 0; j < nWidth; j++) {
            double temp0 = *(pImDataG + i * nWidth + j);
            double temp1 = *(pImDataR + i * nWidth + j);
            double temp2 = *(pImDataNir + i * nWidth + j);
            double temp4 = *(pImDataMIR + i * nWidth + j);
            double L = 1.507599 * temp0 - 0.066392 * temp1 - 1.382209 * temp2 + 1.733790 * temp4 + 11;
            double B = 1.126971 * temp0 + 0.673348 * temp1 + 0.077966 * temp2 - 1.878287 * temp4 + 159;
            double V = 1.636910 * temp0 - 3.396809 * temp1 + 1.915944 * temp2 - 0.156048 * temp4 + 121;
            double temp3 = (double)(temp1 - temp2 - temp4) / (double)(temp1 + temp2 + temp4);
            *(pImDataRVI + i * nWidth + j) = temp3;
        }
    }
    return M;
}

Mat DBI(Mat ImgNIR, Mat ImgFIR) {
    double* pImDataNIR = (double*)ImgNIR.data;
    int nHeight = ImgNIR.rows;//获取图像的高
    int nWidth = ImgNIR.cols;//获取图像的宽
    int nChannels = ImgNIR.channels();//获取图像通道数目 BGR/单通道
    double* pImDataFIR = (double*)ImgFIR.data;
    Mat M = Mat::zeros(nHeight, nWidth, CV_64FC1);
    double* pImDataRVI = (double*)M.data;
    for (int i = 0; i < nHeight; i++) {
        for (int j = 0; j < nWidth; j++) {
            double temp1 = *(pImDataNIR + i * nWidth + j);
            double temp2 = *(pImDataFIR + i * nWidth + j);
            double temp3 = (double)(temp2- temp1);
            *(pImDataRVI + i * nWidth + j) = temp3;
        }
    }
    return M;
}

Mat NDBI(Mat ImgNIR, Mat ImgMIR) {
    double* pImDataNIR = (double*)ImgNIR.data;
    int nHeight = ImgNIR.rows;//获取图像的高
    int nWidth = ImgNIR.cols;//获取图像的宽
    int nChannels = ImgNIR.channels();//获取图像通道数目 BGR/单通道
    double* pImDataMIR = (double*)ImgMIR.data;
    Mat M = Mat::zeros(nHeight, nWidth, CV_64FC1);
    double* pImDataRVI = (double*)M.data;
    for (int i = 0; i < nHeight; i++) {
        for (int j = 0; j < nWidth; j++) {
            double temp1 = *(pImDataNIR + i * nWidth + j);
            double temp2 = *(pImDataMIR + i * nWidth + j);
            double temp3 = (double)(temp2 - temp1)/ (double)(temp2 + temp1);
            *(pImDataRVI + i * nWidth + j) = temp3;
        }
    }
    return M;
}


//线性拉伸
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
    Mat M[7];
    for (int i = 0; i < 7; i++) {
        string strFileName = "C:\\Users\\18440\\Desktop\\Opencv\\ImgRS\\second_data\\tm" + to_string(i + 1) + ".tif";
        M[i] = imread(strFileName, IMREAD_ANYCOLOR);
    }

    for (int i = 0; i < 7; i++) {
        M[i].convertTo(M[i], CV_64FC1);
    }

    //计算指数
    Mat Mnew = Mat::zeros(M[1].rows, M[1].cols, CV_64FC1);
    //Mnew = RVI(M[2], M[3]);
    //Mnew = NDVI(M[2], M[3]);
    //Mnew = SAVI(M[2], M[3]);
    //Mnew = SAVI2(M[2], M[3]);
    
    //Mnew = NDWI(M[1], M[3]);
    //Mnew = NDWI(M[1], M[4]);//MNDWI 使用短波红外波段(4)替换NIR(3)
    //Mnew = EWI(M[1], M[3], M[4]);
  
    //Mnew = DBI(M[3], M[6]);
    Mnew = NDBI(M[3], M[4]);
    //Mnew = NDBI(M[3], M[4]);


    namedWindow("计算指数图", WINDOW_NORMAL);
    imshow("计算指数图", Mnew);
    waitKey();

    //线性拉伸
     Mat Mout = Mat::zeros(M[1].rows, M[1].cols, CV_64FC1);
     Mout=LineTransform(Mnew, 0, 255);
     Mout.convertTo(Mout, CV_8UC1);
     namedWindow("线性拉伸图", WINDOW_NORMAL);
     imshow("线性拉伸图", Mout);
     imwrite("DBI1.jpg", Mout);
     waitKey();

    //阈值分割 迭代 OTUS 最小误差
    //阈值化
    threshold(Mout, Mout, 180, 255, THRESH_BINARY);
    //OSTU法
    //threshold(Mout, Mout, 0, 255, THRESH_OTSU);

    //自适应不同亮暗区域
    //adaptiveThreshold(Mout, Mout, 255, ADAPTIVE_THRESH_MEAN_C, THRESH_BINARY, 11, 2);

    namedWindow("二值化图", WINDOW_NORMAL);
    imshow("二值化图", Mout);
    imwrite("DBI2.jpg", Mout);
    
    waitKey();

    // 腐蚀膨胀操作
    Mat element2 = cv::getStructuringElement(
        cv::MORPH_ELLIPSE, cv::Size(3, 3));
    morphologyEx(Mout, Mout, MORPH_CLOSE,element2);
    //erode(Mout, Mout, element2);
    //applyColorMap(Mnew, Mout, COLORMAP_JET);
    namedWindow("形态学处理图", WINDOW_NORMAL);
    imshow("形态学处理图", Mout);
    imwrite("DBI3.jpg", Mout);
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
