// ImageResection1.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//


#include<windows.h>
#include <iostream>
#include"opencv2/opencv.hpp"
#include <vector>  
using namespace std;
#pragma comment(lib,"opencv_world480d.lib")
using namespace cv;


//内方位元素
const double f = 153.24 / 1000;
const double x0 = 0; 
const double y_0 = 0;
const int n = 4, m = 15000; 
const double ms = 1e-3, mphi = 1e-6;

//定义点
struct  Points {
    double x;
    double y;
    double X;
    double Y;
    double Z;
};


int main()
{
    //输入原始数据
    Points p[4];
    p[0].x = -86.15/1000, p[0].y = -68.99/1000, p[0].X = 36589.41, p[0].Y = 25273.32, p[0].Z = 2195.17;
    p[1].x = -53.40/1000, p[1].y = 82.21/1000, p[1].X = 37631.08, p[1].Y = 31324.51, p[1].Z = 728.69;
    p[2].x = -14.78/1000, p[2].y = -76.63/1000, p[2].X = 39100.97, p[2].Y = 24934.98, p[2].Z = 2386.50;
    p[3].x = 10.46/1000, p[3].y = 64.43/1000, p[3].X = 40426.54, p[3].Y = 30319.81, p[3].Z = 757.31;

    //初值
    double Xs, Ys, Zs, H;
    double phi, omega, kappa; 
    Xs = (p[0].X + p[1].X + p[2].X + p[3].X) / n;
    Ys = (p[0].Y + p[1].Y + p[2].Y + p[3].Y) / n;
    Zs  = m * f;
    phi = omega = kappa = 0;
    //结果数据
    int t = 20;//最多迭代次数
    double a1, a2, a3, b1, b2, b3, c1, c2, c3;//旋转矩阵R
    double zeta;//单位权中误差
    Mat QQ (6, 6, CV_64FC1);//精度矩阵
    Mat vv(2*n, 1, CV_64FC1);//改正数矩阵
 
    do { 
    //利用角元素计算旋转矩阵   
    a1 = cos(phi) * cos(kappa) - sin(phi) * sin(omega) * sin(kappa);
    a2 = -cos(phi) * sin(kappa) - sin(phi) * sin(omega) * cos(kappa);
    a3 = -sin(phi) * cos(omega);
    b1 = cos(omega) * sin(kappa);
    b2 = cos(omega) * cos(kappa);
    b3 = -sin(omega);
    c1 = sin(phi) * cos(kappa) + cos(phi) * sin(omega) * sin(kappa);
    c2 = -sin(phi) * sin(kappa) + cos(phi) * sin(omega) * cos(kappa);
    c3 = cos(phi) * cos(omega);
    
    //循环计算改正数至小于限差
    //计算像点坐标的近似值
    double _X, _Y, _Z, X, Y, Z, x, y,x_j,y_j;
    double* A = new double[2 * n * 6];
    double* L = new double[2 * n ];
    for (int i = 0; i < n; i++) {
        //xy近似值
        X = p[i].X;
        Y = p[i].Y;
        Z = p[i].Z;
        x = p[i].x ;
        y = p[i].y ;
        //_X,_Y,_Z!!!xyo与xy0
        _X = a1 * (X - Xs) + b1 * (Y - Ys) + c1 * (Z - Zs);
        _Y = a2 * (X - Xs) + b2 * (Y - Ys) + c2 * (Z - Zs);
        _Z = a3 * (X - Xs) + b3 * (Y - Ys) + c3 * (Z - Zs);
        //共线方程算xy近似值xj yj
        x_j = x0 - f * _X / _Z;
        y_j= y_0 - f * _Y / _Z;
        //系数矩阵元素数组A 
        A[0 + i * 12] = (a1 * f + a3 * (x - x0)) / _Z;
        A[1 + i * 12] = (b1 * f + b3 * (x - x0)) / _Z;
        A[2 + i * 12] = (c1 * f + c3 * (x - x0)) / _Z;
        A[3 + i * 12] = (y - y_0) * sin(omega) - ((x - x0) / f * ((x - x0) * cos(kappa) - (y - y_0) * sin(kappa)) + f * cos(kappa)) * cos(omega);
        A[4 + i * 12] = -f * sin(kappa) - (x - x0) / f * ((x - x0) * sin(kappa) + (y - y_0) * cos(kappa));
        A[5 + i * 12] = y - y_0;
        A[6 + i * 12] = (a2 * f + a3 * (y - y_0)) / _Z;
        A[7 + i * 12] = (b2 * f + b3 * (y - y_0)) / _Z;
        A[8 + i * 12] = (c2 * f + c3 * (y-y_0)) / _Z;
        A[9 + i * 12] = -(x - x0) * sin(omega) - ((y-y_0) / f * ((x - x0) * cos(kappa) - (y-y_0) * sin(kappa)) - f * sin(kappa)) * cos(omega);
        A[10 + i * 12] = -f * cos(kappa) - (y-y_0) / f * ((x - x0) * sin(kappa) + (y-y_0) * cos(kappa));
        A[11 + i * 12] = -(x - x0);
        //常数项矩阵元素数组L
        L[0 + 2 * i] = x - x_j;
        L[1 + 2 * i] = y - y_j;
    }
    //系数矩阵AA  常数项矩阵L
    Mat AA = Mat(2 * n, 6, CV_64FC1, A);
    Mat LL = Mat(2 * n, 1, CV_64FC1, L);
    //法方程式计算
    Mat At = AA.t();
    Mat AtA = At*AA;
    Mat AtA_1 = AtA.inv();
    QQ = AtA_1.clone();
    Mat AtL = At*LL;
    Mat xx = AtA_1*AtL;
    vv = AA*xx-LL;
    Mat vtv = (vv.t())*vv ;
   //计算外方位元素的新值与精度计算
    Xs += xx.at<double>(0, 0);
    Ys += xx.at<double>(1, 0);
    Zs += xx.at<double>(2, 0);
    phi += xx.at<double>(3, 0);
    omega += xx.at<double>(4, 0);
    kappa += xx.at<double>(5, 0);
    zeta = vtv.at<double>(0, 0);
    zeta = sqrt(zeta / (2 * n - 6));
   //限差检查
    if ((xx.at<double>(0, 0) < ms) && (xx.at<double>(1, 0) < ms) && (xx.at<double>(2, 0) < ms) && (xx.at<double>(3, 0) < mphi) && (xx.at<double>(4, 0) < mphi) && (xx.at<double>(5, 0) < mphi))
        break;
    //最多循环次数
    t--;
    if (t == 0)
        printf("结果不收敛！\n\n\n");
    //释放
    AA.release(); LL.release(); At.release(); AtA.release();
    AtA_1.release(); AtL.release(); xx.release(); vtv.release();
    delete[] A, L;
    
    } while (t > 0);

    //输出结果
    if (t > 0) {
        FILE* fp;//建立一个文件操作指针
        fopen_s(&fp, "Output.txt", "w");
        //"C:\Users\18440\Desktop\Photogrammetry\ImageResection1\Output.txt"

        printf("\n程序一共经过了%d次迭代计算\n", 20 - t);
        for (int i = 0; i < n; i++) {
            printf("第%d个点的坐标分别改正为 x1= %lfmm，y1= %lfmm\n", i + 1, (vv.at<double>(i * 2, 0) + p[i].x) * 1000, (vv.at<double>(i * 2 + 1, 0) + p[i].y) * 1000);
            fprintf(fp, "第%d个点的坐标分别改正为 x1= %lfmm，y1= %lfmm\n", i + 1, vv.at<double>(i * 2, 0) + p[i].x, vv.at<double>(i * 2 + 1, 0) + p[i].y);
        }

        printf("\n解算得到外方位元素为：\nXs=%.2lf\nYs=%.2lf\nZs=%.2lf\n", Xs, Ys, Zs);
        fprintf(fp, "\n解算得到外方位元素为：\nXs=%.2lf\nYs=%.2lf\nZs=%.2lf\n", Xs, Ys, Zs);


        printf("Phi=%lf rad\nOmega=%lf rad\nKappa=%lf rad\n", phi, omega, kappa);
        printf("平均比例尺为 1：%f\n\n", Zs / f);
        fprintf(fp, "Phi=%lf rad\nOmega=%lf rad\nKappa=%lf rad\n", phi, omega, kappa);
        fprintf(fp, "\n平均比例尺为 1：%f\n\n", Zs / f);


        printf("R矩阵为:\n %.5lf %.5lf %.5lf\n%.5lf %.5lf %.5lf\n%.5lf %.5lf %.5lf\n\n", a1, a2, a3, b1, b2, b3, c1, c2, c3);
        fprintf(fp, "R矩阵为:\n %.5lf %.5lf %.5lf\n%.5lf %.5lf %.5lf\n%.5lf %.5lf %.5lf\n\n", a1, a2, a3, b1, b2, b3, c1, c2, c3);

        printf("精度：\n单位权中误差为：%lf\n", zeta);
        printf("Xs精度为：%lfm\n", zeta * sqrt(QQ.at<double>(0, 0)));
        printf("Ys精度为：%lfm\n", zeta * sqrt(QQ.at<double>(1, 1)));
        printf("Zs精度为：%lfm\n", zeta * sqrt(QQ.at<double>(2, 2)));
        printf("Phi精度为：%lf rad\n", zeta * sqrt(QQ.at<double>(3, 3)));
        printf("Omega精度为：%lf rad\n", zeta * sqrt(QQ.at<double>(4, 4)));
        printf("Kappa精度为：%lf rad\n\n", zeta * sqrt(QQ.at<double>(5, 5)));
        fprintf(fp, "单位权中误差为：%lf\n\n", zeta);
        fprintf(fp, "Xs精度为：%lfm\n", zeta * sqrt(QQ.at<double>(0, 0)));
        fprintf(fp, "Ys精度为：%lfm\n", zeta * sqrt(QQ.at<double>(1, 1)));
        fprintf(fp, "Zs精度为：%lfm\n", zeta * sqrt(QQ.at<double>(2, 2)));
        fprintf(fp, "Phi精度为：%lf rad\n", zeta * sqrt(QQ.at<double>(3, 3)));
        fprintf(fp, "Omega精度为：%lf rad\n", zeta * sqrt(QQ.at<double>(4, 4)));
        fprintf(fp, "Kappa精度为：%lf rad\n\n", zeta * sqrt(QQ.at<double>(5, 5)));

        fclose(fp);

    }
    
    QQ.release(); vv.release();

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
