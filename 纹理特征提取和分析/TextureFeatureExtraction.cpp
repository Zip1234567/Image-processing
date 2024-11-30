// TextureFeatureExtraction.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//灰度共生矩阵(GLDM)的统计方法
//熵
//逆差矩
//主函数运行
//彩色 现在已有每个小窗口的5个特征值，我需要提取林地的特征值，分类
//位置i,j 五个特征值
//
//

#include <iostream>
#include "windows.h"
#include "opencv2\core\core.hpp"
#include "opencv2\highgui\highgui.hpp"
#include "dip2023.h"
#include "project2.h"
using namespace std;
using namespace cv;

vector<vector<double>> GLDM(int& TH, int& TW, Mat M);//生成GLDM每个元素对应的特征值



int main()
{
//读取文件：    
B:
    cout << "请输入图像文件名\n";
    string FN;
    cin >> FN;
    cout << "输入图像文件名为：" << FN << "\n";
    Mat M0 = imread(FN, -1);
    if (M0.empty())
    {
        cout << "输入文件不存在" << "\n";
        goto B;
    }
    //1 16位转8位
    cout << "\n是否将图像进行16位到8位的转换，是输入0，否输入-1：";
    int t;
    cin >> t;
    if (t == 0) {
        if (M0.channels() != 1) M0 = pre_process_color(M0);
        if (M0.channels() == 1) M0 = pre_process_pan(M0);
    }
    //初步处理成黑白图像M
    Mat M;
    cvtColor(M0, M, COLOR_BGR2GRAY);//加权公式gray=0.299×R+0.587×G+0.114×B
    namedWindow("图像", WINDOW_NORMAL);  //若图片太大，用WINDOW_NORMAL的方式可以让图片的显示大小随窗口大小缩放
    imshow("图像", M);  // 显示图片 
    waitKey();
    //2   统计灰度矩阵并归一化
    //分割图像设置窗口大小奇数
    int TH, TW = 0;
    cout << "请输入模板高度：";
    cin >> TH;
    cout << "请输入模板长度：";
    cin >> TW;
    //每个像素对应邻域特征值提取到feature

    vector<vector<double>> feature = GLDM(TH, TW, M);
    //生成纹理特征图
    int height = M0.rows;
    int width = M0.cols;
    int i, j, m, n;
    Mat Texture = M.clone();
    int q = 0;
    // feature1生成纹理特征图
    for (i = TH / 2; i < height - TH / 2; i++)
    {
        for (j = TW / 2; j < width - TW / 2; j++)
        {
            Texture.at<uchar>(i, j) = feature[q][1]*255;
            q++;
        }
    }
    namedWindow("纹理特征图", WINDOW_NORMAL);  
    imshow("纹理特征图", Texture);  // 显示图片 
    waitKey();
    imwrite("C:/Users/18440/Desktop/纹理特征图.jpg", Texture);
    
    
    //确定阈值进行分类标示
    int k = 0;   
    for (i = TH / 2; i < height - TH / 2; i++)
    {
        for (j = TW / 2; j < width - TW / 2; j++)
        {
            if (feature[k][0] * 100<20&& feature[k][3] * 130 < 73) {
                M0.at<Vec3b>(i, j)=Vec3b(255,0,0);
            }
            k++;
        }
    }
    namedWindow("", WINDOW_NORMAL);  
    imshow("林地", M0);  // 显示图片 
    waitKey();
    imwrite("C:/Users/18440/Desktop/结果图.jpg", M0);
    feature.clear();
    return 0;
}

vector<vector<double>> GLDM(int& TH, int& TW, Mat M)
{
    int height = M.rows;
    int width = M.cols;
    int i, j, m, n;
    //int b = (height / TH) * (width / TW);
    vector<vector<double>> feature;
    vector<double> B;
    Mat block = Mat::zeros(TH, TW, CV_8UC1);
    Mat dst_horison;
    for (i = TH / 2; i < height - TH / 2; i++)
    {
            for (j = TW / 2; j < width - TW / 2; j++)//窗口中心移动遍历图像给block赋值
            {               
                for (m = 0; m < TH; m++)
                {
                    for (n = 0; n < TW; n++)//卷积
                    {
                        block.at<uchar>(m, n) =  M.at<uchar>(i - TH / 2 + m, j  - TW / 2 + n) ;
                    }
                }
                //block的灰度共生矩阵0度
                getglcm_horison(block, dst_horison);
                //统计block的特征值
                double asm_horison = 0, eng_horison = 0, con_horison = 0, idm_horison = 0,  rel_horison = 0;
                feature_computer(dst_horison, asm_horison, eng_horison, con_horison, idm_horison, rel_horison);//block相当于8位小图
                //填入容器
                B.push_back(asm_horison);
                B.push_back(eng_horison);
                B.push_back(con_horison);
                B.push_back(idm_horison);
                B.push_back(rel_horison);
                //90
                getglcm_vertical(block, dst_horison);
                double eng_horison90 = 0, con_horison90 = 0, idm_horison90 = 0, asm_horison90 = 0, rel_horison90 = 0;
                feature_computer(dst_horison, asm_horison90, eng_horison90, con_horison90, idm_horison90, rel_horison90);//block相当于8位小图
                B.push_back(asm_horison90);
                B.push_back(eng_horison90);
                B.push_back(con_horison90);
                B.push_back(idm_horison90);
                B.push_back(rel_horison90);
                //45
                getglcm_45(block, dst_horison);
                double eng_horison45 = 0, con_horison45 = 0, idm_horison45 = 0, asm_horison45 = 0, rel_horison45 = 0;
                feature_computer(dst_horison, asm_horison45, eng_horison45, con_horison45, idm_horison45, rel_horison45);//block相当于8位小图
                B.push_back(asm_horison45);
                B.push_back(eng_horison45);
                B.push_back(con_horison45);
                B.push_back(idm_horison45);
                B.push_back(rel_horison45);
                //135
                getglcm_135(block, dst_horison);
                double eng_horison135 = 0, con_horison135 = 0, idm_horison135 = 0, asm_horison135 = 0, rel_horison135 = 0;
                feature_computer(dst_horison, asm_horison135, eng_horison135, con_horison135, idm_horison135, rel_horison135);//block相当于8位小图
                B.push_back(asm_horison135);
                B.push_back(eng_horison135);
                B.push_back(con_horison135);
                B.push_back(idm_horison135);
                B.push_back(rel_horison135);          
                feature.push_back(B);
                B.clear();                
            }
    }
    return feature;
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

//feature[k][0] * 100<20&& feature[k][3] * 130 < 73