// dip2023.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <opencv2/core/core.hpp>   
#include <opencv2/highgui/highgui.hpp>


using namespace std;  // 省去屏幕输出函数cout前的std::
using namespace cv;   // 省去opencv函数前面加cv::的必要性

Mat pre_process_color(Mat M)   //预处理返回矩阵
{
    //landsat的bgr分别是234波段，全色是第8波段
    Mat mss = M;   // 读入图片 

    cout << "depth = " << mss.depth() << endl << "channels = " << mss.channels() << endl;
    cout << "Number of rows = " << mss.rows << endl << "Number of columns = " << mss.cols << endl;
    cout << "Dimension = " << mss.dims << endl << "Number of bytes per element = " << mss.elemSize() << endl;
    cout << "Number of bytes per channel per element = " << mss.elemSize1() << endl << "type = " << mss.type() << endl;

    
    Mat b_g_r_ir[4], color_arr[3], color;
    split(mss, b_g_r_ir); //拆分通道

    color_arr[0] = b_g_r_ir[2];
    color_arr[1] = b_g_r_ir[1];
    color_arr[2] = b_g_r_ir[0];
    merge(color_arr, 3, color);//合并color

    //各波段数值范围大概在0-1001之间，所以建议转至8比特时，除以4
    //color.convertTo(color, CV_8UC3, 1.0 / 4.0, 0);//16位2字节/2
	
	//若仅需要裁剪部分做后续处理
	//Mat subset = color(Range(2650, 3150), Range(5400, 5900)); //subset共享color的部分数据，前面为行范围，后面为列范围
	//imwrite("mss_subset.bmp", subset); //16U的数据可以保存至png、jpeg2000或tiff文件，其中png和jpeg2000需要设置压缩参数前面为行范围，后面为列范围
 
    return color;

}
//C: / Users / 18440 / Desktop / 分发学生 / 实习影像 / mss.tif

Mat pre_process_pan(Mat M)
{
    Mat pan = M;    //读入全色波段

    cout << "depth = " << pan.depth() << endl << "channels = " << pan.channels() << endl;
    cout << "Number of rows = " << pan.rows << endl << "Number of columns = " << pan.cols << endl;
    cout << "Dimension = " << pan.dims << endl << "Number of bytes per element = " << pan.elemSize() << endl;
    cout << "Number of bytes per channel per element = " << pan.elemSize1() << endl << "type = " << pan.type() << endl;

    //最大值552， 除以4
    pan.convertTo(pan, CV_8UC1, 1.0 / 4.0, 0);
	

    return pan;

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
