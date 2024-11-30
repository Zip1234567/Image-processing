#include<windows.h>
#include <iostream>
#include"opencv2/opencv.hpp"
#include<opencv2/xfeatures2d/xfeatures2d.hpp>
#include <vector>

using namespace std;
#pragma comment(lib,"opencv_world480d.lib")
using namespace cv;


vector<KeyPoint> Moravec_corner_detection(Mat image, int window_size) {
	int height = image.rows;
	int width = image.cols;
	Mat corners = Mat::zeros(height, width, CV_32FC1);

	int win_size2 = 15;//极值点抑制窗口
	float threshold = 3.0/200.0;
	int k = floor(window_size / 2);
	
	for (int h = window_size / 2; h < height - window_size / 2; h++) {
		for (int w = window_size / 2; w < width - window_size / 2; w++) {
			// 计算兴趣值
			float V[4] = { 0 };
			for (int i = -k; i < k; i++) {
				V[0] += pow(image.at<float>(h, w + i) - image.at<float>(h, w + i + 1), 2);
				V[1] += pow(image.at<float>(h + i, w + i) - image.at<float>(h + i + 1, w + i + 1), 2);
				V[2] += pow(image.at<float>(h + i, w) - image.at<float>(h + i + 1, w), 2);
				V[3] += pow(image.at<float>(h - i, w + i) - image.at<float>(h - i - 1, w + i + 1), 2);
			}
			float iv = *min_element(V, V + 4);
			corners.at<float>(h, w) = iv;
		}
	}
	//候选点
	int countchoosen = 0;
	for (int r = 0; r < height; r++) {
		for (int c = 0; c < width; c++)
		{
			if (corners.at<float>(r, c) > threshold) {
				countchoosen++;
				continue;
			}
			else
				corners.at<float>(r, c) = 0;
		}
	}
	cout << "Number of Choosen keypoints: " << countchoosen << endl;

	//极值点
	vector<KeyPoint> keypoints;
	Mat windowCorners;
	Point maxIdx;
	double maxVal;
	for (int h = win_size2 / 2; h < height - win_size2 / 2; h = h + win_size2) {
		for (int w = win_size2 / 2; w < width - win_size2 / 2; w = w + win_size2) {
			windowCorners = corners(Rect(w - win_size2 / 2, h - win_size2 / 2, win_size2, win_size2));
			minMaxLoc(windowCorners, NULL, &maxVal, NULL, &maxIdx);
			if (maxVal > 0) {
				Point imgLoc(maxIdx.x + (w - win_size2 / 2), maxIdx.y + (h - win_size2 / 2));
				keypoints.push_back(KeyPoint(imgLoc.x, imgLoc.y, 1));
			}
		}
	}
	cout << "Number of Result Keypoints: " << keypoints.size() << endl; 

	return keypoints;
}

vector<KeyPoint> Forstner_corner_detection(Mat image, int window_size ) {
	int win_size2 = 7;//极值点抑制窗口
	int height = image.rows;
	int width = image.cols;
	Mat corners_w = Mat::zeros(height, width, CV_32FC1);
	Mat corners_q = Mat::zeros(height, width, CV_32FC1);

	//Robert梯度
	Mat Ix, Iy;
	Mat kernelX = (Mat_<float>(2, 2) << 1, 0, 0, -1);
	Mat kernelY = (Mat_<float>(2, 2) << 0, 1, -1, 0);

	filter2D(image, Ix, CV_32F, kernelX);
	filter2D(image, Iy, CV_32F, kernelY);
	//窗口协方差矩阵权值圆度
	for (int h = window_size / 2; h < height - window_size / 2; h++) {
		for (int w = window_size / 2; w < width - window_size / 2; w++) {
			// 提取窗口
			Mat windowIx = Ix(Rect(w - window_size / 2, h - window_size / 2, window_size, window_size));
			Mat windowIy = Iy(Rect(w - window_size / 2, h - window_size / 2, window_size, window_size));

			Mat Ixx, Ixy, Iyy;
			multiply(windowIx, windowIx, Ixx);
			multiply(windowIy, windowIy, Iyy);
			multiply(windowIx, windowIy, Ixy);

			double sumIxx = sum(Ixx)[0];
			double sumIyy = sum(Iyy)[0];
			double sumIxy = sum(Ixy)[0];

			double det = (sumIxx * sumIyy) - (sumIxy * sumIxy);
			double trace = sumIxx + sumIyy;

			if (trace == 0)
				continue;

			double ww = det / trace;
			double q = 4 * det / (trace * trace);

			corners_w.at<float>(h, w) = static_cast<float>(ww);
			corners_q.at<float>(h, w) = static_cast<float>(q);

		}
	}
	
	Scalar meanScalar = mean(corners_w, corners_w > 0);
	float meanValue = meanScalar[0];
	cout << "Mean Value of corners_w: " << meanValue << endl;
	float Tq = 0.9;
	float Tw = 8.0 * meanValue;

	//待选点
	int countchoosen = 0;
	for (int r = 0; r < height; r++) {
		for (int c = 0; c < width; c++)
		{
			if (corners_w.at<float>(r, c) > Tw && corners_q.at<float>(r, c) > Tq) {
				countchoosen++;
				continue;
			}
			else
				corners_w.at<float>(r, c) = 0;
		}
	}
	cout << "Number of Choosen keypoints: " << countchoosen << endl;

	//极值点
	vector<KeyPoint> keypoints;
	Mat windowCorners;
	Point maxIdx;
	double maxVal;
	for (int h = win_size2 / 2; h < height - win_size2 / 2; h=h+ win_size2) {
		for (int w = win_size2 / 2; w < width - win_size2 / 2; w=w+ win_size2) {
			windowCorners = corners_w(Rect(w - win_size2 / 2, h - win_size2 / 2, win_size2, win_size2));
			minMaxLoc(windowCorners, NULL, &maxVal, NULL, &maxIdx);
			if (maxVal > 0) {
				Point imgLoc(maxIdx.x + (w - win_size2 / 2), maxIdx.y + (h - win_size2 / 2));
				keypoints.push_back(KeyPoint(imgLoc.x, imgLoc.y, 1));
			}
		}
	}
	cout << "Number of Result keypoints: " << keypoints.size() << endl;  // 0.9  8   300
	return keypoints;
}



// 计算图像块之间的归一化相关系数 (NCC)
double calculateNCC(const Mat& blockA, const Mat& blockB) {
	// 确保两个块大小一致
	if (blockA.size() != blockB.size()) {
		return -1;
	}

	// 均值计算
	Scalar meanA = mean(blockA);
	Scalar meanB = mean(blockB);

	// 去均值
	Mat A = blockA - meanA;
	Mat B = blockB - meanB;

	// 计算分子 (A与B的点积)
	double numerator = sum(A.mul(B))[0];

	// 计算分母 (A和B的范数乘积)
	double denominator = sqrt(sum(A.mul(A))[0]) * sqrt(sum(B.mul(B))[0]);

	// 防止除以零
	if (denominator < 1e-10) {
		return 0;
	}

	return numerator / denominator;
}

// 在给定搜索范围内找到最佳匹配
Point2f findBestMatch(const Mat& imgL, const Mat& imgR, const Point2f& ptL,const Point2f& initialPtR, const Size& windowSize, int searchRadius, double& bestNCC) {
	// 提取左图的图像块
	Rect windowL(ptL.x - windowSize.width / 2, ptL.y - windowSize.height / 2, windowSize.width, windowSize.height);
	windowL &= Rect(0, 0, imgL.cols, imgL.rows);  // 限制窗口在图像范围内
	Mat blockL = imgL(windowL);

	// 搜索范围
	bestNCC = -1;  // 初始为最差匹配
	Point2f bestMatch = initialPtR;

	for (int dx = -searchRadius; dx <= searchRadius; ++dx) {
		for (int dy = -searchRadius; dy <= searchRadius; ++dy) {
			Point2f ptR = initialPtR + Point2f(dx, dy);

			// 提取右图的图像块
			Rect windowR(ptR.x - windowSize.width / 2, ptR.y - windowSize.height / 2, windowSize.width, windowSize.height);
			windowR &= Rect(0, 0, imgR.cols, imgR.rows);  // 限制窗口在图像范围内

			// 如果窗口超出图像范围，跳过
			if (windowR.width != windowSize.width || windowR.height != windowSize.height) {
				continue;
			}

			Mat blockR = imgR(windowR);

			// 计算相关系数
			double ncc = calculateNCC(blockL, blockR);

			// 寻找最大相关系数的匹配点
			if (ncc > bestNCC) {
				bestNCC = ncc;
				bestMatch = ptR;
			}
		}
	}
	//cout << "Testing ptR: " <<bestMatch << " NCC: " << bestNCC << endl;
	return bestMatch;
}

int main() {
	Mat imgL = imread("left_image.tif", IMREAD_GRAYSCALE);
	Mat imgR = imread("right_image.tif", IMREAD_GRAYSCALE);

	if (imgL.empty() || imgR.empty()) {
		cout << "Error: Cannot load images!" << endl;
		return -1;
	}

	imgL.convertTo(imgL, CV_32FC1, 1.0 / 255);
	imgR.convertTo(imgR, CV_32FC1, 1.0 / 255);

	// 从左图提取特征点（Moravec/Forstner角点检测）
	int window_size = 5;  // 特征点检测兴趣值窗口大小
	//vector<KeyPoint> keypointsL = Forstner_corner_detection(imgL, window_size);
	vector<KeyPoint> keypointsL = Moravec_corner_detection(imgL, window_size);
	cout << "Number of keypoints in left image: " << keypointsL.size() << endl;

	Mat Image_Kp;
	imgL.convertTo(Image_Kp, CV_8UC1, 255);
	drawKeypoints(Image_Kp, keypointsL, Image_Kp, Scalar(255, 0, 0), DrawMatchesFlags::DEFAULT);

	namedWindow("imageWithKeypoints",WINDOW_NORMAL);
	imshow("imageWithKeypoints", Image_Kp);
	waitKey(0);

	// 构建图像金字塔
	vector<Mat> pyramidL, pyramidR;
	Mat currentL = imgL, currentR = imgR;

	for (int i = 0; i < 5; i++) {
		pyramidL.push_back(currentL);
		pyramidR.push_back(currentR);
		pyrDown(currentL, currentL);
		pyrDown(currentR, currentR);
	}

	// 在最低分辨率层进行初步匹配
	vector<DMatch> matches1to2;
	vector<KeyPoint> keypointsR;  // 定义右图的特征点集合
	Size windowSize(7, 7);  // 匹配窗口大小
	int searchRadius = 20;     // 搜索范围
	double ncc = 0;

	for (size_t i = 0; i < keypointsL.size(); i++) {
		Point2f ptL = keypointsL[i].pt;  // 左图特征点
		Point2f ptR = findBestMatch(pyramidL.back(), pyramidR.back(), ptL / (1 << (pyramidL.size() - 1)), ptL / (1 << (pyramidL.size() - 1)), windowSize, searchRadius,ncc);
		
		float scaleFactor = 1.0f / (1 << (pyramidL.size() - 1));

		//逐层搜索
		for (int level = pyramidL.size() - 2; level >= 0; --level) {
			// 缩放特征点坐标
			scaleFactor *= 2.0f;    // 放大一倍，ptL逐步缩小分辨率差距
			Point2f scaledPtL = ptL * scaleFactor;//*1/16  *1/8   *1/4
			Point2f scaledPtR = ptR * 2;//放大一倍

			// 在当前层进行匹配，找到更精确的坐标
			ptR = findBestMatch(pyramidL[level],pyramidR[level],scaledPtL, scaledPtR,windowSize,searchRadius / 2,ncc);
		}
		
		// 匹配点位置
		if (ncc > 0.96) {
			KeyPoint kpR(ptR, 1.0f);
			keypointsR.push_back(kpR);
			matches1to2.push_back(DMatch(i, static_cast<int>(keypointsR.size()) - 1, 0));
		}
	}

	// 显示
	Mat imgMatches;
	Mat imgL_8U, imgR_8U;
	imgL.convertTo(imgL_8U, CV_8UC1, 255);// 乘以255是为了将浮点数缩放到0-255范围
	imgR.convertTo(imgR_8U, CV_8UC1, 255);
	drawMatches(imgL_8U, keypointsL, imgR_8U, keypointsR, matches1to2, imgMatches, Scalar(0, 255, 0), Scalar(0, 255, 0), std::vector<char>(),DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS);

	// 显示
	namedWindow("Matches", WINDOW_NORMAL);
	imshow("Matches", imgMatches);
	waitKey(0);

	return 0;
}

//选中需要调整格式的内容（可用CTRL + A），然后再按Ctrl + K 和Ctrl + F 就好了
//#include<opencv2/xfeatures2.hpp>