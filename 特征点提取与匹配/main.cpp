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

	int win_size2 = 15;//��ֵ�����ƴ���
	float threshold = 3.0/200.0;
	int k = floor(window_size / 2);
	
	for (int h = window_size / 2; h < height - window_size / 2; h++) {
		for (int w = window_size / 2; w < width - window_size / 2; w++) {
			// ������Ȥֵ
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
	//��ѡ��
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

	//��ֵ��
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
	int win_size2 = 7;//��ֵ�����ƴ���
	int height = image.rows;
	int width = image.cols;
	Mat corners_w = Mat::zeros(height, width, CV_32FC1);
	Mat corners_q = Mat::zeros(height, width, CV_32FC1);

	//Robert�ݶ�
	Mat Ix, Iy;
	Mat kernelX = (Mat_<float>(2, 2) << 1, 0, 0, -1);
	Mat kernelY = (Mat_<float>(2, 2) << 0, 1, -1, 0);

	filter2D(image, Ix, CV_32F, kernelX);
	filter2D(image, Iy, CV_32F, kernelY);
	//����Э�������ȨֵԲ��
	for (int h = window_size / 2; h < height - window_size / 2; h++) {
		for (int w = window_size / 2; w < width - window_size / 2; w++) {
			// ��ȡ����
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

	//��ѡ��
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

	//��ֵ��
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



// ����ͼ���֮��Ĺ�һ�����ϵ�� (NCC)
double calculateNCC(const Mat& blockA, const Mat& blockB) {
	// ȷ���������Сһ��
	if (blockA.size() != blockB.size()) {
		return -1;
	}

	// ��ֵ����
	Scalar meanA = mean(blockA);
	Scalar meanB = mean(blockB);

	// ȥ��ֵ
	Mat A = blockA - meanA;
	Mat B = blockB - meanB;

	// ������� (A��B�ĵ��)
	double numerator = sum(A.mul(B))[0];

	// �����ĸ (A��B�ķ����˻�)
	double denominator = sqrt(sum(A.mul(A))[0]) * sqrt(sum(B.mul(B))[0]);

	// ��ֹ������
	if (denominator < 1e-10) {
		return 0;
	}

	return numerator / denominator;
}

// �ڸ���������Χ���ҵ����ƥ��
Point2f findBestMatch(const Mat& imgL, const Mat& imgR, const Point2f& ptL,const Point2f& initialPtR, const Size& windowSize, int searchRadius, double& bestNCC) {
	// ��ȡ��ͼ��ͼ���
	Rect windowL(ptL.x - windowSize.width / 2, ptL.y - windowSize.height / 2, windowSize.width, windowSize.height);
	windowL &= Rect(0, 0, imgL.cols, imgL.rows);  // ���ƴ�����ͼ��Χ��
	Mat blockL = imgL(windowL);

	// ������Χ
	bestNCC = -1;  // ��ʼΪ���ƥ��
	Point2f bestMatch = initialPtR;

	for (int dx = -searchRadius; dx <= searchRadius; ++dx) {
		for (int dy = -searchRadius; dy <= searchRadius; ++dy) {
			Point2f ptR = initialPtR + Point2f(dx, dy);

			// ��ȡ��ͼ��ͼ���
			Rect windowR(ptR.x - windowSize.width / 2, ptR.y - windowSize.height / 2, windowSize.width, windowSize.height);
			windowR &= Rect(0, 0, imgR.cols, imgR.rows);  // ���ƴ�����ͼ��Χ��

			// ������ڳ���ͼ��Χ������
			if (windowR.width != windowSize.width || windowR.height != windowSize.height) {
				continue;
			}

			Mat blockR = imgR(windowR);

			// �������ϵ��
			double ncc = calculateNCC(blockL, blockR);

			// Ѱ��������ϵ����ƥ���
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

	// ����ͼ��ȡ�����㣨Moravec/Forstner�ǵ��⣩
	int window_size = 5;  // ����������Ȥֵ���ڴ�С
	//vector<KeyPoint> keypointsL = Forstner_corner_detection(imgL, window_size);
	vector<KeyPoint> keypointsL = Moravec_corner_detection(imgL, window_size);
	cout << "Number of keypoints in left image: " << keypointsL.size() << endl;

	Mat Image_Kp;
	imgL.convertTo(Image_Kp, CV_8UC1, 255);
	drawKeypoints(Image_Kp, keypointsL, Image_Kp, Scalar(255, 0, 0), DrawMatchesFlags::DEFAULT);

	namedWindow("imageWithKeypoints",WINDOW_NORMAL);
	imshow("imageWithKeypoints", Image_Kp);
	waitKey(0);

	// ����ͼ�������
	vector<Mat> pyramidL, pyramidR;
	Mat currentL = imgL, currentR = imgR;

	for (int i = 0; i < 5; i++) {
		pyramidL.push_back(currentL);
		pyramidR.push_back(currentR);
		pyrDown(currentL, currentL);
		pyrDown(currentR, currentR);
	}

	// ����ͷֱ��ʲ���г���ƥ��
	vector<DMatch> matches1to2;
	vector<KeyPoint> keypointsR;  // ������ͼ�������㼯��
	Size windowSize(7, 7);  // ƥ�䴰�ڴ�С
	int searchRadius = 20;     // ������Χ
	double ncc = 0;

	for (size_t i = 0; i < keypointsL.size(); i++) {
		Point2f ptL = keypointsL[i].pt;  // ��ͼ������
		Point2f ptR = findBestMatch(pyramidL.back(), pyramidR.back(), ptL / (1 << (pyramidL.size() - 1)), ptL / (1 << (pyramidL.size() - 1)), windowSize, searchRadius,ncc);
		
		float scaleFactor = 1.0f / (1 << (pyramidL.size() - 1));

		//�������
		for (int level = pyramidL.size() - 2; level >= 0; --level) {
			// ��������������
			scaleFactor *= 2.0f;    // �Ŵ�һ����ptL����С�ֱ��ʲ��
			Point2f scaledPtL = ptL * scaleFactor;//*1/16  *1/8   *1/4
			Point2f scaledPtR = ptR * 2;//�Ŵ�һ��

			// �ڵ�ǰ�����ƥ�䣬�ҵ�����ȷ������
			ptR = findBestMatch(pyramidL[level],pyramidR[level],scaledPtL, scaledPtR,windowSize,searchRadius / 2,ncc);
		}
		
		// ƥ���λ��
		if (ncc > 0.96) {
			KeyPoint kpR(ptR, 1.0f);
			keypointsR.push_back(kpR);
			matches1to2.push_back(DMatch(i, static_cast<int>(keypointsR.size()) - 1, 0));
		}
	}

	// ��ʾ
	Mat imgMatches;
	Mat imgL_8U, imgR_8U;
	imgL.convertTo(imgL_8U, CV_8UC1, 255);// ����255��Ϊ�˽����������ŵ�0-255��Χ
	imgR.convertTo(imgR_8U, CV_8UC1, 255);
	drawMatches(imgL_8U, keypointsL, imgR_8U, keypointsR, matches1to2, imgMatches, Scalar(0, 255, 0), Scalar(0, 255, 0), std::vector<char>(),DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS);

	// ��ʾ
	namedWindow("Matches", WINDOW_NORMAL);
	imshow("Matches", imgMatches);
	waitKey(0);

	return 0;
}

//ѡ����Ҫ������ʽ�����ݣ�����CTRL + A����Ȼ���ٰ�Ctrl + K ��Ctrl + F �ͺ���
//#include<opencv2/xfeatures2.hpp>