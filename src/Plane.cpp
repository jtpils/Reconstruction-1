#include "Plane.h"
#include <iostream>
#include <vector>
using namespace std;

int32_t randomColor() {
	int32_t r = rand() % 255;
	int32_t g = rand() % 255;
	int32_t b = rand() % 255;
	int32_t a = 255;
	a = a << 24;
	r = r << 16;
	g = g << 8;
	return a | r | g | b;
}

void TransformPoint(PointT *point, double angle) {//angleは[rad]
	double x = cos(angle) * point->x - sin(angle) * point->y;
	double y = sin(angle) * point->x + cos(angle) * point->y;
	//cloud_projected->points[i].z = cloud_projected->points[i].z;
	point->x = x;
	point->y = y;
}

int colorType2int(PlaneColor colorType) {
	int32_t color = 255 << 24; // control transparent
	switch (colorType) {
	case (Color_Red):
		color = color | 255 << 16;
		break;
	case Color_Yellow:
		color = color | (255 << 8) | (255 << 16);
		break;
	case Color_Blue:
		color = color | 255;
		break;
	case Color_Green:
		color = color | 255 << 8;
		break;
	case Color_White:
		color = INT32_MAX;
		break;
	case Color_Peach:
		color = color | (255 << 16) | (218 << 16) | 185;
		break;
	case Color_Random:
		color = randomColor();
		break;
	default:
		std::cerr << " color hasnt been decided" << std::endl;
		break;
	}
	return color;
}

//void importTXT(const string path, PointCloudT::Ptr cloud) {
//	ifstream fin(path);
//	if (!fin) {
//		cerr << "Cannot find path " << endl;
//	}
//	PointCloudT::Ptr importCloud(new PointCloudT);
//	string str;
//	string delimiter = ",";
//	while (getline(fin, str)) {
//		string s = str;
//		size_t pos = 0;
//		std::string token;
//		vector<float > cord;
//		while ((pos = s.find(delimiter)) != std::string::npos) {
//			token = s.substr(0, pos);
//			cord.push_back(stof(token));
//			s.erase(0, pos + delimiter.length());
//		}
//		cord.push_back(stof(s));
//		if (cord.size() != 7) {
//			cerr << "error happens when import txt file" << endl;
//		}
//		PointT p;
//		p.x = cord[0];
//		p.y = cord[1];
//		p.z = cord[2];
//		p.a = cord[3];
//		p.r = cord[4];
//		p.g = cord[5];
//		p.b = cord[6];
//		cloud.get()->push_back(p);
//		//cloud.get()->push_back(p);
//	}
//}

/** @brief initialize plane
 */
Plane::Plane()
{
	PointCloudT::Ptr tmp(new PointCloudT);
	this->pointCloud = tmp;
}

/** @brief initialize plane based on raw point and plane vector
 */
Plane::Plane(PointCloudT::Ptr rawPointCloud, Eigen::Vector4d abcd)
{
	PointCloudT::Ptr tmp(new PointCloudT);
	this->pointCloud = tmp;
	pcl::copyPointCloud(*rawPointCloud, *this->pointCloud);
	this->_abcd = abcd;
}
/** @brief initialize plane based on 4 pts a1-a2 => b1-b2
 */
Plane::Plane(PointT a1, PointT a2, PointT b1, PointT b2, float pointPitch, PlaneColor color) {
	PointCloudT::Ptr tmp(new PointCloudT);
	this->pointCloud = tmp;
	generatePlanePointCloud(a1, a2, b1, b2, pointPitch, colorType2int(color));
	updateBoundary();
}

void Plane::setType(PlaneType type) {
	this->_type = type;
}

void Plane::setColor(PlaneColor colorType)
{
	int32_t color = colorType2int(colorType);
	PointCloudT* p = pointCloud.get();
	for (size_t i = 0; i < pointCloud.get()->size(); i++) {
		p->at(i).rgba = color;
	}
}

void Plane::filledPlane(int pointPitch) {
	Eigen::Vector4d planePara = this->_abcd;
	double angle = acos(abs(planePara[0]) / sqrt(planePara[0] * planePara[0] + planePara[1] * planePara[1] + planePara[2] * planePara[2])); //注意！！この方程式は、2平面のなす角度は、0〜90度
	if (planePara[0] / planePara[1] > 0) angle = (-1) * angle;
	Eigen::Matrix4f rotation = Eigen::Matrix4f::Identity();
	rotation(0, 0) = cos(angle); rotation(0, 1) = -sin(angle);
	rotation(1, 0) = sin(angle); rotation(1, 1) = cos(angle);
	pcl::transformPointCloud(*this->pointCloud, *this->pointCloud, rotation);
	PointT proj_min;
	PointT proj_max;
	pcl::getMinMax3D(*this->pointCloud, proj_min, proj_max);
	PointCloudT::Ptr cloud_filled_temp(new PointCloudT);
	int32_t color = colorType2int(Color_White);
	for (int i = 0; i < pointPitch * (proj_max.y - proj_min.y); i++) {
		for (int j = 0; j < pointPitch * (proj_max.z - proj_min.z); j++) {
			PointT pointForFill;
			if (planePara[0] / planePara[1] > 0 && planePara[0] < 0)
				pointForFill.x += planePara[3]; //
			else if (planePara[0] / planePara[1] < 0 && planePara[0] < 0)
				pointForFill.x += planePara[3];
			else if (planePara[0] / planePara[1] > 0 && planePara[0] > 0)
				pointForFill.x -= planePara[3];
			else if (planePara[0] / planePara[1] < 0 && planePara[0] > 0)
				pointForFill.x -= planePara[3];
			pointForFill.y = proj_min.y + ((double)i / pointPitch);
			pointForFill.z = proj_min.z + ((double)j / pointPitch);
			pointForFill.rgba = color;
			cloud_filled_temp->points.push_back(pointForFill);
		}
	}

	rotation(0, 0) = cos(-angle); rotation(0, 1) = -sin(-angle);
	rotation(1, 0) = sin(-angle); rotation(1, 1) = cos(-angle);
	pcl::transformPointCloud(*cloud_filled_temp, *this->pointCloud, rotation);
	PointT leftDown, rightDown, rightUp, leftUp;
	pcl::getMinMax3D(*cloud_filled_temp, leftDown, rightUp); //元に戻す前のy-z平面に平行な面を用いる
	rightDown.x = rightUp.x;
	rightDown.y = rightUp.y;
	rightDown.z = leftDown.z;
	leftUp = leftDown;
	leftUp.z = rightUp.z;
	TransformPoint(&leftDown, (-1) * angle);
	TransformPoint(&rightDown, (-1) * angle);
	TransformPoint(&leftUp, (-1) * angle);
	TransformPoint(&rightUp, (-1) * angle);
	this->_leftDown = leftDown;
	this->_leftUp = leftUp;
	this->_rightDown = rightDown;
	this->_rightUp = rightUp;
}


/** @brief [1] to filter the plane point cloud at x, y, z dimension
 *         [2] update the boundary
 * @param axis
 */
void Plane::applyFilter(const string axis, float min, float max) {
	if (max < min) {
		PCL_ERROR("@setFilter <- min larger than max\n");
		return;
	}
	if (axis != "x" && axis != "y" && axis != "z") {
		PCL_ERROR("@setFilter <- you can only input x,y,z\n");
		return;
	}
	// fixme: filter x and y axis
	if (axis == "x" || axis == "y") {
		PCL_ERROR("@setFilter <- have write the code for x or y filter\n");
		return;
	}
	pcl::PassThrough<PointT> filterHeight;
	filterHeight.setInputCloud(this->pointCloud);
	filterHeight.setFilterFieldName(axis);
	filterHeight.setFilterLimits(min, max);
	filterHeight.setFilterLimitsNegative(false);
	filterHeight.filter(*this->pointCloud);
	// update the boundary points data
	this->updateBoundary();
}


void Plane::extendPlane(PointT a1, PointT a2, PointT b1, PointT b2, float pointPitch) {
	int color = this->pointCloud->points[0].rgba;
	generatePlanePointCloud(a1, a2, b1, b2, pointPitch, color);
	updateBoundary();
}
void Plane::extendPlane(PointT a1, PointT a2, PointT b1, PointT b2, float pointPitch, PlaneColor colorType) {
	generatePlanePointCloud(a1, a2, b1, b2, pointPitch, colorType2int(colorType));
	updateBoundary();
}


void Plane::append(Plane const &plane) {
	for (int i = 0; i < plane.pointCloud->size(); ++i) {
		this->pointCloud->points.push_back(plane.pointCloud->points[i]);
	}
}

float Plane::getEdgeLength(edgeType type) {
	switch (type) {
	case EdgeLeft:
		return pcl::geometry::distance(this->_leftUp, this->_leftDown);
	case EdgeRight:
		return pcl::geometry::distance(this->_rightUp, this->_rightDown);
	case EdgeUp:
		return pcl::geometry::distance(this->_leftUp, this->_rightUp);
	case EdgeDown:
		return pcl::geometry::distance(this->_leftDown, this->_rightDown);
	}
}

Eigen::Vector3d Plane::getNormal() {
	Eigen::Vector3d v0(_leftUp.x, _leftUp.y, _leftUp.z);
	Eigen::Vector3d v1(_rightUp.x, _rightUp.y, _rightUp.z);
	Eigen::Vector3d v2(_leftDown.x, _leftDown.y, _leftDown.z);
	Eigen::Vector3d V1 = -v0 + v1;
	Eigen::Vector3d V2 = -v0 + v2;
	return V2.cross(V1);
}

void Plane::removePointWithin(float xMin, float xMax, float yMin, float yMax, float zMin, float zMax) {

	if (xMin > xMax || yMin > yMax || zMin > zMax) {
		PCL_ERROR("@removePointWithin, the value input max < min");
		return;
	}
	pcl::ConditionOr<PointT>::Ptr range1(new pcl::ConditionOr<PointT>);  // this is not filter, but remain which parts we want
	range1->addComparison(pcl::FieldComparison<PointT>::ConstPtr(new pcl::FieldComparison<PointT>("z", pcl::ComparisonOps::LT, zMin)));
	range1->addComparison(pcl::FieldComparison<PointT>::ConstPtr(new pcl::FieldComparison<PointT>("z", pcl::ComparisonOps::GT, zMax)));

	pcl::ConditionOr<PointT>::Ptr range2(new pcl::ConditionOr<PointT>);
	range2->addComparison(pcl::FieldComparison<PointT>::ConstPtr(new pcl::FieldComparison<PointT>("x", pcl::ComparisonOps::LT, xMin)));
	range2->addComparison(pcl::FieldComparison<PointT>::ConstPtr(new pcl::FieldComparison<PointT>("x", pcl::ComparisonOps::GT, xMax)));

	pcl::ConditionAnd<PointT>::Ptr range3(new pcl::ConditionAnd<PointT>);
	range3->addComparison(pcl::FieldComparison<PointT>::ConstPtr(new pcl::FieldComparison<PointT>("y", pcl::ComparisonOps::LT, yMin)));
	range3->addComparison(pcl::FieldComparison<PointT>::ConstPtr(new pcl::FieldComparison<PointT>("y", pcl::ComparisonOps::GT, yMax)));
	pcl::ConditionOr<PointT>::Ptr range(new pcl::ConditionOr<PointT>);
	range->addCondition(range1);
	range->addCondition(range2);
	range->addCondition(range3);
	pcl::ConditionalRemoval<PointT> condrem;
	condrem.setCondition(range);
	condrem.setInputCloud(this->pointCloud);
	condrem.setKeepOrganized(false);
	condrem.filter(*this->pointCloud);
}

/**\paragraph Private Methods
 */

 /**@bug remember to update the boundary infomation
  */
void Plane::generatePlanePointCloud(PointT a1, PointT a2, PointT b1, PointT b2, float pointPitch, int color) {
	int num1 = (pcl::geometry::distance(a1, a2) * pointPitch);
	int num2 = (pcl::geometry::distance(b1, b2) * pointPitch);
	int num = num1 < num2 ? num1 : num2;
	vector<Eigen::Vector3f> pt1;
	pt1.push_back(Eigen::Vector3f(a1.x, a1.y, a1.z));
	pt1.push_back(Eigen::Vector3f(a2.x, a2.y, a2.z));
	vector<Eigen::Vector3f> pt2;
	pt2.push_back(Eigen::Vector3f(b1.x, b1.y, b1.z));
	pt2.push_back(Eigen::Vector3f(b2.x, b2.y, b2.z));
	Eigen::Vector3f stepLine1 = (-pt1[0] + pt1[1]) / num;
	Eigen::Vector3f stepLine2 = (-pt2[0] + pt2[1]) / num;
	vector<PointT> line1Points;
	vector<PointT> line2Points;
	for (size_t i = 0; i < num; i++) {
		PointT p;
		p.x = pt1[0][0] + i * stepLine1[0];
		p.y = pt1[0][1] + i * stepLine1[1];
		p.z = pt1[0][2] + i * stepLine1[2];
		line1Points.push_back(p);
	}

	for (size_t i = 0; i < num; i++) {
		PointT p;
		p.x = pt2[0][0] + i * stepLine2[0];
		p.y = pt2[0][1] + i * stepLine2[1];
		p.z = pt2[0][2] + i * stepLine2[2];
		line2Points.push_back(p);
	}

	for (size_t i = 0; i < num; i++) {
		generateLinePointCloud(line1Points[i], line2Points[i], pointPitch, color);
	}
}

void Plane::generateLinePointCloud(PointT pt1, PointT pt2, int pointPitch, int color) {
	int numPoints = pcl::geometry::distance(pt1, pt2) * pointPitch;
	float ratioX = (pt1.x - pt2.x) / numPoints;
	float ratioY = (pt1.y - pt2.y) / numPoints;
	float ratioZ = (pt1.z - pt2.z) / numPoints;
	for (size_t i = 0; i < numPoints; i++) {
		PointT p;
		p.x = pt2.x + i * (ratioX);
		p.y = pt2.y + i * (ratioY);
		p.z = pt2.z + i * (ratioZ);
		p.rgba = color;
		this->pointCloud->points.push_back(p);
	}
}

void Plane::updateBoundary() {
	PointT min, max;
	pcl::getMinMax3D(*this->pointCloud, min, max);
	PointT leftUp, leftDown, rightUp, rightDown;
	leftUp.z = max.z;
	leftDown.z = min.z;
	rightUp.z = max.z;
	rightDown.z = min.z;
	if (this->_abcd[0] * this->_abcd[1] > 0) {
		leftUp.x = min.x; leftUp.y = max.y;
		leftDown.x = min.x; leftDown.y = max.y;
		rightUp.x = max.x; rightUp.y = min.y;
		rightDown.x = max.x; rightDown.y = min.y;
	}
	else {
		leftUp.x = min.x; leftUp.y = min.y;
		leftDown.x = min.x; leftDown.y = min.y;
		rightUp.x = max.x; rightUp.y = max.y;
		rightDown.x = max.x; rightDown.y = max.y;
	}
	this->_leftUp = leftUp;
	this->_leftDown = leftDown;
	this->_rightUp = rightUp;
	this->_rightDown = rightDown;
}
