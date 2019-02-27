#ifndef PLANE_H
#define PLANE_H

#include <iostream>
#include <vector>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
#include <pcl/search/search.h>
#include <pcl/search/kdtree.h>
#include <pcl/features/normal_3d.h>
#include <pcl/segmentation/region_growing.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/segmentation/sac_segmentation.h> //RANSACのため
#include <math.h>
#include <pcl/filters/project_inliers.h> //平面に投影するため
#include <pcl/filters/passthrough.h>
#include <pcl/common/pca.h>
#include <pcl/common/common.h>
#include <pcl/io/obj_io.h> //obj形式で保存するため
#include <pcl/filters/voxel_grid.h> //ダウンサンプリングのため
#include <pcl/common/geometry.h>
#include <pcl/filters/conditional_removal.h>
using namespace std;
typedef pcl::PointXYZRGB PointRGB;
typedef pcl::PointXYZRGBNormal PointT;
typedef pcl::PointCloud<PointT> PointCloudT;
enum edgeType {
	EdgeLeft,
	EdgeRight,
	EdgeDown,
	EdgeUp,
};

enum PlaneOrientation {
	Vertical,
	Horizontal,
	UNDEFINED
};

enum PlaneType {
	PlaneType_MainWall,
	PlaneType_Ceiling,
	PlaneType_Ground,
	PlaneType_Other,
};

enum PlaneColor {
	Color_Red,
	Color_Yellow,
	Color_Blue,
	Color_Green,
	Color_White,
	Color_Peach,
	Color_Random,
};


class Plane
{
private:

	PointT _leftDown, _rightDown, _rightUp, _leftUp;
	Eigen::Vector4d _abcd;
	PlaneType _type = PlaneType::PlaneType_Other;

	void generatePlanePointCloud(PointT a1, PointT a2, PointT b1, PointT b2, float pointPitch, int color);
	void generateLinePointCloud(PointT pt1, PointT pt2, int pointPitch, int color);
	void updateBoundary();

public:
	Plane();
	Plane(PointCloudT::Ptr rawPointCloud, Eigen::Vector4d abcd);
	Plane(PointT a1, PointT a2, PointT b1, PointT b2, float pointPitch, PlaneColor color);

	PointCloudT::Ptr pointCloud;
	const PointT& leftDown()         const { return _leftDown; }
	const PointT& rightDown()        const { return _rightDown; }
	const PointT& rightUp()          const { return _rightUp; }
	const PointT& leftUp()           const { return _leftUp; }
	const Eigen::Vector4d& abcd()    const { return _abcd; }
	const PlaneType& type()          const { return _type; }
	Plane *coveredPlane = nullptr;
	PlaneOrientation orientation = PlaneOrientation::UNDEFINED;
	void setType(PlaneType type);
	void setColor(PlaneColor colorType);
	// This is for fill the vertical plane based on their boundary.
	// pros: will fill all the plane, which means ignore same vacant parts of original pts.
	// cons: the high will change to the projection length to y-z plane
	void filledPlane(int pointPitch);
	void applyFilter(const string axis, float min, float max);
	void extendPlane(PointT a1, PointT a2, PointT b1, PointT b2, float pointPitch, PlaneColor colorType);
	void extendPlane(PointT a1, PointT a2, PointT b1, PointT b2, float pointPitch);
	void append(Plane const &plane);
	float getEdgeLength(edgeType type);
	Eigen::Vector3d getNormal();
	void removePointWithin(float xMin, float xMax, float yMin, float yMax, float zMin, float zMax);

};

#endif
