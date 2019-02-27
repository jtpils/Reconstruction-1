#include <iostream>
#include <pcl/visualization/pcl_visualizer.h>
#include <SimpleView.h>
#include <Plane.h>
#include <Reconstruction.h>
using namespace std;
typedef pcl::PointXYZRGBNormal PointT;
typedef pcl::PointCloud<PointT> PointCloudT;
typedef pcl::PointCloud<pcl::PointXYZRGB> PointCloudRGB;

void copyOnlyRgba(PointCloudT::Ptr input, PointCloudRGB::Ptr output) {
	for (auto &i : input->points) {
		pcl::PointXYZRGB p;
		p.x = i.x;
		p.y = i.y;
		p.z = i.z;
		p.rgb = i.rgb;
		output->points.push_back(p);
	}
}

void simpleView(const string &title, const pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr &cloud) {
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer(title));
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_noNormal(new pcl::PointCloud<pcl::PointXYZRGB>);
	copyOnlyRgba(cloud, cloud_noNormal);
	viewer->setBackgroundColor(0, 0, 0);
	viewer->addPointCloud<pcl::PointXYZRGB>(cloud_noNormal, "1", 0);
	viewer->addCoordinateSystem(1.0);
	viewer->initCameraParameters();
	while (!viewer->wasStopped())
	{
		viewer->spinOnce(100);
		boost::this_thread::sleep(boost::posix_time::microseconds(100000));
	}
}

void simpleView(const string &title, const pcl::PointCloud<pcl::PointXYZ>::Ptr &cloud) {
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer(title));
	viewer->setBackgroundColor(0, 0, 0);
	viewer->addPointCloud<pcl::PointXYZ>(cloud, "1", 0);
	viewer->addCoordinateSystem(1.0);
	viewer->initCameraParameters();
	while (!viewer->wasStopped())
	{
		viewer->spinOnce(100);
		boost::this_thread::sleep(boost::posix_time::microseconds(100000));
	}
}

void simpleView(const string& title, vector<Plane> &planes) {

	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer(title));
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_noNormal(new pcl::PointCloud<pcl::PointXYZRGB>);
	for (auto plane : planes) {
		PointCloudT::Ptr cloud = plane.pointCloud;
		copyOnlyRgba(cloud, cloud_noNormal);
	}
	viewer->setBackgroundColor(0, 0, 0);
	viewer->addPointCloud<pcl::PointXYZRGB>(cloud_noNormal, "1", 0);
	viewer->addCoordinateSystem(1.0);
	viewer->initCameraParameters();
	while (!viewer->wasStopped())
	{
		viewer->spinOnce(100);
		boost::this_thread::sleep(boost::posix_time::microseconds(100000));
	}
}

void simpleView(const string& title, Reconstruction &re) {
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer(title));
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_noNormal(new pcl::PointCloud<pcl::PointXYZRGB>);
	copyOnlyRgba(re.pointCloud, cloud_noNormal);
	viewer->setBackgroundColor(0, 0, 0);
	viewer->addPointCloud<pcl::PointXYZRGB>(cloud_noNormal, "1", 0);
	viewer->addCoordinateSystem(1.0);
	viewer->initCameraParameters();
	while (!viewer->wasStopped())
	{
		viewer->spinOnce(100);
		boost::this_thread::sleep(boost::posix_time::microseconds(100000));
	}
}