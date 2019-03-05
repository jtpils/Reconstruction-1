#include <iostream>
#include <pcl/visualization/pcl_visualizer.h>
#include <SimpleView.h>
#include <Plane.h>
#include <Reconstruction.h>
#include <extract_walls.h>
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


void simpleView(const string& title, vector<PointCloudT::Ptr> &clusters) {
	vector<Eigen::Vector3i> colors;
	for (int k = 0; k <= clusters.size()/6; ++k) {
		colors.push_back(Eigen::Vector3i(255,0,0)); // red
		colors.push_back(Eigen::Vector3i(255,255,0)); // yellow
		colors.push_back(Eigen::Vector3i(0,0,255)); // blue
		colors.push_back(Eigen::Vector3i(0,255,0)); // green
		colors.push_back(Eigen::Vector3i(0,255,255)); // cyan
		colors.push_back(Eigen::Vector3i(255,0,255)); // pink
	}

	PointCloudT::Ptr coloredClusterPts(new PointCloudT);
	for (int l = 0; l < clusters.size(); ++l) {
		int color = 255 <<24 | colors[l][0] << 16 | colors[l][1] << 8 | colors[l][2];
		for (auto &p:clusters[l]->points) {
			p.rgba = color;
			coloredClusterPts->push_back(p);
		}
	}

	simpleView(title, coloredClusterPts);
}

void simpleView(const string& title, vector<EdgeLine> &lines){
	vector<Eigen::Vector3i> colors;
	for (int k = 0; k <= lines.size()/6; ++k) {
		colors.push_back(Eigen::Vector3i(255,0,0)); // red
		colors.push_back(Eigen::Vector3i(255,255,0)); // yellow
		colors.push_back(Eigen::Vector3i(0,0,255)); // blue
		colors.push_back(Eigen::Vector3i(0,255,0)); // green
		colors.push_back(Eigen::Vector3i(0,255,255)); // cyan
		colors.push_back(Eigen::Vector3i(255,0,255)); // pink
	}

	PointCloudT::Ptr linePts(new PointCloudT);
	for (int i = 0; i < lines.size(); ++i) {
		int color = 255 <<24 | colors[i][0] << 16 | colors[i][1] << 8 | colors[i][2];
		generateLinePointCloud(lines[i].p, lines[i].q, 20,color, linePts);
	}
	simpleView(title , linePts);
}
