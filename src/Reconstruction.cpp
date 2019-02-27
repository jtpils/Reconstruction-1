#include <iostream>
#include <vector>
#include <stdexcept>
#include <math.h>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
#include <pcl/io/obj_io.h>
#include <pcl/search/search.h>
#include <pcl/search/kdtree.h>
#include <pcl/features/normal_3d.h>
#include <pcl/segmentation/region_growing.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/common/pca.h>
#include <pcl/common/common.h>
#include <pcl/common/geometry.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/conditional_removal.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/filters/project_inliers.h>
#include <pcl/filters/passthrough.h>
#include "Reconstruction.h"
#include "Plane.h"
using namespace std;


void Reconstruction::debugPrint(stringstream& ss) {
	if (!this->isPrintDebugInfo) return;
	cout << ss.str() << endl;
}


std::vector<std::string> split(const std::string& s, char delimiter)
{
	std::vector<std::string> tokens;
	std::string token;
	std::istringstream tokenStream(s);
	while (std::getline(tokenStream, token, delimiter))
	{
		tokens.push_back(token);
	}
	return tokens;
}


Reconstruction::Reconstruction(const string filePath) {
	stringstream ss;
	ss << "Input File: " << filePath;
	debugPrint(ss);
	PointCloudT::Ptr tmp(new PointCloudT);
	this->pointCloud = tmp;

	string fileType = filePath.substr(filePath.length() - 3);
	if (!(fileType == "ply" || fileType == "obj" || fileType == "txt" || fileType == "pcd"))
		throw invalid_argument("the file type is not allowed");
	if (fileType == "ply") {
		if (pcl::io::loadPLYFile <PointT>(filePath, *this->pointCloud) == -1) { // the file doesnt exist
			throw invalid_argument("Cannot load the input file, please check and try again!\n");
		}
	}
	else if (fileType == "obj") {
		if (pcl::io::loadOBJFile <PointT>(filePath, *this->pointCloud) == -1) { // the file doesnt exist
			throw invalid_argument("Cannot load the input file, please check and try again!\n");
		}
	}
	else if (fileType == "pcd") {
		if (pcl::io::loadPCDFile <PointT>(filePath, *this->pointCloud) == -1) { // the file doesnt exist
			throw invalid_argument("Cannot load the input file, please check and try again!\n");
		}
	}
	else if (fileType == "txt") {
		std::ifstream file(filePath);
		std::string str;
		while (std::getline(file, str))
		{
			vector<string> tokens;
			tokens = split(str, ',');
			if (tokens.size() != 7) {
				PCL_WARN("@s IMPORT FAILIURE\n", str);
				continue;
			}
			PointT p;
			p.x = stof(tokens[0]); p.y = stof(tokens[1]); p.z = stof(tokens[2]);
			p.r = stof(tokens[4]); p.g = stof(tokens[5]); p.b = stof(tokens[6]);
			p.a = stof(tokens[3]);
			this->pointCloud->push_back(p);
		}
	}
}

Reconstruction::Reconstruction(PointCloudT::Ptr cloud) {
    this->pointCloud = cloud;
}

void Reconstruction::downSampling(float leafSize)
{
	if (leafSize == -1) {
		throw invalid_argument("please set leafSize parameters");
	}
	stringstream ss;
	ss << "\nDownSampling...: leafSize-> " << leafSize << "\n";

	ss << "Before-> " << this->pointCloud->points.size();

	pcl::VoxelGrid<PointT> sor;
	sor.setInputCloud(this->pointCloud);

	sor.setLeafSize(leafSize, leafSize, leafSize); //5mm“_ŒQ‚É‘Î‚µ‚ÄA0.01‚¾‚ÆA’l‚ª¬‚³‚·‚¬‚é‚Æ‚¢‚¤ƒGƒ‰[‚ªo‚é
	sor.filter(*this->pointCloud);
	ss << "  After-> " << this->pointCloud->points.size();
	debugPrint(ss);
}

void Reconstruction::applyRegionGrow(int NumberOfNeighbours, int SmoothnessThreshold, int CurvatureThreshold, int MinSizeOfCluster, int KSearch)
{
	stringstream ss;
	ss << "\nRegionGrowing..." << "\n" << "MaxClusterSize: 100000 - NumberOfNeighbours: " << NumberOfNeighbours << "\n";
	ss << "SmoothnessThreshold: " << SmoothnessThreshold << "\n" << "CurvatureThreshold: " << CurvatureThreshold << "\n";
	ss << "Min size of Cluster: " << MinSizeOfCluster << "\n";
	std::vector <pcl::PointIndices> clustersIndices;
	pcl::search::Search<PointT>::Ptr tree = boost::shared_ptr<pcl::search::Search<PointT> >(new pcl::search::KdTree<PointT>);
	pcl::PointCloud <pcl::Normal>::Ptr normals_all(new pcl::PointCloud <pcl::Normal>);
	calculateNormals(normals_all, KSearch);
	pcl::RegionGrowing<PointT, pcl::Normal> reg;
	reg.setMinClusterSize(0);
	reg.setMaxClusterSize(100000);
	reg.setSearchMethod(tree);
	reg.setNumberOfNeighbours(NumberOfNeighbours);
	reg.setInputCloud(this->pointCloud);
	reg.setInputNormals(normals_all);
	reg.setSmoothnessThreshold(static_cast<float>(SmoothnessThreshold / 180.0 * M_PI));
	reg.setCurvatureThreshold(CurvatureThreshold);
	reg.extract(clustersIndices);
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr colored_cloud = reg.getColoredCloud();
	for (size_t i = 0; i < clustersIndices.size(); ++i) {
		if (clustersIndices[i].indices.size() < MinSizeOfCluster) continue;
		PointCloudT::Ptr singleCluster(new PointCloudT);
		for (auto &p : clustersIndices[i].indices) {
			singleCluster->points.push_back(this->pointCloud->points[p]);
		}
		this->clusters.push_back(singleCluster);
	}
	ss << "num of Clusters: " << this->clusters.size();
	debugPrint(ss);
}
void Reconstruction::getClusterPts(PointCloudT::Ptr output){
	output->resize(0);

	for(auto &cluster:clusters){
		for (auto &p: cluster->points) {
			output->push_back(p);
		}
	}

}
void Reconstruction::applyRANSACtoClusters(float RANSAC_DistThreshold, float RANSAC_PlaneVectorThreshold,
	float RANSAC_MinInliers) {
	stringstream ss;
	ss << "\nApplying RANSAC to all Clusters...";
	ss << "\nRANSAC DistThreshold: " << RANSAC_DistThreshold;
	ss << "\nRANSAC PlaneVectorThreshold: " << RANSAC_PlaneVectorThreshold;
	ss << "\nRANSAC RANSAC MinInliers: " << RANSAC_MinInliers*100 << "%";
	
	if (this->clusters.size() == 0) {
		throw invalid_argument("Cluster Size == 0!\n");
	}
	for (auto &cluster : this->clusters)
	{
		// mark:  apply ransac
		pcl::ModelCoefficients::Ptr sacCoefficients(new pcl::ModelCoefficients);
		pcl::PointIndices::Ptr sacInliers(new pcl::PointIndices);
		calculateRANSAC_plane(cluster, sacInliers, sacCoefficients, RANSAC_DistThreshold);
		double a1, b1, c1, d1;
		a1 = sacCoefficients->values[0];
		b1 = sacCoefficients->values[1];
		c1 = sacCoefficients->values[2];
		d1 = sacCoefficients->values[3];
		Eigen::Vector4d abcd(a1, b1, c1, d1);
		PointCloudT::Ptr extracted_cloud(new PointCloudT);
		pcl::ExtractIndices<PointT> extract;
		extract.setInputCloud(cluster);
		extract.setIndices(sacInliers);
		extract.setNegative(false);
		extract.filter(*extracted_cloud);
		// control the num of RANSAC plane
		float threshold = RANSAC_PlaneVectorThreshold;
		if (sacInliers->indices.size() / cluster->size() >= RANSAC_MinInliers) {
			Plane plane(extracted_cloud, abcd);
			if (threshold > abs(a1) && threshold > abs(b1)) {
				plane.orientation = PlaneOrientation::Horizontal;
			}
			else {
				plane.orientation = PlaneOrientation::Vertical;
			}
			this->ransacPlanes.push_back(plane);
		}
	}

	ss << "\nInput nums of clusters: " << this->clusters.size();
	ss << "\nOutput nums of ransac clusters: " << this->ransacPlanes.size();
	debugPrint(ss);
}


void Reconstruction::outputFile(const string path)
{
	pcl::io::savePLYFile(path, *this->pointCloud);
}

void Reconstruction::getPlane(PlaneOrientation ori, vector<Plane>& planes) {
	if (this->ransacPlanes.size() < 1) {
		throw invalid_argument("The ransac Planes equals to 0");
	}
	planes.resize(0);
	for (Plane& plane : this->ransacPlanes) {
		if (plane.orientation == ori)
		{
			planes.push_back(plane);
		}
	}
}


// private methods

void Reconstruction::calculateNormals(pcl::PointCloud <pcl::Normal>::Ptr &normals_all, int KSearch)
{

	//1-1. generating the normal for each point
	if (this->pointCloud->points[0].normal_x == 0 &&
		this->pointCloud->points[0].normal_y == 0 &&
		this->pointCloud->points[0].normal_z == 0) {
		stringstream ss;
		ss << "The point you input doesn't contain normals, calculating normals...";
		debugPrint(ss);
		pcl::search::Search<PointT>::Ptr tree = boost::shared_ptr<pcl::search::Search<PointT> >(new pcl::search::KdTree<PointT>);
		pcl::NormalEstimation<PointT, pcl::Normal> normal_estimator;
		normal_estimator.setSearchMethod(tree);
		normal_estimator.setInputCloud(this->pointCloud);
		normal_estimator.setKSearch(KSearch);
		normal_estimator.compute(*normals_all);

		for (size_t i = 0; i < normals_all->points.size(); ++i)
		{
			this->pointCloud->points[i].normal_x = normals_all->points[i].normal_x;
			this->pointCloud->points[i].normal_y = normals_all->points[i].normal_y;
			this->pointCloud->points[i].normal_z = normals_all->points[i].normal_z;
		}
	}
	else
	{
		for (size_t i = 0; i < this->pointCloud->points.size(); ++i)
		{
			pcl::Normal normal_temp;
			normal_temp.normal_x = this->pointCloud->points[i].normal_x;
			normal_temp.normal_y = this->pointCloud->points[i].normal_y;
			normal_temp.normal_z = this->pointCloud->points[i].normal_z;
			normals_all->push_back(normal_temp);
		}
	}
}

void Reconstruction::calculateRANSAC_plane(PointCloudT::Ptr cloud_cluster, pcl::PointIndices::Ptr sacInliers,
	pcl::ModelCoefficients::Ptr sacCoefficients, double distanceFromRANSACPlane) {
	pcl::SACSegmentation<PointT> seg;
	seg.setOptimizeCoefficients(true);
	seg.setModelType(pcl::SACMODEL_PLANE);
	seg.setMethodType(pcl::SAC_RANSAC);
	seg.setDistanceThreshold(distanceFromRANSACPlane);
	seg.setInputCloud(cloud_cluster);
	seg.segment(*sacInliers, *sacCoefficients);
}