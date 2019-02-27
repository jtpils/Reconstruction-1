#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H

#include <iostream>
#include <vector>
#include "Plane.h"
typedef pcl::PointXYZRGB PointRGB;
typedef pcl::PointXYZRGBNormal PointT;
typedef pcl::PointCloud<PointT> PointCloudT;



class Reconstruction
{

public:
	Reconstruction(const string filePath);
    Reconstruction(PointCloudT::Ptr cloud);
	bool isPrintDebugInfo = true;
	bool isOutputEachStep = true;
	string outputPath = "OutputData/";
	//reconstructParas paras;
	PointCloudT::Ptr pointCloud;
	void downSampling(float leafSize);
	void applyRegionGrow(int NumberOfNeighbours, int SmoothnessThreshold, int CurvatureThreshold, int MinSizeOfCluster, int KSearch);
	void getClusterPts(PointCloudT::Ptr output);

	void applyRANSACtoClusters(float RANSAC_DistThreshold, float RANSAC_PlaneVectorThreshold,
		float RANSAC_MinInliers);
	// output methods
	void outputFile(const string path);
	void getPlane(PlaneOrientation ori, vector<Plane>& planes);
	vector<PointCloudT::Ptr> clusters; // store the result of region grow
	vector<Plane> ransacPlanes;
	
private:
	void debugPrint(stringstream& ss);
	void calculateRANSAC_plane(PointCloudT::Ptr cloud_cluster, pcl::PointIndices::Ptr sacInliers,
		pcl::ModelCoefficients::Ptr sacCoefficients, double distanceFromRANSACPlane);
	void calculateNormals(pcl::PointCloud <pcl::Normal>::Ptr &normals_all, int KSearch);
};



#endif
