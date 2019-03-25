//
// Created by czh on 10/17/18.
//
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <chrono>
#include <vector>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
#include <pcl/io/vtk_io.h>
#include <pcl/search/search.h>
#include <pcl/search/kdtree.h>
#include <pcl/features/normal_3d.h>
#include <pcl/segmentation/region_growing.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/segmentation/sac_segmentation.h> //RANSACのため
#include <cmath>
#include <pcl/filters/project_inliers.h> //平面に投影するため
#include <pcl/filters/passthrough.h>
#include <pcl/common/pca.h>
#include <pcl/common/common.h>
#include <pcl/io/obj_io.h> //obj形式で保存するため
#include <pcl/filters/voxel_grid.h> //ダウンサンプリングのため
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/common/geometry.h>
#include <pcl/filters/conditional_removal.h>
#include <jsoncpp/json/value.h>
#include <jsoncpp/json/json.h>
#include "main.h"
#include "Plane.h"
#include "SimpleView.h"
#include "Reconstruction.h"
#include "DxfExporter.h"
#include <pcl/surface/concave_hull.h>
#include <pcl/surface/convex_hull.h>
#include <yaml-cpp/yaml.h>
#include "DebugFileExporter.h"
#include <pcl/console/print.h>

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arr_simple_point_location.h>
using namespace std;
using namespace pcl;
typedef pcl::PointXYZRGB PointRGB;
typedef pcl::PointXYZRGBNormal PointT;
typedef pcl::PointCloud<PointT> PointCloudT;

YAML::Node config ;

// parse the parameters set in config.yaml file
template<typename T>
T Paras(string a)
{
    return config[a].as<T>();
}

template<typename T>
T Paras(string a, string b)
{
    YAML::Node aa = config[a];
    return aa[b].as<T>();
}

int main(int argc, char** argv) {

	PCL_WARN("This program is based on assumption that ceiling and ground on the X-Y  \n");
	pcl::console::setVerbosityLevel(pcl::console::L_ALWAYS); // for not show the warning
	#ifdef _WIN32
		//string fileName = argv[2];
		string fileName = "TestData/Room_A.ply";
	#elif defined __unix__
		string fileName = "";
		string configPath = "../config.yaml";
		if (argc == 3) {
		    fileName = argv[1];
            configPath = argv[2];
		}else{
			// pre-set 3 files for debug
            cout << "1. RoomA \n";
		    cout << "2. RoomE\n";
            cout << "3. Hall\n";
		    int index = 0;
		    cout << "Please choose test dataset " << " \n";
            cin >> index;
            if (index == 1) fileName = "/home/czh/Desktop/PointCloudDataset/LaserScanner/old/Room_A.ply";
            else if (index == 2) fileName = "/home/czh/Desktop/PointCloudDataset/LaserScanner/roomE_empty/roomE_empty.ply";
            else if (index == 3) fileName = "/home/czh/Desktop/PointCloudDataset/LaserScanner/cooridor/Downsampled_coorridor_all.ply";
        }
    #endif
    config = YAML::LoadFile(configPath);
	Reconstruction re(fileName);
	if (Paras<bool>("View","original")) simpleView("input original point clouds", re.pointCloud);

	smoothNoise(re.pointCloud, Paras<int>("SmoothNoise","K"), Paras<int>("SmoothNoise","alpha"));
	if (Paras<bool>("View","smoothNoise")) simpleView("[smooth Noise]", re.pointCloud);

    densityFilter(re.pointCloud);
	if (Paras<bool>("View","densityFilter")) simpleView("[densityFilter]", re.pointCloud);

	PointCloudT::Ptr twoDimPts(new PointCloudT);

    float heightLow,heightHigh;
    detectHeightRange(re.pointCloud,heightHigh, heightLow);

	convert2D(re.pointCloud,twoDimPts);
	if (Paras<bool>("View","2D")) simpleView("[convert2D]", twoDimPts);

    PointCloudT::Ptr largestComp(new PointCloudT);
    findBiggestComponent2D(twoDimPts, largestComp);
	if (Paras<bool>("View","BiggestComponent")) simpleView("[findBiggestComponent2D]", largestComp);

    PointCloudT::Ptr hullOutput(new PointCloudT);
    computeHull(largestComp, hullOutput, Paras<float>("ComputeHull","alpha"));

	if (Paras<bool>("View","hull")) simpleView("[compute hull] ", hullOutput);

    // extract lines from edge point clouds
    vector<EdgeLine> edgeLines;
    extractLineFromEdge(hullOutput, edgeLines);
	if (Paras<bool>("View","extractLine")) simpleView("[extractLineFromEdge] result" , edgeLines);

	// extract lines from edge point clouds
	vector<EdgeLine> edgeLines2;
	extractLineFromEdge2(hullOutput, edgeLines2);
	if (Paras<bool>("View","extractLine")) simpleView("[extractLineFromEdge2] result " , edgeLines2);

	//minimumCut(edgeLines, largestComp);

	vector< pair<vector<PointT>, float> > beamPlanes;
	vector< pair<EdgeLine, float>> beamLines;
	if (Paras<bool>("BeamRecons","isExist")) BeamRANSAC(re.pointCloud,heightHigh, beamPlanes, beamLines);

	findLinkedLines(edgeLines2);
	cout << "after find linked lines, extract " << edgeLines.size() << " lines\n";
	if (Paras<bool>("View","linkedEdges")) simpleView("line pts after findLinkedLines" , edgeLines2);
	exportToDxf("./output.dxf", edgeLines2, beamPlanes, beamLines, heightLow, heightHigh);
    return 0;
}


void generateLinePointCloud(PointT pt1, PointT pt2, int pointPitch, int color, PointCloudT::Ptr output) {
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
		output->points.push_back(p);
	}
}

void generateLinePointCloud(PointT a, PointT b, PointT c, PointT d, int pointPitch, int color, PointCloudT::Ptr output){
	generateLinePointCloud(a,b,pointPitch, color, output);
	generateLinePointCloud(b,c,pointPitch, color, output);
	generateLinePointCloud(c,d,pointPitch, color, output);
	generateLinePointCloud(d,a,pointPitch, color, output);
}

bool onSegment(PointT p, PointT q, PointT r)
{
	// Given three colinear points p, q, r, the function checks if
	// point q lies on line segment 'pr'
	if (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) &&
		q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y) &&
		q.z <= max(p.z, r.z) && q.z >= min(p.z, r.z))
		return true;
	return false;
}

float orientation(PointT p, PointT q, PointT r) {
	// To find orientation of ordered triplet (p, q, r).
	// The function returns following values
	// 0 --> p, q and r are colinear
	// 1 --> Clockwise
	// 2 --> Counterclockwise
	float val = (q.y - p.y) * (r.x - q.x) -
		(q.x - p.x) * (r.y - q.y);

	if (val == 0) return 0;  // colinear

	return (val > 0) ? 1 : 2; // clock or counterclock wise
}

bool isIntersect(PointT p1, PointT q1, PointT p2, PointT q2)
{
	// Find the four orientations needed for general and
	// special cases
	float o1 = orientation(p1, q1, p2);
	float o2 = orientation(p1, q1, q2);
	float o3 = orientation(p2, q2, p1);
	float o4 = orientation(p2, q2, q1);

	// General case
	if (o1 != o2 && o3 != o4)
		return true;

	// Special Cases
	// p1, q1 and p2 are colinear and p2 lies on segment p1q1
	if (o1 == 0 && onSegment(p1, p2, q1)) return true;

	// p1, q1 and q2 are colinear and q2 lies on segment p1q1
	if (o2 == 0 && onSegment(p1, q2, q1)) return true;

	// p2, q2 and p1 are colinear and p1 lies on segment p2q2
	if (o3 == 0 && onSegment(p2, p1, q2)) return true;

	// p2, q2 and q1 are colinear and q1 lies on segment p2q2
	if (o4 == 0 && onSegment(p2, q1, q2)) return true;

	return false; // Doesn't fall in any of the above cases
}

void smoothNoise(PointCloudT::Ptr input, int K, float alpha) {
    int n = input->size();

    pcl::KdTreeFLANN<PointT> kdtree;
    kdtree.setInputCloud (input);
    unordered_set<int> outliers;
    for (int i = 0; i < n; ++i) {
        //if (outliers.count(i)) continue;  // if points are noise, we ignore
        std::vector<int> pointIdxNKNSearch(K);
        std::vector<float> pointNKNSquaredDistance(K);
        if ( kdtree.nearestKSearch (input->points[i], K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0 )
        {
        	for (auto&value : pointNKNSquaredDistance) value = sqrt(value);
            float mean, stdev;
            calculateMeanStandardDev(pointNKNSquaredDistance, mean, stdev);
            for (int i = 0; i < pointNKNSquaredDistance.size(); ++i) {
            	if (pointNKNSquaredDistance[i] > (mean + alpha*stdev) ||
            	    pointNKNSquaredDistance[i] < (mean - alpha*stdev) ) {
            		// if out out range, put them into outlier
					outliers.insert(pointIdxNKNSearch[i]);
				}
            }

        }
    }
	pcl::PointIndices::Ptr indices(new pcl::PointIndices());
	indices->indices.insert(indices->indices.end(), outliers.begin(), outliers.end());

	PointCloudT::Ptr filtered(new PointCloudT);
	pcl::ExtractIndices<PointT> extract;
	extract.setInputCloud(input);
	extract.setIndices(indices);
	extract.setNegative(true);
	extract.filter(*input);
	cout << "[smooth noise] remove " << indices->indices.size() << " points from " <<n << ". remind " << input->size() << endl;

}

void densityFilter(PointCloudT::Ptr input){
	float alpha2   = Paras<float>("DensityFilter","alpha2");
    float leafSize = Paras<float>("DensityFilter","leafSize");
	int n = input->size();
	PointCloudT::Ptr filtered(new PointCloudT);

	// voxel grid filterd is used
	pcl::VoxelGrid<PointT> sor;
	sor.setInputCloud (input);
	sor.setLeafSize (leafSize, leafSize, leafSize);
	sor.filter(*filtered);

	// get the divisions in i,j,k three dimensions.
	Eigen::Vector3i divisions = sor.getNrDivisions();


	vector<vector<int>> voxels(divisions[0]*divisions[1]*divisions[2], vector<int>());

	// group points based on grid id
	for (int i = 0; i < input->points.size(); ++i) {
		auto cord = sor.getGridCoordinates(input->points[i].x, input->points[i].y, input->points[i].z);
		cord = cord - sor.getMinBoxCoordinates();
		voxels[cord[0]*divisions[1]*divisions[2] + cord[1]*divisions[2] + cord[2]].push_back(i);
	}

	// calculate num of points in each voxel
	vector<float> voxelsNums;
	for (auto &voxel: voxels) {
		if (voxel.size()!=0) voxelsNums.push_back(voxel.size());
	}
	float mean, stdev;
	calculateMeanStandardDev(voxelsNums, mean, stdev);
	cout <<  "[densityFilter] mean " << mean << ", stdev " << stdev << endl;

	// discard those points in the voxel that density is lower than the threshold
	PointCloudT::Ptr output(new PointCloudT);
	for (auto &voxel: voxels) {
		if (voxel.size() < (mean - stdev*alpha2)) continue;
		for (auto &idx :  voxel) output->push_back(input->points[idx]);
	}
	copyPointCloud(*output, *input);
	cout << "[density Filter] left " << input->size() << " from " << n << endl;
}

void detectHeightRange(PointCloudT::Ptr input, float& high, float& low){
    unordered_map<float, int> hist;

    // count the points in z value
    for(auto&p : input->points) {
    	// only want to remain first place decimal
		hist[roundf(p.z * 10) / 10]++;
	}

    // we need to find two z values that have most points
    int count1 = 0, count2 = 0; // count1 >= count2
    int maxNum = 0;
    for(auto h:hist) {
    	maxNum = max(h.second, maxNum);
		if (h.second > count1) {
			count2 = count1;
			low = high;
			count1 = h.second;
			high = h.first;
		} else if (h.second > count2) {
			count2 = h.second;
			low = h.first;
		}
	}
    // height difference is expected to be larger than 2 meters
    assert(abs(high-low) >= 2);
    if (high < low) swap(high, low);
    cout << "[detectHeightRange] [" << low << " " << high << "] max num "<< maxNum << "\n";
}

void extractTopPts(PointCloudT::Ptr input, PointCloudT::Ptr output, float highest, float dimension, reconstructParas paras){
	//PointCloudT::Ptr topTemp(new PointCloudT);

	pcl::PassThrough<PointT> filterZ;
	filterZ.setInputCloud(input);
	filterZ.setFilterFieldName("z");

	filterZ.setFilterLimits(highest - 2 * dimension, highest);
	filterZ.filter(*output);

	Reconstruction topTmp_re(output);
	topTmp_re.applyRegionGrow(paras.NumberOfNeighbours, paras.SmoothnessThreshold,
							  paras.CurvatureThreshold, paras.MinSizeOfCluster, paras.KSearch);
	topTmp_re.getClusterPts(output);
}

void convert2D(PointCloudT::Ptr input,PointCloudT::Ptr output) {
    pcl::copyPointCloud(*input, *output);
    for(auto& p : output->points)  p.z = 0;
	cout << "[convert2D] points num " << output->size() << endl;
}

void computeHull(PointCloudT::Ptr input, PointCloudT::Ptr output, float alpha){
	// since concavehull only support pointXYZ, we need to convert input into pointXYZ type
	pcl::ConcaveHull<pcl::PointXYZ> concave_hull;
	pcl::PointCloud<pcl::PointXYZ>::Ptr tmpInput(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr tmpOutput(new pcl::PointCloud<pcl::PointXYZ>);

	for (auto& p: input->points) {
		pcl::PointXYZ q;
		q.x = p.x;
		q.y = p.y;
		q.z = p.z;
		tmpInput->push_back(q);
	}
	concave_hull.setInputCloud (tmpInput);
	// if alpha is smaller, the output is with more details
	concave_hull.setAlpha (alpha);
	concave_hull.reconstruct (*tmpOutput);
	for (auto& p: tmpOutput->points) {
		pcl::PointXYZRGBNormal q;
		q.x = p.x;
		q.y = p.y;
		q.z = p.z;
		q.rgba = INT32_MAX;
		output->push_back(q);
	}

	cout << "[computeHull]: before " << input->size() << " -> after " << output->size() << "\n";
}

void findBiggestComponent2D(PointCloudT::Ptr input, PointCloudT::Ptr output){
	// the size of input point cloud should be larger than 0
    assert(input->points[0].z == 0);
    PointCloudT::Ptr tmp(new PointCloudT);
    pcl::copyPointCloud(*input, *tmp);
    vector<PointCloudT::Ptr> groups;
    // calculate the time cost
    auto start = std::chrono::system_clock::now();
    // separate point clouds into several groups
	separatePtsToGroups(tmp, 0.1, groups);
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    cout << "[findBiggestComponent2D] found " << groups.size() << " groups from " << tmp->size() << " pts. ";

    // find the one which contains most of points
    int maxNum = 0;
    int id = -1;
    for (int i = 0; i < groups.size(); ++i) {
        if (groups[i]->points.size() > maxNum) {
            maxNum = groups[i]->points.size();
            id = i;
        }
    }
    cout << "the max size of groups is  " << maxNum << "  ";
	std::cout << "Elapsed time: " << elapsed.count() << " s\n";
    assert(id >= 0);
    for (auto p : groups[id]->points) output->push_back(p);

}

void extractLineFromEdge(PointCloudT::Ptr input, vector<EdgeLine>& edgeLines){
	PointCloudT::Ptr inputTemp(new PointCloudT);
	copyPointCloud(*input, *inputTemp);
	vector<PointCloudT::Ptr> roofEdgeClusters; // record the different lines
	vector<Eigen::VectorXf>  roofEdgeClustersCoffs; // record the coefficients of lines
	assert(inputTemp->size() > 0);
	int iter = 0;
	while (iter++ < 1000 && inputTemp->size() > 1) {
		pcl::ModelCoefficients::Ptr sacCoefficients(new pcl::ModelCoefficients);
		pcl::PointIndices::Ptr sacInliers(new pcl::PointIndices);
		pcl::SACSegmentation<PointT> seg;
		seg.setOptimizeCoefficients(true);
		seg.setModelType(pcl::SACMODEL_LINE);
		seg.setMethodType(pcl::SAC_RANSAC);
		seg.setDistanceThreshold(0.05);
		seg.setInputCloud(inputTemp);
		seg.segment(*sacInliers, *sacCoefficients);
		Eigen::VectorXf coff(6);
		coff << sacCoefficients->values[0], sacCoefficients->values[1], sacCoefficients->values[2],
				sacCoefficients->values[3], sacCoefficients->values[4], sacCoefficients->values[5];

		PointCloudT::Ptr extracted_cloud(new PointCloudT);
		pcl::ExtractIndices<PointT> extract;
		extract.setInputCloud(inputTemp);
		extract.setIndices(sacInliers);
		extract.setNegative(false); // we only want filtered part
		extract.filter(*extracted_cloud);

		vector<PointCloudT::Ptr> tmpClusters;
		separatePtsToGroups(extracted_cloud, 0.5, tmpClusters); // use this method may make several parts as one part since RANSAC agrees
		for (auto &c:tmpClusters) {
			roofEdgeClusters.push_back(c);
			roofEdgeClustersCoffs.push_back(coff);
		}
		extract.setNegative(true); // remove the found parts. otherwise the program won't stop
		extract.filter(*inputTemp);
	}

	assert(roofEdgeClusters.size() > 0);
	//simpleView("[extractLineFromEdge] clusters",roofEdgeClusters);

	for (int i = 0; i < roofEdgeClusters.size(); ++i) {
		EdgeLine line;
		ptsToLine(roofEdgeClusters[i], roofEdgeClustersCoffs[i], line);
		if (pcl::geometry::distance(line.p,line.q) > 0.1) edgeLines.push_back(line);
	}
	cout << "[extractLineFromEdge] extract " << edgeLines.size() << " lines"<< endl;

}

void extractLineFromEdge2(PointCloudT::Ptr input, vector<EdgeLine>& edgeLines){
	vector<PointCloudT::Ptr> reginGrow2DOutput;
	vector<Eigen::VectorXf> roofEdgeClustersCoffs;

	// run region grow to extract clusters of boundary points
	regionGrow2D(input, reginGrow2DOutput);

	assert(reginGrow2DOutput.size() > 0);
	if (Paras<bool>("View","extractLine")) simpleView("[extractLineFromEdge2] regionGrow2DOutput",reginGrow2DOutput);

	// using new Generate Pts
	for (auto &Pts : reginGrow2DOutput) {
	    // for each cluster run RANSAC to extract Line model
		pcl::ModelCoefficients::Ptr sacCoefficients(new pcl::ModelCoefficients);
		pcl::PointIndices::Ptr sacInliers(new pcl::PointIndices);
		pcl::SACSegmentation<PointT> seg;
		seg.setOptimizeCoefficients(true);
		seg.setModelType(pcl::SACMODEL_LINE);
		seg.setMethodType(pcl::SAC_RANSAC);
		seg.setDistanceThreshold(0.05);
		seg.setInputCloud(Pts);
		seg.segment(*sacInliers, *sacCoefficients);
		Eigen::VectorXf coff(6);
		coff << sacCoefficients->values[0], sacCoefficients->values[1], sacCoefficients->values[2],
				sacCoefficients->values[3], sacCoefficients->values[4], sacCoefficients->values[5];
		roofEdgeClustersCoffs.push_back(coff);
	}

	for (int i = 0; i < reginGrow2DOutput.size(); ++i) {
		EdgeLine line;
		// find the start and end points of each line
		ptsToLine(reginGrow2DOutput[i], roofEdgeClustersCoffs[i], line);
		if (pcl::geometry::distance(line.p,line.q) > 0.1) edgeLines.push_back(line);
	}

	cout << "[extractLineFromEdge2] extract " << edgeLines.size() << " lines"<< endl;

}

void separatePtsToGroups(PointCloudT::Ptr input, float radius, vector<PointCloudT::Ptr>& output){

	// we use K-d tree to iteratively find the points around seed point in radius
	// and set surrounding points as seed points to iteratively group points
	//TODO this parts can only return indices instead of points, which saves more memory
	assert(input->size() > 0);
    pcl::KdTreeFLANN<PointT> kdtree;
    kdtree.setInputCloud(input);
    unordered_set<int> found;
    int k = 0;
	while(found.size() < input->size()) {
        vector<int> group;  // use to record found points
        stack<PointT> seeds;

        while(found.count(k)) k++; // if this point has been already been found, skip it
        if (k >= input->size()) break;
        seeds.push(input->points[k]);
        found.insert(k);
        group.push_back(k);

        while(!seeds.empty()) {
            PointT seed = seeds.top();
            seeds.pop();
            vector<int> pointIdx;
            vector<float> dist;
            kdtree.radiusSearch(seed, radius, pointIdx,dist);
			for(auto& ix : pointIdx) {
                if (!found.count(ix))  {
                    found.insert(ix);
                    seeds.push(input->points[ix]);
                    group.push_back(ix); // put them into the group
                }
            }
        }
        PointCloudT::Ptr foundPts(new PointCloudT);
        for (auto& idx:group) {
            foundPts->points.push_back(input->points[idx]);
        }
        output.push_back(foundPts);
	}
}

void ptsToLine(PointCloudT::Ptr input, Eigen::VectorXf& paras, EdgeLine& output) {
	// for simplicity, only consider two dimension
	// need to improve

	PCL_WARN("@ptsToLine need to be improved, it is shorter than the real length ");
	float k1 = paras[4] / paras[3];
	//float b = paras[1] - k * paras[0];
	float b1 = input->points[0].y - k1 * input->points[0].x;


//	PointT min,max;
//	pcl::getMinMax3D(*input, min, max);

	float minX = INT_MAX, maxX = INT_MIN, minY = INT_MAX, maxY = INT_MIN;
	for (auto &p:input->points) {
		float k2 = -1/k1;
		float b2 = p.y - k2 * p.x;
		Eigen::Matrix2f A;
		Eigen::Vector2f b;
		A << k1,-1, k2, -1;
		b << -b1,-b2;
		Eigen::Vector2f x = A.colPivHouseholderQr().solve(b);
		minX = min (minX , x[0]);
		maxX = max (maxX , x[0]);

		minY = min (minY , x[1]);
		maxY = max (maxY , x[1]);
	}

	PointT p,q;
	p.z = 0; q.z = 0;
	if (k1>0) {
		p.x = minX;
		p.y = minY;
		q.x = maxX;
		q.y = maxY;
	}else {
		p.x = minX;
		p.y = maxY;
		q.x = maxX;
		q.y = minY;
	}


	Eigen::VectorXf coff(6);
	output.paras = coff;
	output.p = p;
	output.q = q;
}

void findLinkedLines(vector<EdgeLine>& edgeLines) {
	vector<PointT> allPts;
	for(auto& l: edgeLines) {
		allPts.push_back(l.p);
		allPts.push_back(l.q);
	}
	int n = allPts.size();
	std::unordered_map<int,int> map;

	for (int i = 0; i < n; i++) {
		float minDist = 100;
		int index = -1;
		for (int j = 0; j < n; j++) {
			if (i==j) continue; // self connot connect self
			if (i%2 == 0 && j == i+1) continue; // avoid self p connect self q
			if (i%2 != 0 && j == i-1) continue;
			if (pcl::geometry::distance( allPts[i], allPts[j] ) < minDist) {
				minDist = pcl::geometry::distance( allPts[i], allPts[j] ); // record the min dist
				index = j; // record index which has min dist
			}
		}

		// make sure dist between two points less than 1m
		if (minDist < 1) {
			assert(index >= 0);
			if (map.count(index) && map[index] == i) continue;
			EdgeLine line; // construct new line
			line.p = allPts[i];
			line.q = allPts[index];
			edgeLines.push_back(line);
			map[i] = index;
		}

	}
}

void calculateMeanStandardDev(vector<float> v, float& mean, float& stdev) {
    float sum = std::accumulate(v.begin(), v.end(), 0.0);
    mean = sum / v.size();

    std::vector<float> diff(v.size());
    std::transform(v.begin(), v.end(), diff.begin(), [mean](double x) { return x - mean; });
    float sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    stdev = std::sqrt(sq_sum / v.size());
}

void regionGrow(PointCloudT::Ptr input, int NumberOfNeighbours, int SmoothnessThreshold, int CurvatureThreshold,
		int MinSizeOfCluster, int KSearch, vector<PointCloudT::Ptr>& outputClusters) {
	std::vector<pcl::PointIndices> clustersIndices;
	pcl::search::Search<PointT>::Ptr tree = boost::shared_ptr<pcl::search::Search<PointT> >(
			new pcl::search::KdTree<PointT>);
	pcl::PointCloud<pcl::Normal>::Ptr normals_all(new pcl::PointCloud<pcl::Normal>);
	calculateNormals(input, normals_all, KSearch);
	pcl::RegionGrowing<PointT, pcl::Normal> reg;
	reg.setMinClusterSize(0);
	reg.setMaxClusterSize(100000);
	reg.setSearchMethod(tree);
	reg.setNumberOfNeighbours(NumberOfNeighbours);
	reg.setInputCloud(input);
	reg.setInputNormals(normals_all);
	reg.setSmoothnessThreshold(static_cast<float>(SmoothnessThreshold / 180.0 * M_PI));
	reg.setCurvatureThreshold(CurvatureThreshold);
	reg.extract(clustersIndices);
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr colored_cloud = reg.getColoredCloud();
	for (size_t i = 0; i < clustersIndices.size(); ++i) {
		if (clustersIndices[i].indices.size() < MinSizeOfCluster) continue;
		PointCloudT::Ptr singleCluster(new PointCloudT);
		for (auto &p : clustersIndices[i].indices) {
			singleCluster->points.push_back(input->points[p]);
		}
		outputClusters.push_back(singleCluster);
	}
}

void calculateNormals(PointCloudT::Ptr input, pcl::PointCloud <pcl::Normal>::Ptr &normals_all, int KSearch)
{
	//1-1. generating the normal for each point
    stringstream ss;
    cout << "The point you input doesn't contain normals, calculating normals...\n";
    pcl::search::Search<PointT>::Ptr tree = boost::shared_ptr<pcl::search::Search<PointT> >(new pcl::search::KdTree<PointT>);
    pcl::NormalEstimation<PointT, pcl::Normal> normal_estimator;
    normal_estimator.setSearchMethod(tree);
    normal_estimator.setInputCloud(input);
    normal_estimator.setKSearch(KSearch);
    normal_estimator.compute(*normals_all);

    for (size_t i = 0; i < normals_all->points.size(); ++i)
    {
        input->points[i].normal_x = normals_all->points[i].normal_x;
        input->points[i].normal_y = normals_all->points[i].normal_y;
        input->points[i].normal_z = normals_all->points[i].normal_z;
    }
}


void removePtsAroundLine(PointCloudT::Ptr input, PointCloudT::Ptr output, vector<EdgeLine>& lines, float dist){
	PointCloudT::Ptr allPts(new PointCloudT);
	copyPointCloud(*input, *allPts);
	for (auto &p : allPts->points) p.z = 0;
	vector<PointCloudT::Ptr> linesPts;
	for (auto &line: lines) {
		PointCloudT::Ptr tmp(new PointCloudT);
		generateLinePointCloud(line.p, line.q, 100, 255 <<24 | 1255 << 16 | 0 << 8 | 0, tmp);
		linesPts.push_back(tmp);
		for (auto& p: tmp->points) allPts->push_back(p);
	}

	simpleView("[removePtsAroundLine] allPts" ,allPts);

	PointCloudT::Ptr filtered(new PointCloudT);
	VoxelGrid<PointT> sor;
	sor.setInputCloud(allPts);
	sor.setLeafSize(dist,dist,20);
	sor.filter(*filtered);
	Eigen::Vector3i divisions = sor.getNrDivisions();
	cout << divisions.transpose() << endl;
	unordered_set<int> filterGridIdx;
	for (auto &linePts : linesPts) {
		for (auto &p : linePts->points) {
			auto cord = sor.getGridCoordinates(p.x, p.y, p.z);
			cord = cord - sor.getMinBoxCoordinates();
			filterGridIdx.insert(cord[0]*divisions[1]*divisions[2] + cord[1]*divisions[2] + cord[2]);
		}
	}

	for (auto &p : input->points) {
		auto cord = sor.getGridCoordinates(p.x, p.y, p.z);
		cord = cord - sor.getMinBoxCoordinates();
		if (filterGridIdx.count(cord[0]*divisions[1]*divisions[2] + cord[1]*divisions[2] + cord[2])) continue;
		output->push_back(p);
	}
	cout << "[removePtsAroundLine] before " << input->size() << "  after " << output->size() << endl;
	simpleView("[removePtsAroundLine] ",output);
}

void regionGrow2D(PointCloudT::Ptr input, vector<PointCloudT::Ptr>& output){
    // calculate Normal
    PointCloud<Normal>::Ptr cloud_normals(new PointCloud<Normal>);
    calculateNormal2D(input, cloud_normals);

    // 2 thresholds for region grow
    float radius = Paras<float>("RegionGrow2D","serarchRadius");
    int normalTh = Paras<int>("RegionGrow2D","normalTh");
    pcl::KdTreeFLANN<PointT> kdtree;
    kdtree.setInputCloud(input);
    unordered_set<int> found;
    int k = 0;

    // start from the seed point, and check normal difference between seed points and neighbor points is less than thresholds.
    while(found.size() < input->size()) {
        vector<int> group;
        stack<PointT> seeds;
        while(found.count(k)) k++;
        if (k >= input->size()) break;
        seeds.push(input->points[k]);
        found.insert(k);
        group.push_back(k);
        Eigen::Vector4f seedNormal = cloud_normals->points[k].getNormalVector4fMap();
        while(!seeds.empty()) {
            PointT seed = seeds.top();
            seeds.pop();
            vector<int> pointIdx;
            vector<float> dist;
            kdtree.radiusSearch(seed, radius, pointIdx,dist);
            for(auto& ix : pointIdx) {
                if (found.count(ix)) continue;
                Eigen::Vector4f q = cloud_normals->points[ix].getNormalVector4fMap();
                auto r = seedNormal.dot(q);
                float angle_in_radiant = acos(r);
                if (angle_in_radiant < M_PI/180*normalTh) {
                    found.insert(ix);
                    seeds.push(input->points[ix]);
                    group.push_back(ix);
                }
            }
        }
        if (group.size() < 2) continue;
        PointCloudT::Ptr foundPts(new PointCloudT);
        for (auto& idx:group) {
            foundPts->points.push_back(input->points[idx]);
        }
        output.push_back(foundPts);
    }

    cout << "[regionGrow2D] " << output.size() << " clusters are found" << endl;
}

void calculateNormal2D(PointCloudT::Ptr input, PointCloud<Normal>::Ptr cloud_normals){
    // using SVD to calculate the normal
    KdTreeFLANN<PointT> kdtree;
    kdtree.setInputCloud(input);
    // K nearest points are used to compute matrix
    int K = 10;
    for (auto &p:input->points) {
        Eigen::Vector3f v; v << p.x, p.y, p.z ;
        std::vector<int> pointIdxNKNSearch(K);
        std::vector<float> pointNKNSquaredDistance(K);
        Eigen::MatrixXf conv = Eigen::MatrixXf::Zero(3,3);
        Normal normal; normal.normal_x = 0; normal.normal_y = 0; normal.normal_z = 0;
        if ( kdtree.nearestKSearch (p, K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0 )
        {
            for (int & idx: pointIdxNKNSearch) {
                Eigen::Vector3f x_i;
                x_i << input->points[idx].x, input->points[idx].y, input->points[idx].z;
                conv = conv +  ( (v - x_i) * (v - x_i).transpose() );
            }
            Eigen::JacobiSVD<Eigen::MatrixXf> svd(conv,  Eigen::ComputeThinU | Eigen::ComputeThinV);
            // the second normal vector is perpendicular with line
            normal.normal_x = svd.matrixU()(1,0);
            normal.normal_y = svd.matrixU()(1,1);
            normal.normal_z = svd.matrixU()(1,2);
            // cout << p.x << "," << p.y << "," << p.z << "," << svd.matrixU()(1,0) << ", " << svd.matrixU()(1,1) << ", " << svd.matrixU()(1,2) << ";\n";
        }
        cloud_normals->push_back(normal);
    }
}

void BeamRANSAC(PointCloudT::Ptr input, float high, vector< pair < vector<PointT>, float >> &planes,
				vector<pair<EdgeLine, float>> &lines){
	float possibleHeight = Paras<float>("BeamRecons","possibleHeight");; // we consider possible height of

	PointCloudT::Ptr inputTmp(new PointCloudT);
	pcl::copyPointCloud(*input,*inputTmp);

	// filter the input point clouds
	pcl::PassThrough<PointT> filterZ;
	filterZ.setInputCloud(inputTmp);
	filterZ.setFilterFieldName("z");
	filterZ.setFilterLimits(high-possibleHeight, high - 0.1);
	filterZ.filter(*inputTmp);
	simpleView("[BeamRANSAC] extract Z", inputTmp);

	vector<PointCloudT::Ptr> outputClusters;
	// use region grow to separete into different clusters
	regionGrow(inputTmp, 30, 2, 2, 30, 10, outputClusters);

	vector<PointCloudT::Ptr> possibleBeamPlanes;
	vector<Eigen::Vector4f> possibleBeamPlanesCoff;
	for (auto& cluster: outputClusters) {
	    // run RANSAC to check each cluster is possible beam or not
		pcl::ModelCoefficients::Ptr sacCoefficients(new pcl::ModelCoefficients);
		pcl::PointIndices::Ptr sacInliers(new pcl::PointIndices);
		pcl::SACSegmentation<PointT> seg;
		seg.setOptimizeCoefficients(true);
		seg.setModelType(pcl::SACMODEL_PLANE);
		seg.setMethodType(pcl::SAC_RANSAC);
		seg.setDistanceThreshold(0.05);
		seg.setInputCloud(cluster);
		seg.segment(*sacInliers, *sacCoefficients);
		Eigen::Vector4f  coff;
		coff << sacCoefficients->values[0], sacCoefficients->values[1], sacCoefficients->values[2], sacCoefficients->values[3];
		// we only need to remain the planes which are very parallel with xy plane
		if (abs(coff[0]) < 0.1 && abs(coff[1]) < 0.1){
			possibleBeamPlanes.push_back(cluster);
			possibleBeamPlanesCoff.push_back(coff);
		}
	}
	cout << "[BeamRANSAC] " << possibleBeamPlanes.size() << " are found from " << outputClusters.size() << endl;

	// next step, we want to construct lines around found points and check these lines up space whether contain pts
	//simpleView("[BeamRANSAC] possible Beam Planes", possibleBeamPlanes);
	//vector<PointCloudT::Ptr> tmp;

	for (int i = 0; i < possibleBeamPlanes.size(); ++i){
		PointCloudT::Ptr resTmp(new PointCloudT);
		// find 4 edges for each plane
		vector<PointT> edges;
		bool isBeam = findEdgeForPlane(possibleBeamPlanes[i], high, inputTmp, edges);
		// construct one horizontal plane and 4 vertical planes
		vector<PointT> plane{edges[0],edges[1],edges[2],edges[3]};
		planes.push_back( make_pair(plane,-possibleBeamPlanesCoff[i][3]));
		// if it is beam, construct 4 planes
		if (isBeam) {
			pair<EdgeLine, float> l1,l2,l3,l4; // construct line and height
			EdgeLine line1; line1.p = edges[0]; line1.q = edges[1]; l1 = make_pair(line1,-possibleBeamPlanesCoff[i][3]);
			EdgeLine line2; line2.p = edges[1]; line2.q = edges[2]; l2 = make_pair(line2,-possibleBeamPlanesCoff[i][3]);
			EdgeLine line3; line3.p = edges[2]; line3.q = edges[3]; l3 = make_pair(line3,-possibleBeamPlanesCoff[i][3]);
			EdgeLine line4; line4.p = edges[3]; line4.q = edges[0]; l4 = make_pair(line4,-possibleBeamPlanesCoff[i][3]);
			lines.push_back(l1); lines.push_back(l2); lines.push_back(l3); lines.push_back(l4);
		}

/*		generateLinePointCloud(edges[0],edges[1],100, INT32_MAX, resTmp);
		generateLinePointCloud(edges[1],edges[2],100, INT32_MAX, resTmp);
		generateLinePointCloud(edges[2],edges[3],100, INT32_MAX, resTmp);
		generateLinePointCloud(edges[3],edges[0],100, INT32_MAX, resTmp);
		tmp.push_back(resTmp);*/
	}
	//for (auto& t : tmp) possibleBeamPlanes.push_back(t);

	simpleView("[BeamRANSAC] possible Beam Planes with check boundary", possibleBeamPlanes);

}

bool findEdgeForPlane(PointCloudT::Ptr input, float maxRoomZ, PointCloudT::Ptr roofPart, vector<PointT>& edges){
	bool isBeam = true;
    // 1 - determine the direction of the plane
	pcl::ModelCoefficients::Ptr sacCoefficients(new pcl::ModelCoefficients);
	pcl::PointIndices::Ptr sacInliers(new pcl::PointIndices);
	pcl::SACSegmentation<PointT> seg;
	seg.setOptimizeCoefficients(true);
	seg.setModelType(pcl::SACMODEL_LINE);
	seg.setMethodType(pcl::SAC_RANSAC);
	seg.setDistanceThreshold(Paras<float>("BeamRecons","RANSACDistTh"));
	seg.setInputCloud(input);
	seg.segment(*sacInliers, *sacCoefficients);
	Eigen::VectorXf coff(6);
	coff << sacCoefficients->values[0], sacCoefficients->values[1], sacCoefficients->values[2],
			sacCoefficients->values[3], sacCoefficients->values[4], sacCoefficients->values[5];
	//cout << coff[3]<< " " << coff[4] << " " << coff[5]<< " " <<endl;

	// 2 - rotate the plane according to the angle found
	PointCloudT::Ptr transformed(new PointCloudT);
	Eigen::Matrix4f transform = Eigen::Matrix4f::Identity();
	float theta = acos(abs(coff[4]) / sqrt(coff[3]*coff[3] + coff[4]*coff[4] + coff[5]*coff[5]));
	if (coff[3]*coff[4] < 0) theta = -theta;
	transform(0,0) = cos(theta); transform(0,1) = -sin(theta); transform(1,0) = sin(theta); transform(1,1) = cos(theta);
	transformPointCloud(*input, *transformed, transform);

    PointCloudT::Ptr roofPartTrans(new PointCloudT);
    transformPointCloud(*roofPart, *roofPartTrans, transform);
	//simpleView("[findEdgeForPlane] transformed pts", transformed);

    // 3 - define edges : transform so it can be parallel with x or y axis
	PointT minP, maxP, p1, p2, p3, p4;
	getMinMax3D(*transformed, minP, maxP);
	p1.x = minP.x; p1.y = minP.y; p1.z = maxP.z;
	p2.x = minP.x; p2.y = maxP.y; p2.z = maxP.z;
	p3.x = maxP.x; p3.y = maxP.y; p3.z = maxP.z;
	p4.x = maxP.x; p4.y = minP.y; p4.z = maxP.z;

	// 4 - determine the points in the rectangular shape (to verify whether it is a beam)
    ConditionAnd<PointT>::Ptr range_cond (new ConditionAnd<PointT> ());
    range_cond->addComparison (FieldComparison<PointT>::ConstPtr (new FieldComparison<PointT> ("x", ComparisonOps::GT, minP.x)));
    range_cond->addComparison (FieldComparison<PointT>::ConstPtr (new FieldComparison<PointT> ("x", ComparisonOps::LT, maxP.x)));
    range_cond->addComparison (FieldComparison<PointT>::ConstPtr (new FieldComparison<PointT> ("y", ComparisonOps::GT, minP.y)));
    range_cond->addComparison (FieldComparison<PointT>::ConstPtr (new FieldComparison<PointT> ("y", ComparisonOps::LT, maxP.y)));
    range_cond->addComparison (FieldComparison<PointT>::ConstPtr (new FieldComparison<PointT> ("z", ComparisonOps::GT, maxP.z+0.05)));
    range_cond->addComparison (FieldComparison<PointT>::ConstPtr (new FieldComparison<PointT> ("z", ComparisonOps::LT, maxRoomZ-0.05)));
    //simpleView("roofPartTrans", roofPartTrans);

    ConditionalRemoval<PointT> condrem;
    condrem.setCondition (range_cond);
    condrem.setInputCloud (roofPartTrans);
	PointCloudT::Ptr filtered(new PointCloudT);
    condrem.filter (*filtered);
    // todo: we can calculate how many points insides the cube. if it is more than a certain num, then it is not beam
    //cout << "num of points left " << filtered->size() << " from" <<  roofPart->size() << endl;
	//simpleView("trans -  filtered",filtered);

	// 6 - check lines in Z direction whether has num of points
	PointCloudT::Ptr debugCube(new PointCloudT);
	copyPointCloud(*roofPartTrans,*debugCube);
	int numPtsZ1 = checkNumofPointofLineinZ(p1,p2, roofPartTrans, maxP.z, maxRoomZ,debugCube);
	int numPtsZ2 = checkNumofPointofLineinZ(p2,p3, roofPartTrans, maxP.z, maxRoomZ,debugCube);
	int numPtsZ3 = checkNumofPointofLineinZ(p3,p4, roofPartTrans, maxP.z, maxRoomZ,debugCube);
	int numPtsZ4 = checkNumofPointofLineinZ(p4,p1, roofPartTrans, maxP.z, maxRoomZ,debugCube);
	if (Paras<bool>("View","BeamCube")) simpleView("[findEdgeForPlane] BeamCube", debugCube);
	// test Part
/*	PointT testA,testB,testC,testD;
	testA.x = minP.x; testA.y = minP.y; testA.z = maxP.z;
	testB.x = minP.x; testB.y = maxP.y; testB.z = maxP.z;
	testC.x = maxP.x; testC.y = maxP.y; testC.z = maxP.z;
	testD.x = maxP.x; testD.y = minP.y; testD.z = maxP.z;
	generateLinePointCloud(testA ,testB, 100, INT32_MAX, roofPartTrans);
	generateLinePointCloud(testB ,testC, 100, INT32_MAX, roofPartTrans);
	generateLinePointCloud(testC ,testD, 100, INT32_MAX, roofPartTrans);
	generateLinePointCloud(testD ,testA, 100, INT32_MAX, roofPartTrans);
	simpleView("trans -  with 4 edges",roofPartTrans);*/

	// 7 - check whether it is valid beam by check ratio of num. pts over cube size
	float distTh = Paras<float>("BeamRecons","EdgeZDist");
	assert( abs(p1.z-p2.z) < 0.1 && abs(p2.z-p3.z) < 0.1 && abs(p3.z-p4.z) < 0.1 && abs(p4.z-p1.z) < 0.1);
	float cubeSize1 = (geometry::distance(p1,p2)) * (maxRoomZ - p1.z) *  distTh * 2 * 1000;
	float cubeSize2 = (geometry::distance(p2,p3)) * (maxRoomZ - p2.z) *  distTh * 2 * 1000;
	float cubeSize3 = (geometry::distance(p3,p4)) * (maxRoomZ - p3.z) *  distTh * 2 * 1000;
	float cubeSize4 = (geometry::distance(p4,p1)) * (maxRoomZ - p4.z) *  distTh * 2 * 1000;
	assert( cubeSize1 > 0 && cubeSize2 > 0 && cubeSize3 > 0 && cubeSize4 > 0);
    //cout << numPtsZ1/cubeSize1 << " "<< numPtsZ2/cubeSize2 << " "<< numPtsZ3/cubeSize3 << " "<< numPtsZ4/cubeSize4 << endl;
	// just make sure this way works in different environment
	assert( ((numPtsZ1 > -1) + (numPtsZ2 > -1) + (numPtsZ3 > -1) + (numPtsZ4 > -1)) == 4 );
	if ( ( (numPtsZ1/cubeSize1 >= Paras<float>("BeamRecons","cubeNumPtsRatio")) +
		   (numPtsZ2/cubeSize2 >= Paras<float>("BeamRecons","cubeNumPtsRatio")) +
		   (numPtsZ3/cubeSize3 >= Paras<float>("BeamRecons","cubeNumPtsRatio")) +
		   (numPtsZ4/cubeSize4 >= Paras<float>("BeamRecons","cubeNumPtsRatio"))) < 2 ) isBeam = false;
	// 5 - reverse
	Eigen::Affine3f reverseT = Eigen::Affine3f::Identity();
	reverseT.rotate(Eigen::AngleAxisf(-theta,Eigen::Vector3f::UnitZ()));
	p1 = transformPoint(p1,reverseT);
	p2 = transformPoint(p2,reverseT);
	p3 = transformPoint(p3,reverseT);
	p4 = transformPoint(p4,reverseT);

/*	generateLinePointCloud(p1,p2,100, INT32_MAX, transformed);
	generateLinePointCloud(p2,p3,100, INT32_MAX, transformed);
	generateLinePointCloud(p3,p4,100, INT32_MAX, transformed);
	generateLinePointCloud(p4,p1,100, INT32_MAX, transformed);*/
	edges = vector<PointT>{p1,p2,p3,p4};
	return isBeam;

}

int checkNumofPointofLineinZ(PointT p, PointT q, PointCloudT::Ptr input, float minZ, float maxZ, PointCloudT::Ptr debug_cube){

	float distTh = Paras<float>("BeamRecons","EdgeZDist");  // distance threshold
	float k = (p.y - q.y) / (p.x - q.x);
	float theta = -atan(k);

	// transform points to be parallel with x-axis
	PointCloudT::Ptr transformed(new PointCloudT);
	Eigen::Matrix4f transform = Eigen::Matrix4f::Identity();
	transform(0,0) = cos(theta); transform(0,1) = -sin(theta); transform(1,0) = sin(theta); transform(1,1) = cos(theta);
	transformPointCloud(*input, *transformed, transform);

	PointT r_p, r_q;
	Eigen::Affine3f T = Eigen::Affine3f::Identity();
	T.rotate(Eigen::AngleAxisf(theta,Eigen::Vector3f::UnitZ()));
	r_p = transformPoint(p,T); r_q = transformPoint(q,T);

	//cout << p.y << " " << q.y << " " << p.x << " " << q.x << "  " << k << " theta:" << theta*180/M_PI <<  endl;
	if (abs(r_p.y-r_q.y) > 0.1) { cout << "r_p.y:" << r_p.y << "  r_q.y:" << r_q.y  << endl; abort();}
	ConditionAnd<PointT>::Ptr range_cond (new ConditionAnd<PointT> ());

	// set the box to count number of points within in this box
	range_cond->addComparison (FieldComparison<PointT>::ConstPtr (new FieldComparison<PointT> ("x", ComparisonOps::GT, min(r_p.x,r_q.x))));
	range_cond->addComparison (FieldComparison<PointT>::ConstPtr (new FieldComparison<PointT> ("x", ComparisonOps::LT, max(r_p.x,r_q.x))));
	range_cond->addComparison (FieldComparison<PointT>::ConstPtr (new FieldComparison<PointT> ("y", ComparisonOps::GT, r_p.y-distTh)));
	range_cond->addComparison (FieldComparison<PointT>::ConstPtr (new FieldComparison<PointT> ("y", ComparisonOps::LT, r_p.y+distTh)));
	range_cond->addComparison (FieldComparison<PointT>::ConstPtr (new FieldComparison<PointT> ("z", ComparisonOps::GT, minZ)));
	range_cond->addComparison (FieldComparison<PointT>::ConstPtr (new FieldComparison<PointT> ("z", ComparisonOps::LT, maxZ)));
	ConditionalRemoval<PointT> condrem;
	condrem.setCondition (range_cond);
	condrem.setInputCloud (transformed);
	PointCloudT::Ptr filtered(new PointCloudT);
	condrem.filter (*filtered);

	// debug part, to visualize the result
	if (Paras<bool>("View","BeamCube")) {
		PointT a = r_p, b = r_q, c = r_q, d = r_p, a2 = r_p, b2 = r_q, c2 = r_q, d2 = r_p;
		a.y -= 0.1;  b.y -= 0.1; c.y += 0.1;  d.y += 0.1;
		a.z = minZ; b.z =minZ; c.z = minZ; d.z = minZ;
		a2.y -= 0.1;  b2.y -= 0.1; c2.y += 0.1;  d2.y += 0.1;
		a2.z = maxZ; b2.z =maxZ; c2.z = maxZ; d2.z = maxZ;
		Eigen::Affine3f reverseT = Eigen::Affine3f::Identity();
		reverseT.rotate(Eigen::AngleAxisf(-theta,Eigen::Vector3f::UnitZ()));
		r_p = transformPoint(r_p,reverseT); r_q = transformPoint(r_q,reverseT);
		a = transformPoint(a,reverseT); b = transformPoint(b,reverseT);c = transformPoint(c,reverseT); d = transformPoint(d,reverseT);
		a2 = transformPoint(a2,reverseT); b2 = transformPoint(b2,reverseT);c2 = transformPoint(c2,reverseT); d2 = transformPoint(d2,reverseT);
		generateLinePointCloud(r_p,r_q,100,Colors.White, debug_cube);
		generateLinePointCloud(a,b,c,d,100,Colors.Red, debug_cube);
		generateLinePointCloud(a2,b2,c2,d2,100,Colors.Blue, debug_cube);
		if (Paras<bool>("View","BeamCubeEach")) simpleView("[checkNumofPointofLineinZ] test check num of pts inside cube", debug_cube);
	}

	return filtered->size();

}

void exportToDxf(string outputPath, const vector<EdgeLine>& wallLines, const vector< pair<vector<PointT>, float> > &beamPlanes,
				 const vector<pair<EdgeLine, float>> &beamLines, float minZ, float maxZ){
    KKRecons::DxfExporter exporter;
    for (auto& line: wallLines) {
        PointXYZ a,b,c,d;
        a.x = line.p.x; a.y = line.p.y; a.z = minZ;
        b.x = line.q.x; b.y = line.q.y; b.z = minZ;
        c.x = line.q.x; c.y = line.q.y; c.z = maxZ;
        d.x = line.p.x; d.y = line.p.y; d.z = maxZ;
        DxfFace face(a,b,c,d);
        exporter.insert(face);
    }

	for (auto& plane: beamPlanes) {
		PointXYZ a,b,c,d;
		auto pts_v = plane.first;
		a.x = pts_v[0].x; a.y = pts_v[0].y; a.z = plane.second;
		b.x = pts_v[1].x; b.y = pts_v[1].y; b.z = plane.second;
		c.x = pts_v[2].x; c.y = pts_v[2].y; c.z = plane.second;
		d.x = pts_v[3].x; d.y = pts_v[3].y; d.z = plane.second;
		DxfFace face(a,b,c,d);
		exporter.insert(face);
	}

	for (auto& bline: beamLines) {
		auto line = bline.first;
		PointXYZ a,b,c,d;
		a.x = line.p.x; a.y = line.p.y; a.z = bline.second;
		b.x = line.q.x; b.y = line.q.y; b.z = bline.second;
		c.x = line.q.x; c.y = line.q.y; c.z = maxZ;
		d.x = line.p.x; d.y = line.p.y; d.z = maxZ;
		DxfFace face(a,b,c,d);
		exporter.insert(face);
	}


    exporter.exportDXF(outputPath);
}

void minimumCut(const vector<EdgeLine>& lines, PointCloudT::Ptr input){

	Arrangement_2   arr;

	PointT minP, maxP;
	getMinMax3D(*input, minP, maxP);

	// 1. set the min-cut boundary
	PointT P_a, P_b, P_c, P_d;
    P_a.x = minP.x - 2; P_a.y = minP.y - 2; P_a.z = 0;
    P_b.x = minP.x - 2; P_b.y = maxP.y + 2; P_b.z = 0;
    P_c.x = maxP.x + 2; P_c.y = maxP.y + 2; P_c.z = 0;
    P_d.x = maxP.x + 2; P_d.y = minP.y - 2; P_d.z = 0;
	KKRecons::DebugFileExporter lineOutput("./CGAL_Line_Output.txt");
	KKRecons::DebugFileExporter wallOutput("./CGAL_Wall_Output.txt");

	int kk = 1;
	for (auto& line: lines) {
		PointT p,q;
		p = line.p;
		q = line.q;
		float k = -1 / ((p.y - q.y) / (p.x - q.x));
		float b1 = p.y - (k*p.x);
		float b2 = q.y - (k*q.x);
		PointT line1_s1, line1_s2, line2_s1, line2_s2; // two lines are vertical with the edge
		line1_s1.x = minP.x-2; line1_s1.y = k*line1_s1.x + b1; line1_s1.z = 0;
		line1_s2.x = maxP.x+2; line1_s2.y = k*line1_s2.x + b1; line1_s2.z = 0;

		line2_s1.x = minP.x-2; line2_s1.y = k*line2_s1.x + b2; line2_s1.z = 0;
		line2_s2.x = maxP.x+2; line2_s2.y = k*line2_s2.x + b2; line2_s2.z = 0;
		//if (k < 0 ) continue;
		Segment_2 s0(Point_2(p.x, p.y), Point_2(q.x, q.y));
		Segment_2 s1(Point_2(line1_s1.x, line1_s1.y), Point_2(line1_s2.x, line1_s2.y));
		Segment_2 s2(Point_2(line2_s1.x, line2_s1.y), Point_2(line2_s2.x, line2_s2.y));
		insert(arr, s0);
		insert(arr, s1);
		insert(arr, s2);

//		int test = 0;
//		if (isIntersect(line1_s1,line1_s2,P_a,P_b)) test++;
//        if (isIntersect(s1,s2,P_b,P_c)) test++;
//        if (isIntersect(s1,s2,P_c,P_d)) test++;
//        if (isIntersect(s1,s2,P_a,P_d)) test++;
//        if (test!=2) {
//        	cout << test << endl;
//        	abort();
//        }
        //if (kk++ > 10) break;
		string ss;
		ss = to_string(line1_s1.x) + " " + to_string(line1_s1.y) + " " + to_string(line1_s2.x) + " " + to_string(line1_s2.y);
		lineOutput.insertLine( ss );
		ss = to_string(line2_s1.x) + " " + to_string(line2_s1.y) + " " + to_string(line2_s2.x) + " " + to_string(line2_s2.y);
		lineOutput.insertLine( ss );
		ss = to_string(p.x) + " " + to_string(p.y) + " " + to_string(q.x) + " " + to_string(q.y);
		wallOutput.insertLine( ss );
	}
	lineOutput.exportToPath();
	wallOutput.exportToPath();
	Arrangement_2::Vertex_const_iterator  vit;
	std::cout << arr.number_of_vertices() << " vertices:" << std::endl;
	KKRecons::DebugFileExporter CGALOutput("./CGAL_Vertex_Output.txt");
	for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit)
	{
		Traits_2::Point_2 point = vit->point();
		double x = CGAL::to_double(point.x());
		double y = CGAL::to_double(point.y());
		string ss;
		ss = to_string(x) + " " + to_string(y);
		CGALOutput.insertLine(ss);
	}
	CGALOutput.exportToPath();
    vector<CGALFace> allFaces;
    get_arrangement_face(arr, allFaces);
    cout << "[minimumCut] num of faces: " << allFaces.size() << endl;

    unordered_set<CGALFace> queryFaces;
	PointLocation query;
	query.attach(arr);
	int time = 0;
	cout << "input size " << input->size() << endl;
	for (auto& p : input->points) {
		Point_2 _p(p.x, p.y);
		vector<Point_2> queryRes;
		queryRes = locate_point(query, _p);
		if (queryRes.size() != 0) {
			CGALFace face;
			face.vertex = queryRes;
			queryFaces.insert(face);
		}
		if ((time++ % 15000) == 0) cout << time << endl;
	}

	cout << "size of queryFaces " << queryFaces.size() << endl;

	KKRecons::DebugFileExporter facesOutput("./CGAL_Face_Output.txt");
	for (auto& face : allFaces) {
		if (queryFaces.count(face)) {
			face.type = 1;
		}

		string ss = "";
		for (auto &v: face.vertex) {
			ss = ss + to_string(CGAL::to_double(v.x())) + " " + to_string(CGAL::to_double(v.y())) + " ";
		}

		ss += to_string(face.type);
		facesOutput.insertLine(ss);
	}
	facesOutput.exportToPath();
}


void get_arrangement_face (CGAL::Arrangement_2<Traits_2>& arr, vector<CGALFace> &faces)
{

	CGAL_precondition (arr.is_valid());
	typename Arrangement_2::Face_const_iterator    fit;

    assert(faces.size() == 0);
	std::cout << arr.number_of_faces() << " faces:" << std::endl;
	for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
		CGALFace face;
		if (fit->is_unbounded()) std::cout << "Unbounded face. " << std::endl;
		else {
			Arrangement_2::Ccb_halfedge_const_circulator circ = fit->outer_ccb();
			Arrangement_2::Ccb_halfedge_const_circulator curr = circ;
			Arrangement_2::Halfedge_const_handle          he;
			string ss = "";
			do
			{
				he = curr;
				ss = ss + to_string(CGAL::to_double(he->target()->point().x())) + " " + to_string(CGAL::to_double(he->target()->point().y())) + " ";
				Point_2 p(he->target()->point().x(), he->target()->point().y());
				face.vertex.push_back(p);
				++curr;
			} while (curr != circ);

			faces.push_back(face);
		}

	}


	return;
}

vector<Point_2> locate_point(const PointLocation& pl,
				  const PointLocation::Arrangement_2::Point_2& q)
{
	typedef PointLocation                                 Point_location;
	typedef typename Point_location::Arrangement_2        Arrangement_2;
	typename CGAL::Arr_point_location_result<Arrangement_2>::Type obj =
			pl.locate(q);
	return print_point_location(q, obj);
}

vector<Point_2> print_point_location
		(const Arrangement_2::Point_2& q,
		 typename CGAL::Arr_point_location_result<Arrangement_2>::Type obj)
{
	typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
	typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
	typedef typename Arrangement_2::Face_const_handle     Face_const_handle;
	const Face_const_handle*     f;
	vector<Point_2> vertex;
	if (f = boost::get<Face_const_handle>(&obj)) {
		if ((*f)->is_unbounded()) return vertex;
		else{
			Arrangement_2::Ccb_halfedge_const_circulator  circ =  (*f)->outer_ccb();
			Arrangement_2::Ccb_halfedge_const_circulator  curr = circ;
			Arrangement_2::Halfedge_const_handle          he;
			//string ss = "";
			do
			{
				he = curr;
				Point_2 p(he->target()->point().x(), he->target()->point().y());
				//ss = ss + to_string(CGAL::to_double(he->target()->point().x())) + " " + to_string(CGAL::to_double(he->target()->point().y())) + " ";
				vertex.push_back(p);
				++curr;
			} while (curr != circ);
			//cout << ss << endl;
		}
	}
	return vertex;
}