//
// Created by czh on 11/4/18.
//

#ifndef TEST_PCL_EXTRACTWALL_H
#define TEST_PCL_EXTRACTWALL_H

#include <iostream>
#include <vector>
#include <pcl/point_types.h>
#include <pcl/search/search.h>
#include <pcl/search/kdtree.h>
#include <pcl/features/normal_3d.h>
#include <pcl/segmentation/region_growing.h>
#include <pcl/segmentation/sac_segmentation.h> //RANSAC‚Ì‚½‚ß
#include <math.h>
#include <pcl/filters/project_inliers.h> //•½–Ê‚É“Š‰e‚·‚é‚½‚ß
#include <pcl/filters/passthrough.h>
#include <pcl/common/pca.h>
#include <pcl/common/common.h>
#include <pcl/filters/voxel_grid.h> //ƒ_ƒEƒ“ƒTƒ“ƒvƒŠƒ“ƒO‚Ì‚½‚ß
#include <pcl/common/geometry.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/filters/conditional_removal.h>
#include "Plane.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arr_simple_point_location.h>
using namespace std;
using namespace pcl;

typedef pcl::PointXYZRGBNormal PointT;
typedef pcl::PointCloud<PointT> PointCloudT;

typedef CGAL::Quotient<CGAL::MP_Float>             Number_type;
typedef CGAL::Cartesian<Number_type>                  Kernel;
typedef CGAL::Arr_linear_traits_2<Kernel>            Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::X_monotone_curve_2                 Segment_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;
typedef Arrangement_2::Vertex_handle                  Vertex_handle;
typedef Arrangement_2::Halfedge_handle                Halfedge_handle;
typedef CGAL::Arr_simple_point_location<Arrangement_2> PointLocation;

struct reconstructParas
{
    // Downsampling
    int KSearch = 10;
    float leafSize = 0.05; // unit is meter -> 5cm

    // Clustering
    int MinSizeOfCluster = 800;
    int NumberOfNeighbours = 30;
    int SmoothnessThreshold = 5; // angle 360 degree
    int CurvatureThreshold = 10;
    // RANSAC
    double RANSAC_DistThreshold = 0.25; //0.25;
    float RANSAC_MinInliers = 0.5; // 500 todo: should be changed to percents
    float RANSAC_PlaneVectorThreshold = 0.2;

    // Fill the plane
    int pointPitch = 20; // number of point in 1 meter

    // Combine planes
    int minimumEdgeDist = 1; //we control the distance between two edges and the height difference between two edges
    float minHeightDiff = 0.5;
    int minAngle_normalDiff = 5;// when extend smaller plane to bigger plane, we will calculate the angle between normals of planes

    int roof_NumberOfNeighbours = 1;
    int roof_SmoothnessThreshold = 2;
    int roof_CurvatureThreshold = 1;
    int roof_MinSizeOfCluster = 1;
};

struct EdgeLine{
    Eigen::VectorXf paras;
    PointT p;
    PointT q;
};

struct CGALFace{
    vector<Point_2> vertex;
    vector<pair<float,float>> vertexOrigin;
    int type = 0;
    bool operator==(const CGALFace& f) const {
        return vertex == f.vertex;
    }
};

struct Colors{
    const int32_t Red    = 255 <<24 | 255 << 16 | 0   << 8 | 0;
    const int32_t Yellow = 255 <<24 | 255 << 16 | 255 << 8 | 0;
    const int32_t Blue   = 255 <<24 | 0   << 16 | 0   << 8 | 255;
    const int32_t Green  = 255 <<24 | 0   << 16 | 255 << 8 | 0;
    const int32_t Cyan   = 255 <<24 | 0   << 16 | 255 << 8 | 255;
    const int32_t Pink   = 255 <<24 | 255 << 16 | 0   << 8 | 255;
    const int32_t White  = INT32_MAX;
}Colors;

// override the hash value in order to use CGALFace as the key of unordered_map
namespace std {
    template <>
    struct hash<CGALFace>
    {
        std::size_t operator()(const CGALFace& f) const
        {
            using std::size_t;
            using std::hash;
            using std::string;

            // Compute individual hash values for first,
            // second and third and combine them using XOR
            // and bit shifting:
            size_t res = 0;
            for (int i = 0; i < f.vertex.size(); ++i) {
                res = ( ( res + ((int)CGAL::to_double(f.vertex[i].x()) << 1) ) >> 1 )+ CGAL::to_double(f.vertex[i].x());
            }
            return res;
        }
    };

}

template <typename Container> // we can make this generic for any container [1]
struct container_hash {
    std::size_t operator()(Container const& c) const {
        return boost::hash_range(c.begin(), c.end());
    }
};


// generate points between pt1 and pt2. pointPitch is the number of points within 1 meter
void generateLinePointCloud(PointT pt1, PointT pt2, int pointPitch, int color, PointCloudT::Ptr output);
void generateLinePointCloud(PointT a, PointT b, PointT c, PointT d, int pointPitch, int color, PointCloudT::Ptr output);

bool onSegment(PointT p, PointT q, PointT r);
float orientation(PointT p, PointT q, PointT r);
// check whether two lines are intersected
bool isIntersect(PointT p1, PointT q1, PointT p2, PointT q2);


void smoothNoise(PointCloudT::Ptr input, int K, float alpha);
void densityFilter(PointCloudT::Ptr input);

void detectHeightRange(PointCloudT::Ptr input, float& high, float& low);

void extractTopPts(PointCloudT::Ptr input, PointCloudT::Ptr output, float highest, float dimension, reconstructParas paras);
void convert2D(PointCloudT::Ptr input,PointCloudT::Ptr output);
void findBiggestComponent2D(PointCloudT::Ptr input, PointCloudT::Ptr output);

void computeHull(PointCloudT::Ptr input, PointCloudT::Ptr output, float alpha);
void extractLineFromEdge(PointCloudT::Ptr input, vector<EdgeLine>& edgeLines);
void extractLineFromEdge2(PointCloudT::Ptr input, vector<EdgeLine>& edgeLines);
void separatePtsToGroups(PointCloudT::Ptr input, float radius, vector<PointCloudT::Ptr>& output);
void ptsToLine(PointCloudT::Ptr input, Eigen::VectorXf& paras, EdgeLine& output);
void findLinkedLines(vector<EdgeLine>& edgeLines);

void calculateMeanStandardDev(vector<float> v, float& mean, float& stdev);

void calculateNormals(PointCloudT::Ptr input, pcl::PointCloud <pcl::Normal>::Ptr &normals_all, int KSearch);
void regionGrow(PointCloudT::Ptr input, int NumberOfNeighbours, int SmoothnessThreshold, int CurvatureThreshold,
                int MinSizeOfCluster, int KSearch, vector<PointCloudT::Ptr>& outputClusters);

void removePtsAroundLine(PointCloudT::Ptr input, PointCloudT::Ptr output, vector<EdgeLine>& lines, float dist);

void regionGrow2D(PointCloudT::Ptr input, vector<PointCloudT::Ptr>& output);
void calculateNormal2D(PointCloudT::Ptr input, PointCloud<Normal>::Ptr cloud_normals);

void BeamRANSAC(PointCloudT::Ptr input, float high, vector< pair < vector<PointT>, float >> &planes, vector<pair<EdgeLine, float>> &lines);
bool findEdgeForPlane(PointCloudT::Ptr input, float maxRoomZ, PointCloudT::Ptr roofPart, vector<PointT>& edges);

void exportToDxf(string outputPath, const vector<EdgeLine>& wallLines, const vector< pair < vector<PointT>, float >> &beamPlanes,
                 const vector<pair<EdgeLine, float>> &beamLines, float minZ, float maxZ);

int checkNumofPointofLineinZ(PointT p, PointT q, PointCloudT::Ptr input, float minZ, float maxZ, PointCloudT::Ptr debug_cube);

void minimumCut(const vector<EdgeLine>& lines, PointCloudT::Ptr input);
void get_arrangement_face (CGAL::Arrangement_2<Traits_2>& arr, vector<CGALFace> &faces);
void print_ccb (Arrangement_2::Ccb_halfedge_const_circulator circ);

vector<Point_2> print_point_location
        (const typename Arrangement_2::Point_2& q,
         typename CGAL::Arr_point_location_result<Arrangement_2>::Type obj);

vector<Point_2> locate_point(const PointLocation& pl,
                  const PointLocation::Arrangement_2::Point_2& q);
#endif //TEST_PCL_EXTRACTWALL_H
