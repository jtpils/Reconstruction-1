#ifndef SIMPLEVIEW_H
#define SIMPLEVIEW_H

#include <iostream>
#include <pcl/visualization/pcl_visualizer.h>
#include <Plane.h>
#include <Reconstruction.h>
#include <main.h>
using namespace std;
typedef pcl::PointXYZRGBNormal PointT;
typedef pcl::PointCloud<PointT> PointCloudT;

void simpleView(const string &title, const pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr &cloud);

void simpleView(const string &title, const pcl::PointCloud<pcl::PointXYZ>::Ptr &cloud);

void simpleView(const string& title, vector<Plane> &planes);

void simpleView(const string& title, Reconstruction &re);

void simpleView(const string& title, vector<PointCloudT::Ptr> &clusters);

void simpleView(const string& title, vector<EdgeLine> &lines);

#endif