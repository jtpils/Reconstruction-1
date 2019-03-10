#include <iostream>
#include <vector>
#include <DxfExporter.h>
#include <fstream>
#include <sstream>

using namespace std;


KKRecons::DxfExporter::DxfExporter() {
}

void KKRecons::DxfExporter::insert(DxfFace face) {
    this->faces.push_back(face);
}

void KKRecons::DxfExporter::exportDXF(string path) {
    std::ofstream output(path);
    addHeaderPart(output);
    int numCompleteFace = 0;
    for (int i = 0; i < this->faces.size(); ++i) {
        DxfFace& face = this->faces[i];
        if (face.vacants.size() == 0){
            numCompleteFace++;
            addRectFace(output,face.a,face.b,face.c,face.d);
        }
    }
    for (int j = 0; j < numCompleteFace; ++j) {
        output << "0" << "\n";
        output << "VERTEX" << "\n";
        output << "5" << "\n";
        output << "AcDbEntity" << "\n";
        output << "100" << "\n";
        output << "AcDbFaceRecord" << "\n";
        output << "70" << "\n"; output << "128" << "\n";
        output << "71" << "\n"; output << 4*j+1 << "\n";
        output << "72" << "\n"; output << 4*j+2 << "\n";
        output << "73" << "\n"; output << 4*j+3 << "\n";
        output << "74" << "\n"; output << 4*j+4 << "\n";
    }

    vector<vector<int>> hides;
    for (int i = 0; i < this->faces.size(); ++i){
        DxfFace& face = this->faces[i];
        if(face.vacants.size() == 0) continue;
        DxfFace subFace = face.vacants[0];
        // the algorithm is not complete. the closest point might be not a.
        assert(pcl::geometry::distance(face.a,subFace.a) < pcl::geometry::distance(face.a,subFace.b));
        if (pcl::geometry::distance(face.a,subFace.a) < pcl::geometry::distance(face.a,subFace.b)) {
            addTriangleFace(output,face.a,subFace.a,subFace.b);
            hides.push_back(vector<int>{1,0,1});
            addTriangleFace(output,face.a,subFace.b,face.b);
            hides.push_back(vector<int>{1,1,0});
            addTriangleFace(output,face.b,subFace.b,subFace.c);
            hides.push_back(vector<int>{1,0,1});
            addTriangleFace(output,face.b,subFace.c,face.c);
            hides.push_back(vector<int>{1,1,0});
            addTriangleFace(output,face.c,subFace.c,subFace.d);
            hides.push_back(vector<int>{1,0,1});
            addTriangleFace(output,face.c,subFace.d,face.d);
            hides.push_back(vector<int>{1,1,0});
            addTriangleFace(output,face.d,subFace.d,subFace.a);
            hides.push_back(vector<int>{1,0,1});
            addTriangleFace(output,face.d,subFace.a,face.a);
            hides.push_back(vector<int>{1,1,0});
        }
    }
    int index = numCompleteFace*4;
    for (int k = 0; k < hides.size(); ++k) {
        output << "0" << "\n";
        output << "VERTEX" << "\n";
        output << "5" << "\n";
        output << "AcDbEntity" << "\n";
        output << "100" << "\n";
        output << "AcDbFaceRecord" << "\n";
        output << "70" << "\n"; output << "128" << "\n";
        output << "71" << "\n"; output << ( hides[k][0] == 1 ? -(k*3+1+index) : (k*3+1+index) )<< "\n";
        output << "72" << "\n"; output << ( hides[k][1] == 1 ? -(k*3+2+index) : (k*3+2+index) )<< "\n";
        output << "73" << "\n"; output << ( hides[k][2] == 1 ? -(k*3+3+index) : (k*3+3+index) )<< "\n";
    }
    addEndPart(output);
    output.close();
}

void KKRecons::DxfExporter::addHeaderPart(ofstream& output) {
    output << "999\n";
    output << "This DXF file is exported from KKReconstruct\n";
    output << "0\n";
    output << "SECTION\n";
    output << "2\n";
    output << "HEADER\n";
    output << "9\n";
    output << "$INSUNITS\n";
    output << "70\n";
    output << "6\n";
    output << "0\n";
    output << "SECTION\n";
    output << "2\n";
    output << "ENTITIES\n";
    output << "0\n";
    output << "POLYLINE\n";
    output << "100\n";
    output << "AcDbEntity\n";
    output << "100\n";
    output << "AcDbPolyFaceMesh\n";
    output << "66\n";
    output << "1\n";
    output << "10\n";
    output << "0.0\n";
    output << "20\n";
    output << "0.0\n";
    output << "30\n";
    output << "0.0\n";
    output << "70\n";
    output << "64\n";
    output << "71\n";
    output << "8\n";
    output << "72\n";
    output << "8\n";
}

void KKRecons::DxfExporter::addEndPart(ofstream &output) {
    output << "0\n";
    output << "SEQEND\n";
    output << "0\n";
    output << "ENDSEC\n";
    output << "0\n";
    output << "EOF\n";
}

void KKRecons::DxfExporter::addRectFace(ofstream &output, Point a, Point b, Point c, Point d) {
    output << "0\n";
    output << "VERTEX\n";
    output << "100\n";
    output << "AcDbEntity\n";
    output << "100\n";
    output << "AcDbVertex\n";
    output << "100\n";
    output << "AcDbPolyFaceMeshVertex\n";
    output << "10\n";
    output << a.x << "\n";
    output << "20\n";
    output << a.y << "\n";
    output << "30\n";
    output << a.z << "\n";
    output << "70\n";
    output << "192\n";
    output << "0\n";
    output << "VERTEX\n";
    output << "100\n";
    output << "AcDbEntity\n";
    output << "100\n";
    output << "AcDbVertex\n";
    output << "100\n";
    output << "AcDbPolyFaceMeshVertex\n";
    output << "10\n";
    output << b.x << "\n";
    output << "20\n";
    output << b.y << "\n";
    output << "30\n";
    output << b.z << "\n";
    output << "70\n";
    output << "192\n";
    output << "0\n";
    output << "VERTEX\n";
    output << "100\n";
    output << "AcDbEntity\n";
    output << "100\n";
    output << "AcDbVertex\n";
    output << "100\n";
    output << "AcDbPolyFaceMeshVertex\n";
    output << "10\n";
    output << c.x << "\n";
    output << "20\n";
    output << c.y << "\n";
    output << "30\n";
    output << c.z << "\n";
    output << "70\n";
    output << "192\n";
    output << "0\n";
    output << "VERTEX\n";
    output << "100\n";
    output << "AcDbEntity\n";
    output << "100\n";
    output << "AcDbVertex\n";
    output << "100\n";
    output << "AcDbPolyFaceMeshVertex\n";
    output << "10\n";
    output << d.x << "\n";
    output << "20\n";
    output << d.y << "\n";
    output << "30\n";
    output << d.z << "\n";
    output << "70\n";
    output << "192\n";
}

void KKRecons::DxfExporter::addTriangleFace(ofstream &output, Point a, Point b, Point c) {
    output << "0\n";
    output << "VERTEX\n";
    output << "100\n";
    output << "AcDbEntity\n";
    output << "100\n";
    output << "AcDbVertex\n";
    output << "100\n";
    output << "AcDbPolyFaceMeshVertex\n";
    output << "10\n";
    output << a.x << "\n";
    output << "20\n";
    output << a.y << "\n";
    output << "30\n";
    output << a.z << "\n";
    output << "70\n";
    output << "192\n";
    output << "0\n";
    output << "VERTEX\n";
    output << "100\n";
    output << "AcDbEntity\n";
    output << "100\n";
    output << "AcDbVertex\n";
    output << "100\n";
    output << "AcDbPolyFaceMeshVertex\n";
    output << "10\n";
    output << b.x << "\n";
    output << "20\n";
    output << b.y << "\n";
    output << "30\n";
    output << b.z << "\n";
    output << "70\n";
    output << "192\n";
    output << "0\n";
    output << "VERTEX\n";
    output << "100\n";
    output << "AcDbEntity\n";
    output << "100\n";
    output << "AcDbVertex\n";
    output << "100\n";
    output << "AcDbPolyFaceMeshVertex\n";
    output << "10\n";
    output << c.x << "\n";
    output << "20\n";
    output << c.y << "\n";
    output << "30\n";
    output << c.z << "\n";
    output << "70\n";
    output << "192\n";
}

string KKRecons::DxfExporter::getPlanesPart() {

}