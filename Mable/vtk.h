//
//  vtk.h
//  Mable
//
//  Created by 須之内俊樹 on 2022/11/21.
//

#ifndef vtk_h
#define vtk_h

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <Eigen/Core>
#endif
void inputVTK(const char* InputFileName,std::vector<Eigen::Vector3d> &Vertices,std::vector<std::vector<int>> &Faces){
    std::cout << "InputFileName:" << InputFileName << std::endl;
    //std::cout << "OutputFileName:" << OutputFileName << std::endl;
    FILE *ifp = fopen(InputFileName,"r");
    int num_vertices, num_faces,dummy;
    fscanf(ifp, "OFF %d %d %d", &num_vertices, &num_faces, &dummy);
    for(int i=0;i<num_vertices;i++){
        double x,y,z;//点の入力
        fscanf(ifp, "%lf %lf %lf", &x, &y, &z);
        Eigen::Vector3d v(x,y,z);
        Vertices.push_back(v);
    }
    for(int i=0;i<num_faces;i++){//面の入力
        int num_sides, v0, v1, v2;
        fscanf(ifp, "%d %d %d %d", &num_sides, &v0, &v1, &v2);
        Faces.push_back({v0,v1,v2});
    }
    fclose(ifp);
}
void outputVTK(const char* OutputFileName,std::vector<Eigen::Vector3d> &Vertices){
    FILE *ofp = fopen(OutputFileName,"w");
    fprintf(ofp,"# vtk DataFile Version 2.0\n");
    fprintf(ofp,"Title Data\n");
    fprintf(ofp,"ASCII\n");
    fprintf(ofp,"DATASET UNSTRUCTURED_GRID\n");
    fprintf(ofp,"POINTS %d float\n",(int)Vertices.size());
    for(auto &x:Vertices){
        fprintf(ofp,"%lf %lf %lf\n",x(0),x(1),x(2));
    }
    fprintf(ofp,"\n");
    fprintf(ofp,"CELL_TYPES %d\n",(int)Vertices.size());
    for(int i=0;i<Vertices.size();i++)fprintf(ofp,"1\n");
    fprintf(ofp,"\n");
    fprintf(ofp,"POINT_DATA %d\n",(int)Vertices.size());
    fprintf(ofp,"SCALARS radius float\n");
    fprintf(ofp,"LOOKUP_TABLE default\n");
    for(int i=0;i<Vertices.size();i++)fprintf(ofp,"0.5\n");
    fclose(ofp);
}
