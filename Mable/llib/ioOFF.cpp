//
//  ioOFF.cpp
//  Mable
//
//  Created by 須之内俊樹 on 2023/03/16.
//

#include "ioOFF.hpp"
void inputOFF(const char* InputFileName,std::vector<Eigen::Vector3d> &Vertices,std::vector<std::vector<int>> &Faces){
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
void outputOFF(std::string OutputFileName,std::vector<Eigen::Vector3d> &Vertices,std::vector<std::vector<int>> &Faces){
    std::ofstream writing_file;
    writing_file.open(OutputFileName, std::ios::out);
    writing_file << "OFF"<< std::endl;
    writing_file << Vertices.size() << " " << Faces.size() << " 0" <<std::endl;
    for(auto &x:Vertices){
        writing_file << x(0) << " " << x(1) << " " << x(2) << std::endl;
    }
    for(int i=0;i<Faces.size();i++){
        writing_file << "3 ";
        for(int j=0;j < 3;j++){
            writing_file << Faces[i][j];
            if(j != 2)writing_file << " ";
            else writing_file << std::endl;
        }
    }
    writing_file.close();
}
