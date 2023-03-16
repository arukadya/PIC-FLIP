//
//  marching_cubes.hpp
//  marchining_cubes
//
//  Created by 須之内俊樹 on 2023/03/04.
//

#ifndef marching_cubes_hpp
#define marching_cubes_hpp

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <string>
#include <math.h>
#include <sstream>
#include <iomanip>
#include <set>
#include <unordered_map>
#include <queue>
#include "Hasher.h"
#include "Eigen/Core"

using VertexID = std::tuple<int,int,int>;
using edge = std::pair<int,int>;
using worldEdge = std::pair<VertexID,VertexID>;
const long long INF=1LL<<60;
template <typename T>
struct ImplicitFunction{
    int nx,ny,nz;
    T dx,dy,dz;
    int num_cells;
    std::vector<std::vector<std::vector<T>>> values;
    ImplicitFunction(int num_cells,int nx,int ny,int nz, T dx,T dy,T dz,std::vector<T> &data){
        this->num_cells = num_cells;
        this->nx = nx;
        this->ny = ny;
        this->nz = nz;
        this->dx = dx;
        this->dy = dy;
        this->dz = dz;
        values.resize(nx);
        for(int i=0;i<nx;++i){
            values[i].resize(ny);
            for(int j=0;j<ny;++j){
                values[i][j].resize(nz);
                for(int k=0;k<nz;++k){
                    //values[i][j][k] = data[i + j*nx + k*nx*ny];
                    values[i][j][k] = data[k + j*nz + i*nz*ny];
                }
            }
        }
    }
    void print_values(){
        for(int i=0;i<nx;++i){
            for(int j=0;j<ny;++j){
                for(int k=0;k<nz;++k){
                    std::cout << values[i][j][k] << ",";
                }std::cout << std::endl;
            }std::cout << std::endl;
        }std::cout << std::endl;
    }
    T getValue(VertexID id){
        auto [x,y,z] = id;
        return values[x][y][z];
    }
};

void marching_cubes(std::vector<Eigen::Vector3d> &Vertices,std::vector<std::vector<int>> &Faces,Eigen::Vector3d origin,Eigen::Vector3d &dist,ImplicitFunction<double> &imp,double threshold);

void addTriangle(VertexID glidID, Eigen::Vector3d origin, Eigen::Vector3d dist, double threshold, ImplicitFunction<double> &imp,  std::vector<Eigen::Vector3d> &Vertices, std::vector<std::vector<int>> &Faces,std::unordered_map<std::vector<int>, int, ArrayHasher<2>> &map, int &vert_index);

void addVertex(int l_x0ID, int l_x1ID, VertexID glidID, Eigen::Vector3d origin, Eigen::Vector3d dist, double threshold, ImplicitFunction<double> &imp,  std::vector<Eigen::Vector3d> &Vertices, std::unordered_map<std::vector<int>, int, ArrayHasher<2>> &map, int &vert_index);

Eigen::Vector3d linearInterporation(double isoValue,Eigen::Vector3d x0,Eigen::Vector3d x1,double value0,double value1);

//変換する関数群
Eigen::Vector3d getWorldPosition(VertexID glidID, Eigen::Vector3d origin, Eigen::Vector3d dist);

int getVertexTable(VertexID v_id,ImplicitFunction<double> &implicitfunction,double threshold);

int transIndexL2W(VertexID glidID, int localID, int nx, int ny, int nz);

std::pair<int,int> edge2Vertex(int edgeID);

std::vector<int> getKey(VertexID w_id0,VertexID w_id1, int size);
//テスト用
void cubetest(VertexID id,ImplicitFunction<double> &implicitfunction);

#endif /* marching_cubes_hpp */
