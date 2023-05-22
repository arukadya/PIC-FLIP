//
//  Array3d.cpp
//  Mable
//
//  Created by 須之内俊樹 on 2023/02/26.
//

#include "Array3d.hpp"
    
myArray3d::myArray3d(int size_x,int size_y,int size_z){
    nx = size_x;
    ny = size_y;
    nz = size_z;
    size = nx*ny*nz;
    value.resize(nx);
    for(int i=0;i<nx;i++){
        value[i].resize(ny);
        for(int j=0;j<ny;j++){
            value[i][j].resize(nz);
        }
    }
}
myArray3d::myArray3d(int size_x,int size_y,int size_z,double val){
    nx = size_x;
    ny = size_y;
    nz = size_z;
    size = nx*ny*nz;
    value.resize(nx);
    for(int i=0;i<nx;i++){
        value[i].resize(ny);
        for(int j=0;j<ny;j++){
            value[i][j].resize(nz);
            for(int k=0;k<nz;k++)value[i][j][k] = val;
        }
    }
}
void myArray3d::reset(double val){
    for(int i=0;i<nx;i++){
        for(int j=0;j<ny;j++){
            for(int k=0;k<nz;k++)value[i][j][k] = val;
        }
    }
}
void myArray3d::print(){
    for(int i=0;i<nx;i++){
        for(int j=0;j<ny;j++){
            for(int k=0;k<nz;k++)std::cout << value[i][j][k] << ",";
        }std::cout << std::endl;
    }std::cout << std::endl;
}

std::vector<double> myArray3d::convert2Vector(){
    std::vector<double> data(nx*ny*nz);
    for(int i=0;i<nx;++i){
        for(int j=0;j<ny;++j){
            for(int k=0;k<nz;++k){
                data[i + j*nx + k*nx*ny] = value[i][j][k];
            }
        }
    }
    return data;
}
myMap::myMap(){}
myMap::myMap(int size_x,int size_y,int size_z){
    nx = size_x;
    ny = size_y;
    nz = size_z;
    std::vector<int> zero = {};
    value.resize(size_x);
    //size.resize(size_x);
    for(int i=0;i<size_x;++i){
        value[i].resize(size_y);
        //size[i].resize(size_y);
        for(int j=0;j<size_y;j++){
            value[i][j].resize(size_z);
            //size[i][j].resize(size_z);
            for(int k=0;k<size_z;k++){
                value[i][j][k] = zero;
                //size[i][j][k] = 0;
            }
        }
    }
}
void myMap::reset(){
    for(int i=0;i<nx;i++){
        for(int j=0;j<ny;j++){
            for(int k=0;k<nz;k++)value[i][j][k] = std::vector<int>();
        }
    }
}
void myMap::print(){
    for(int i=0;i<nx;i++){
        for(int j=0;j<ny;j++){
            for(int k=0;k<nz;k++){
                if(!value[i][j][k].empty())std::cout << value[i][j][k].size() << " ";
                else std::cout << 0 << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
bool myMap::contains(std::vector<int> &key){
    //std::cout << key[0] << "," << key[1] << "," << key[2] << "," << value[key[0]][key[1]][key[2]].size() << std::endl;
    //if(size[key[0]][key[1]][key[2]] == 0)return false;
    if(value[key[0]][key[1]][key[2]].empty())return false;
    else return true;
}
void myMap::push_back_particles(std::vector<int> &key,int i){
    //std::cout << "p" << key[0] << "," << key[1] << "," << key[2] << "," << i << "," << value[key[0]][key[1]][key[2]].size() << std::endl;
    //++size[key[0]][key[1]][key[2]];
    value[key[0]][key[1]][key[2]].push_back(i);
}
std::vector<int> myMap::at(std::vector<int> &key){
    //std::cout << key[0] << "," << key[1] << "," << key[2] << "," << value[key[0]][key[1]][key[2]].size() << std::endl;
    return value[key[0]][key[1]][key[2]];
}
