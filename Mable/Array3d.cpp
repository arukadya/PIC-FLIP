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
