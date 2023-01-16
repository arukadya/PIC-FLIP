//
//  Array3d.h
//  Mable
//
//  Created by 須之内俊樹 on 2023/01/16.
//

#ifndef Array3d_h
#define Array3d_h
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <math.h>
#include <sstream>
#include <iomanip>
#include <queue>
#include <Eigen/Core>
#include <iostream>
struct myArray3d{
    std::vector<std::vector<std::vector<double>>>value;
    int nx,ny,nz;
    myArray3d();
    myArray3d(int size_x,int size_y,int size_z){
        nx = size_x;
        ny = size_y;
        nz = size_z;
        value.resize(nx);
        for(int i=0;i<nx;i++){
            value[i].resize(ny);
            for(int j=0;j<ny;j++){
                value[i][j].resize(nz);
            }
        }
    }
    myArray3d(int size_x,int size_y,int size_z,double val){
        nx = size_x;
        ny = size_y;
        nz = size_z;
        value.resize(nx);
        for(int i=0;i<nx;i++){
            value[i].resize(ny);
            for(int j=0;j<ny;j++){
                value[i][j].resize(nz);
                for(int k=0;k<nz;k++)value[i][j][k] = val;
            }
        }
    }
    void reset(double val){
        for(int i=0;i<nx;i++){
            for(int j=0;j<ny;j++){
                for(int k=0;k<nz;k++)value[i][j][k] = val;
            }
        }
    }
};
#endif /* Array3d_h */
