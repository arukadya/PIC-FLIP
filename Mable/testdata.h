//
//  testdata.h
//  Mable
//
//  Created by 須之内俊樹 on 2022/11/13.
//

#ifndef testdata_h
#define testdata_h
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <math.h>
#include <sstream>
#include <iomanip>
#include <set>
#include <map>
#include <unordered_map>
#include <queue>
#include <iostream>

#endif /* testdata_h */

void twoWave(int Nx,int Ny,double dx,double** u,double** v,double** p){
    double v0 = 1.0;//初期値流速
    double p0 = 1.0;//初期値圧力
    
    //初期化
    for(unsigned int i=0;i<Ny;i++)for(unsigned int j=0;j<Nx+1;j++)u[i][j] = 0.0;
    for(unsigned int i=0;i<Ny+1;i++)for(unsigned int j=0;j<Nx;j++)v[i][j] = 0.0;
    for(unsigned int i=0;i<Ny;i++)for(unsigned int j=0;j<Nx;j++)p[i][j] = 0.0;
    //2本の縦線の計算
    for(unsigned int i=0;i<Ny+1;i++){
        v[Ny/3][i] = v0;
        v[Ny/3*2][i] = v0;
        p[Ny/3][i] = p0;
        p[Ny/3*2][i] = p0;
    }/*
    for(unsigned int i=0;i<Ny;i++)for(unsigned int j=0;j<Nx+1;j++){
        std::cout << v[i][j] << ",";
    }
      */
          
}
