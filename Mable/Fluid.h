//
//  Flip.h
//  Mable
//
//  Created by 須之内俊樹 on 2022/11/14.
//

#ifndef Fluid_h
#define Fluid_h

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
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <iostream>
#include "Array3d.hpp"
#define Nx 16
#define Ny 16
#define Nz 16//グリッドの数
#define g0 9.8
using ScalarType = double;
using IndexType = int64_t;
using Triplet = Eigen::Triplet<ScalarType,IndexType>;
using SparseMatrix = Eigen::SparseMatrix<ScalarType>;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor, int64_t> SpMat;
struct Fluid{
    double dx;//セルの大きさ
    double dt;//時間の刻み幅
    double rho;
    myArray3d u = myArray3d(Nx+1,Ny,Nz,0);//水平
    myArray3d old_u = myArray3d(Nx+1,Ny,Nz,0);//水平
    myArray3d umi = myArray3d(Nx+1,Ny,Nz,0);//Gridの重さ
    myArray3d ufi = myArray3d(Nx+1,Ny,Nz,0);//Gridに加わる外力
    
    myArray3d v = myArray3d(Nx,Ny+1,Nz,0);//鉛直
    myArray3d old_v = myArray3d(Nx,Ny+1,Nz,0);//水平
    myArray3d vmi = myArray3d(Nx,Ny+1,Nz,0);//Gridの重さ
    myArray3d vfi = myArray3d(Nx,Ny+1,Nz,0);//Gridに加わる外力
    
    myArray3d w = myArray3d(Nx,Ny,Nz+1,0);
    myArray3d old_w = myArray3d(Nx,Ny,Nz+1,0);
    myArray3d wmi = myArray3d(Nx,Ny,Nz+1,0);//Gridの重さ
    myArray3d wfi = myArray3d(Nx,Ny,Nz+1,0);//Gridに加わる外力
    
    myArray3d p = myArray3d(Nx,Ny,Nz,0);//圧力
    std::vector<double> weights;//粒子の重み
    double L;
    Eigen::Vector3d f0 = {0.0,-g0,0.0};
    Fluid(double x,double t,double density){
        dx = x;
        dt = t;
        rho = density;
        L = dx*Nx;
        vfi.reset(f0.y());
    }
    std::vector<int>DirichletBoundaryCondition(int i,int j,int k,std::unordered_map<std::vector<int>,std::vector<int>,ArrayHasher<3>>&map){
        std::vector<int>ret(6,1);
        if(i == Nx-1)ret[0] = 0;
        else if(map.find({i+1,j,k}) == map.end())ret[0] = 0;
        if(j == Ny-1)ret[1] = 0;
        else if(map.find({i,j+1,k}) == map.end())ret[1] = 0;

        if(i == 0)ret[2] = 0;
        else if(map.find({i-1,j,k}) == map.end())ret[2] = 0;
        if(j == 0)ret[3] = 0;
        else if(map.find({i,j-1,k}) == map.end())ret[3] = 0;

        if(k == 0)ret[4] = 0;
        else if(map.find({i,j,k-1}) == map.end())ret[4] = 0;
        if(k == Nz-1)ret[5] = 0;
        else if(map.find({i,j,k+1}) == map.end())ret[5] = 0;
        return ret;
    }
    void project(std::unordered_map<std::vector<int>,std::vector<int>,ArrayHasher<3>>&map){
        SparseMatrix A(Nx*Ny*Nz,Nx*Ny*Nz),B(Nx*Ny*Nz,Nx*Ny*Nz);
        Eigen::VectorXd b = Eigen::VectorXd::Zero(Nx*Ny*Nz);
        Eigen::VectorXd px;
        std::set<int> DirichletKey;
        std::vector<std::vector<int>>keys;
        //Tripletの計算
        std::vector<Triplet> triplets;
        for(int i=0;i<Nx;i++){
            for(int j=0;j<Ny;j++){
                for(int k=0;k<Nz;k++){
                    std::vector<int>key = {i,j,k};
                    if(map.find(key) == map.end()){
                        //前処理でAが変更されてしまうので，境界条件として別で無理矢理設定する．
                        triplets.emplace_back(i+j*Nx+k*Nx*Ny,i+j*Nx+k*Nx*Ny,1);
                        keys.push_back(key);
                        DirichletKey.insert(i+j*Nx+k*Nx*Ny);
                        continue;
                    }
                    double scale = dt/(rho*dx*dx);
                    //std::cout << i << "," << j << std::endl;
                    double D[6] = {1.0,1.0,-1.0,-1.0,-1.0,1.0};//周囲6方向に向かって働く、圧力の向き
                    //double F[4] = {(double)(i<Nx-1),(double)(j<Ny-1),(double)(i>0),(double)(j>0)};//境界条件。壁なら0,流体なら1
                    
                    std::vector<int> F = {i<Nx-1,j<Ny-1,i>0,j>0,k>0,k<Nz-1};
                    double U[6] = {
                        u.value[i+1][j][k],
                        v.value[i][j+1][k],
                        u.value[i][j][k],
                        v.value[i][j][k],
                        w.value[i][j][k],
                        w.value[i][j][k+1]};
                    double sumP = 0;
                    for(int n=0;n<6;n++){
                        sumP += -F[n]*scale;
                        //sumP += scale;
                        b(i+j*Nx+k*Nx*Ny) += D[n]*F[n]*U[n]/(dx);
                    }
                    F = DirichletBoundaryCondition(i,j,k,map);
//                    for(int n=0;n<6;n++){
//                        std::cout << F_pri[n] << "," << F[n] << std::endl;
//                    }
                    triplets.emplace_back(i+j*Nx+k*Nx*Ny,i+j*Nx+k*Nx*Ny, sumP);
                    if(F[0])triplets.emplace_back(i+j*Nx+k*Nx*Ny,i+1+j*Nx+k*Nx*Ny, F[0]*scale);
                    if(F[1])triplets.emplace_back(i+j*Nx+k*Nx*Ny,i+(j+1)*Nx+k*Nx*Ny, F[1]*scale);
                    if(F[2])triplets.emplace_back(i+j*Nx+k*Nx*Ny,i-1+j*Nx+k*Nx*Ny, F[2]*scale);
                    if(F[3])triplets.emplace_back(i+j*Nx+k*Nx*Ny,i+(j-1)*Nx+k*Nx*Ny, F[3]*scale);
                    if(F[4])triplets.emplace_back(i+j*Nx+k*Nx*Ny,i+j*Nx+(k-1)*Nx*Ny, F[4]*scale);
                    if(F[5])triplets.emplace_back(i+j*Nx+k*Nx*Ny,i+j*Nx+(k+1)*Nx*Ny, F[5]*scale);
                }
            }
        }
        A.setFromTriplets(triplets.begin(), triplets.end());
        Eigen::ConjugateGradient<SparseMatrix> solver;
        
        //Eigen::BiCGSTAB<SparseMatrix> solver;
        
        for(int i=0;i<A.outerSize();++i){
            for(SparseMatrix::InnerIterator it(A,i);it;++it){
                if(it.row() == *DirichletKey.begin()){
                    it.valueRef() = 1;
                    DirichletKey.erase(DirichletKey.begin());
                }
            }
        }
        solver.compute(A);
        px = solver.solve(b);
        for(int i=0;i<Nx;i++){
            for(int j=0;j<Ny;j++){
                for(int k=0;k<Nz;k++){
                    p.value[i][j][k] = px(i+j*Nx+k*Nx*Ny);
                }
            }
        }
        for(auto x:keys){
            if(px(x[0] + x[1]*Ny + x[2]*Nz*Nz) > 1.0e-4)std::cout << x[0] + x[1]*Ny + x[2]*Nz*Nz << "," << px(x[0] + x[1]*Ny + x[2]*Nz*Nz) <<std::endl;
            //else std::cout << "success:" << x[0] + x[1]*Ny + x[2]*Nz*Nz << std::endl;
        }
        
        for(int i=1; i<Nx;i++){
            for(int j=0;j<Ny;j++){
                for(int k=0;k<Nz;k++)u.value[i][j][k] = u.value[i][j][k] - dt/rho * (p.value[i][j][k]-p.value[i-1][j][k])/dx;
            }
        }
        for(int i=0;i<Nx;i++){
            for(int j=1;j<Ny;j++){
                for(int k=0;k<Nz;k++)v.value[i][j][k] = v.value[i][j][k] - dt/rho * (p.value[i][j][k]-p.value[i][j-1][k])/dx;
            }
        }
        for(int i=0;i<Nx;i++){
            for(int j=0;j<Ny;j++){
                for(int k=1;k<Nz;k++)w.value[i][j][k] = w.value[i][j][k] - dt/rho * (p.value[i][j][k]-p.value[i][j][k-1])/dx;
            }
        }
    }
};
#endif /* Fluid_h */
