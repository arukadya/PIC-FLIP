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
#define Nx 64
#define Ny 64 //グリッドの数
#define g0 9.8
using ScalarType = double;
using IndexType = int64_t;
using Triplet = Eigen::Triplet<ScalarType, IndexType>;
using SparseMatrix = Eigen::SparseMatrix<ScalarType>;
struct Fluid{
    double dx;//セルの大きさ
    double dt;//時間の刻み幅
    double rho;
    std::vector<std::vector<double>>u;//水平
    std::vector<std::vector<double>>old_u;//水平
    std::vector<std::vector<double>>v;//鉛直
    std::vector<std::vector<double>>old_v;//水平
    //std::vector<std::vector<double>>delta_u;//水平
    //std::vector<std::vector<double>>delta_v;//鉛直
    std::vector<std::vector<double>>p;//圧力
    std::vector<std::vector<double>>p_Eigen;//圧力
    std::vector<std::vector<double>>umi;//Gridの重さ
    std::vector<std::vector<double>>vmi;//Gridの重さ
    //std::vector<std::vector<Eigen::Vector2d>>fi;//Gridに加わる外力
    std::vector<std::vector<double>>ufi;//Gridに加わる外力
    std::vector<std::vector<double>>vfi;//Gridに加わる外力
    std::vector<double> weights;//粒子の重み
    double L;
    
    Fluid(double x,double t,double density,std::vector<std::vector<double>> &horizontal_v,std::vector<std::vector<double>> &vertical_v,std::vector<std::vector<double>> &pressure,std::vector<std::vector<double>>&gridUM,std::vector<std::vector<double>>&gridVM,std::vector<std::vector<double>>&gridUF,std::vector<std::vector<double>>&gridVF){
        dx = x;
        dt = t;
        rho = density;
        u = horizontal_v;//v[nx+1][ny]
        old_u = horizontal_v;
        //delta_u = horizontal_v;
        v = vertical_v;//v[nx][ny+1]
        old_v = vertical_v;
        //delta_v = vertical_v;
        p = pressure;//p[nx][ny]
        p_Eigen = pressure;
        umi = gridUM;
        vmi = gridVM;
        ufi = gridUF;
        vfi = gridVF;
        L = dx*Nx;
        initPressure();
        std::cout << "initpressure" << std::endl;
        initForce();
        std::cout << "initforce" << std::endl;
    }
    std::vector<int>DirichletBoundaryCondition(int i,int j,std::unordered_map<std::vector<int>,std::vector<int>,ArrayHasher<2>>&map){
        std::vector<int>ret(4,1);
        if(i == 0)ret[2] = 0;
        else if(map.find({i-1,j}) == map.end())ret[2] = 0;
        if(i == Nx-1)ret[0] = 0;
        else if(map.find({i+1,j}) == map.end())ret[0] = 0;
        if(j == 0)ret[3] = 0;
        else if(map.find({i,j-1}) == map.end())ret[3] = 0;
        if(j == Ny-1)ret[1] = 0;
        else if(map.find({i,j+1}) == map.end())ret[1] = 0;
        return ret;
    }
    void project(std::unordered_map<std::vector<int>,std::vector<int>,ArrayHasher<2>>&map){
        SparseMatrix A(Nx*Ny,Nx*Ny);
        Eigen::VectorXd b = Eigen::VectorXd::Zero(Nx*Ny);
        Eigen::VectorXd px;
        //Tripletの計算
        std::vector<Triplet> triplets;
        for(int i=0;i<Nx;i++){
            for(int j=0;j<Ny;j++){
                std::vector<int>key = {i,j};
                if(map.find(key) == map.end()){
                    triplets.emplace_back(i*Nx+j,i*Nx+j, 1);
                    //p[i][j] = 0;
                    continue;
                }
                double scale = dt/(rho*dx*dx);
                //std::cout << i << "," << j << std::endl;
                double D[4] = {1.0,1.0,-1.0,-1.0};//周囲4方向に向かって働く、圧力の向き
                //double F[4] = {(double)(i<Nx-1),(double)(j<Ny-1),(double)(i>0),(double)(j>0)};//境界条件。壁なら0,流体なら1
                std::vector<int> F = {i<Nx-1,j<Ny-1,i>0,j>0};
                double U[4] = {u[i+1][j],v[i][j+1],u[i][j],v[i][j]};
                double sumP = 0;
                for(int n=0;n<4;n++){
                    sumP += -F[n]*scale;
                    b(i*Nx+j) += D[n]*F[n]*U[n]/dx;
                }
                //F = DirichletBoundaryCondition(i,j,map);
                triplets.emplace_back(i*Nx+j,i*Nx+j, sumP);
                if(F[0])triplets.emplace_back(i*Nx+j,(i+1)*Nx+j, F[0]*scale);
                if(F[1])triplets.emplace_back(i*Nx+j,i*Nx+j+1, F[1]*scale);
                if(F[2])triplets.emplace_back(i*Nx+j,(i-1)*Nx+j, F[2]*scale);
                if(F[3])triplets.emplace_back(i*Nx+j,i*Nx+j-1, F[3]*scale);
            }
        }
        A.setFromTriplets(triplets.begin(), triplets.end());
//        Eigen::ConjugateGradient<SparseMatrix> solver;
        Eigen::BiCGSTAB<SparseMatrix> solver;
        solver.compute(A);
        px = solver.solve(b);
        for(int i=0;i<Nx;i++){
            for(int j=0;j<Ny;j++){
                //p_Eigen[i][j] = px(i*Nx+j);
                p[i][j] = px(i*Nx+j);
            }
        }
        for(int i=1; i<Nx;i++)for(int j=0;j<Ny;j++){
            u[i][j] = u[i][j] - dt/rho * (p[i][j]-p[i-1][j])/dx;
            //if(fabs(u[i][j] > 1.0e-10))std::cout << i << "," << j << ":"<< u[i][j] << std::endl;
            //delta_u[i][j] = -dt/rho * (p[i][j]-p[i-1][j])/dx;
        }
        for(int i=0;i<Nx;i++)for(int j=1;j<Ny;j++){
            v[i][j] = v[i][j] - dt/rho * (p[i][j]-p[i][j-1])/dx;
        }
    }
    void print_pressure(){
        for(int i=0;i<Nx;i++){
            for(int j=0;j<Ny;j++){
                std::cout << i <<","<< j << " " <<p[i][j] << "," << p_Eigen[i][j] << std::endl;
            }
        }
    }
    void print_velocity(){
        //for(int i=0;i<Nx;i++){
            for(int j=0;j<Ny;j++){
            std::cout << u[1][j] << " ";
            }
            std::cout << std::endl;
        for(int j=0;j<Ny;j++){
        std::cout << u[Nx-1][j] << " ";
        }
        std::cout << std::endl;
        //}
        //std::cout << std::endl;
    }
    void initPressure(){
        for(unsigned int i=0;i<Nx;i++)for(unsigned int j=0;j<Ny;j++){
            if(i > Nx/2)p[i][j] = 0;
            //else p[i][j] = (double)(Nx-1-i)/Ny;
            else p[i][j] = 0;
        }
    }
    void initForce(){
        Eigen::Vector2d f0 = {0.0,-g0};
        Eigen::Vector2d f1 = {0.0,0.0};
        Eigen::Vector2d f2 = {0.0,g0};
        for(unsigned int i=0;i<Nx+1;i++)for(unsigned int j=0;j<Ny;j++){
            //std::cout << f0.x() << "," << f0.y() << std::endl;
//            if(i == 0 || i == Nx || j == 0 || j == Ny-1)ufi[i][j] = f2.x();
//            else ufi[i][j] = f0.x();
            ufi[i][j] = f0.x();
        }
        for(unsigned int i=0;i<Nx;i++)for(unsigned int j=0;j<Ny+1;j++){
            //std::cout << f0.x() << "," << f0.y() << std::endl;
//            if(i == 0 || i == Nx || j == 0 || j == Ny-1)vfi[i][j] = f1.y();
//            else vfi[i][j] = f0.y();
            vfi[i][j] = f0.y();
        }
    }
};
#endif /* Fluid_h */
