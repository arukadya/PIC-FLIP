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
#include <iostream>
#define Nx 50
#define Ny 50 //グリッドの数
#define g0 9.8
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
        double scale = dt/(rho*dx*dx);//左辺の係数部分
        double eps = 1.0e-4;//ガウスザイデル法の精度
        double err;//ガウスザイデル法の残差
//        std::cout << "inputP" << std::endl;
//        print_pressure();
        do{
            err = 0.0;
        //圧力場を求める
            for(int j=0;j<Ny;j++)for(int i=0;i<Nx;i++){
                std::vector<int>key = {i,j};
                if(map.find(key) == map.end()){
                    p[i][j] = 0;
                    continue;
                }
                //std::cout << i << "," << j << std::endl;
                double D[4] = {1.0,1.0,-1.0,-1.0};//周囲4方向に向かって働く、圧力の向き
                //double F[4] = {(double)(i<Nx-1),(double)(j<Ny-1),(double)(i>0),(double)(j>0)};//境界条件。壁なら0,流体なら1
                std::vector<int> F = DirichletBoundaryCondition(i,j,map);
                double P[4] = {//pn。周囲4つのセルの圧力値
                    (F[0] ? p[i+1][j] : 0.0),
                    (F[1] ? p[i][j+1] : 0.0),
                    (F[2] ? p[i-1][j] : 0.0),
                    (F[3] ? p[i][j-1] : 0.0)};
                
                F = {i<Nx-1,j<Ny-1,i>0,j>0};
                //std::cout << P[0] << " " << P[1] << " " << P[2] << " " << P[3] << std::endl;
                //std::cout << i << "," << j << std::endl;
                //std::cout << i << "," << j << std::endl;
                double U[4] = {u[i+1][j],v[i][j+1],u[i][j],v[i][j]};
                //セルの圧力値を求める
                double det = 0.0;//左辺のdeterminant
                int sumF = 0;
                double sum_L = 0.0;
                double sum_R = 0.0;
                for(int n=0;n<4;n++){
                    sumF += F[n];
                    det += F[n]*scale;
                    //sum_L += F[n]*P[n]*scale;
                    sum_L += P[n]*scale;
                    sum_R += F[n]*D[n]*U[n]/dx;
                }
                //std::cout << "det:" << det << " sum_L:" << sum_L << " sum_R:" << sum_R << std::endl;
                err = fmax(err,fabs(det*p[i][j]-sum_L+sum_R));
                if(sumF == 0)p[i][j] = 0;
                else p[i][j] = (sum_L-sum_R)/det;
                //std::cout << p[i][j] <<std::endl;
            }
        }while(eps<err);//反復回数は初期値と収束速度？に依存。定数回なら全体で計算量はO(n)
        //新しい流速を求める。u(t+1) = u* - (dt/rho)grad(p)
        //std::cout << "finGauss" << std::endl;
//        print_velocity();
        for(int i=1; i<Nx;i++)for(int j=0;j<Ny;j++){
            u[i][j] = u[i][j] - dt/rho * (p[i][j]-p[i-1][j])/dx;
            //delta_u[i][j] = -dt/rho * (p[i][j]-p[i-1][j])/dx;
//---------多分質量０のグリッドの速さは０である------------------------------------------------------------
//            if(umi[i][j] < eps){
//                //u[i][j] =u[i][j] = u[i][j] - dt/rho * (p[i][j]-p[i-1][j])/dx;
//                u[i][j] = 0;
//                delta_u[i][j] = -u[i][j];
//            }
//            else{
//               u[i][j] = u[i][j] - dt/rho * (p[i][j]-p[i-1][j])/dx + dt*ufi[i][j]/umi[i][j];
//                delta_u[i][j] = -dt/rho * (p[i][j]-p[i-1][j])/dx + dt*ufi[i][j]/umi[i][j];
//                u[i][j] = u[i][j] - dt/rho * (p[i][j]-p[i-1][j])/dx;
//                delta_u[i][j] = -dt/rho * (p[i][j]-p[i-1][j])/dx;
//            }
//--------------------------------------------------------------------------------------------------
        }
//        for(int j=0;j<Ny;j++){
//            u[0][j] = -u[1][j];
//        }
        for(int i=0;i<Nx;i++)for(int j=1;j<Ny;j++){
            v[i][j] = v[i][j] - dt/rho * (p[i][j]-p[i][j-1])/dx;
            //delta_v[i][j] = -dt/rho * (p[i][j]-p[i][j-1])/dx;
//            if(vmi[i][j] < eps){
//               v[i][j] = v[i][j] - dt/rho * (p[i][j]-p[i][j-1])/dx;
//                v[i][j] = 0;
//                delta_v[i][j] = -v[i][j];
//            }
//            else {
//                v[i][j] = v[i][j] - dt/rho * (p[i][j]-p[i][j-1])/dx + dt*vfi[i][j]/vmi[i][j];
//                delta_v[i][j] = -dt/rho * (p[i][j]-p[i][j-1])/dx + dt*vfi[i][j]/vmi[i][j];
//                v[i][j] = v[i][j] - dt/rho * (p[i][j]-p[i][j-1])/dx;
//                delta_v[i][j] = -dt/rho * (p[i][j]-p[i][j-1])/dx;
//            }
        }
//        for(int i=0;i<Nx;i++){
//            v[i][0] = -v[i][1];
//        }
    }
    void print_pressure(){
        for(int i=0;i<Nx;i++){
            for(int j=0;j<Ny;j++){
            std::cout << p[i][j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    void print_velocity(){
        for(int i=0;i<Nx;i++){
            for(int j=0;j<Ny;j++){
            std::cout << v[i][j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    void initPressure(){
        for(unsigned int i=0;i<Nx;i++)for(unsigned int j=0;j<Ny;j++){
            if(i > Nx/2)p[i][j] = 0;
            //else p[i][j] = (double)(Nx-1-i)/Ny;
            else p[i][j] = 0;
        }
    }
    void initForce(){
        Eigen::Vector2d f0 = {0.0,g0*dx};
        Eigen::Vector2d f1 = {0.0,0.0};
        Eigen::Vector2d f2 = {0.0,g0*dx};
        for(unsigned int i=0;i<Nx+1;i++)for(unsigned int j=0;j<Ny;j++){
            //std::cout << f0.x() << "," << f0.y() << std::endl;
            if(i == 0 || i == Nx || j == 0 || j == Ny-1)ufi[i][j] = f2.x();
            else ufi[i][j] = f0.x();
        }
        for(unsigned int i=0;i<Nx;i++)for(unsigned int j=0;j<Ny+1;j++){
            //std::cout << f0.x() << "," << f0.y() << std::endl;
            if(i == 0 || i == Nx || j == 0 || j == Ny-1)vfi[i][j] = f1.y();
            else vfi[i][j] = f0.y();
        }
    }
};
#endif /* Fluid_h */
