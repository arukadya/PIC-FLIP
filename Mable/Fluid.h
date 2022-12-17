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
    std::vector<std::vector<double>>v;//鉛直
    std::vector<std::vector<double>>p;//圧力
    std::vector<std::vector<double>>mi;//Gridの重さ
    std::vector<std::vector<Eigen::Vector2d>>fi;//Gridに加わる外力
    std::vector<double> weights;//粒子の重み
    double L;
    
    Fluid(double x,double t,double density,std::vector<std::vector<double>> &horizontal_v,std::vector<std::vector<double>> &vertical_v,std::vector<std::vector<double>> &pressure,std::vector<std::vector<double>>&gridM,std::vector<std::vector<Eigen::Vector2d>>&gridF){
        dx = x;
        dt = t;
        rho = density;
        u = horizontal_v;//v[nx+1][ny]
        v = vertical_v;//v[nx][ny+1]
        p = pressure;//p[nx][ny]
        mi = gridM;
        fi = gridF;
        L = dx*Nx;
        initPressure();
        std::cout << "initpressure" << std::endl;
        initForce();
        std::cout << "initforce" << std::endl;
    }
    void advect(){
        std::vector<std::vector<double>>old_u = u;
        std::vector<std::vector<double>>old_v = v;
        //水平成分
        for(int i=1;i<Nx;i++)for(int j=0;j<Ny;j++){
            //MAC格子座標でのサンプリング点の計算。MACグリッドの速度成分は、水平成分が(0.0,0.5)鉛直成分が(0.5,0.0)に設定してある。
            double x = i*dx;
            double y = (j + 0.5)*dx;
            //サンプリング点の移流
            double newx = x - dt*interpolation(Nx+1, Ny, x,y-0.5*dx, old_u);
            double newy = y - dt*interpolation(Nx, Ny+1, x-0.5*dx,y, old_v);
            
            u[i][j] = interpolation(Nx+1, Ny, newx, newy-0.5*dx, old_u);
        }
        for(int i=0;i<Nx;i++)for(int j=1;j<Ny;j++){
            //MAC格子座標でのサンプリング点の計算。MACグリッドの速度成分は、水平成分が(0.0,0.5)鉛直成分が(0.5,0.0)に設定してある。
            double x = (i+0.5)*dx;
            double y = j*dx;
            //サンプリング点の移流
            double newx = x - dt*interpolation(Nx+1, Ny, x-0.5*dx,y, old_u);
            double newy = y - dt*interpolation(Nx, Ny+1, x,y-0.5*dx, old_v);
            
            v[i][j] = interpolation(Nx, Ny+1, newx-0.5*dx, newy, old_v);
        }
    }
    //point=計算する点の座標
    //q=補完する離散値の二次元配列
    double interpolation(int nx,int ny,double x,double y,std::vector<std::vector<double>>&q){
        double s = fmax(0.0,fmin(nx-1-1e-6,x/dx));
        double t = fmax(0.0,fmin(ny-1-1e-6,y/dx));
    
        int i = s;
        int j = t;
        Eigen::Vector4d f = {q[i][j],q[i][j+1],q[i+1][j],q[i+1][j+1]};
        s -= i;
        t -= j;
        Eigen::Vector4d c = { (1-s)*(1-t) , (1-s)*t ,
                                s*(1-t) , s*t};
        return f.dot(c);
    }
    std::vector<int>BoundaryCondition(int i,int j,std::unordered_map<std::vector<int>,std::vector<int>,ArrayHasher<2>>&map){
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
        //128以上でオーバーフローする．
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
                if(map.find(key) == map.end())continue;
                //std::cout << i << "," << j << std::endl;
                double D[4] = {1.0,1.0,-1.0,-1.0};//周囲4方向に向かって働く、圧力の向き
                //double F[4] = {(double)(i<Nx-1),(double)(j<Ny-1),(double)(i>0),(double)(j>0)};//境界条件。壁なら0,流体なら1
                std::vector<int> F = BoundaryCondition(i,j,map);
                double P[4] = {//pn。周囲4つのセルの圧力値
                    (F[0] ? p[i+1][j] : 0.0),
                    (F[1] ? p[i][j+1] : 0.0),
                    (F[2] ? p[i-1][j] : 0.0),
                    (F[3] ? p[i][j-1] : 0.0)};
                
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
                    sum_L += F[n]*P[n]*scale;
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
        std::cout << "finGauss" << std::endl;
//        print_velocity();
        for(int i=1; i<Nx;i++)for(int j=0;j<Ny;j++){
            if(mi[i][j] < eps)u[i][j] =u[i][j] = u[i][j] - dt/rho * (p[i][j]-p[i-1][j])/dx;
            else u[i][j] = u[i][j] - dt/rho * (p[i][j]-p[i-1][j])/dx + dt*fi[i][j].x()/mi[i][j];
        }
        for(int i=0; i<Nx;i++)for(int j=1;j<Ny;j++){
            if(mi[i][j] < eps)v[i][j] = v[i][j] - dt/rho * (p[i][j]-p[i][j-1])/dx;
            else v[i][j] = v[i][j] - dt/rho * (p[i][j]-p[i][j-1])/dx + dt*fi[i][j].y()/mi[i][j];
        }
//        std::cout << "outputP" << std::endl;
//        print_pressure();
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
            if(i < Nx/2)p[i][j] = 0;
            //else p[i][j] = (double)(Nx-1-i)/Ny;
            else p[i][j] = 0;
        }
    }
    void initForce(){
        for(unsigned int i=0;i<Nx;i++)for(unsigned int j=0;j<Ny;j++){
            Eigen::Vector2d f0 = {0.0,-g0*dx};
            Eigen::Vector2d f1 = {0.0,-g0*dx/6};
            //std::cout << f0.x() << "," << f0.y() << std::endl;
            if(i>Nx/3 && i < Nx/3*2)fi[i][j] = f0;
            else fi[i][j] = f0;
        }
    }
};
#endif /* Fluid_h */
