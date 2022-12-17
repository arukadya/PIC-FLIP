//
//  mable.h
//  Mable
//
//  Created by 須之内俊樹 on 2022/11/08.
//
#ifndef mable_h
#define mable_h
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

struct Mable{
    double dx;//セルの大きさ
    double dt;//時間の刻み幅
    double rho;
    double** u;//水平
    double** v;//鉛直
    double** p;//圧力
    
    Mable(double x,double t,double density,double** horizontal_v,double** vertical_v,double** pressure){
        dx = x;
        dt = t;
        rho = density;
        u = horizontal_v;//v[nx+1][ny]
        v = vertical_v;//v[nx][ny+1]
        p = pressure;//p[nx][ny]
    }
    void advect(int Nx,int Ny){
        double** old_u = u;
        double** old_v = v;
        //水平成分
        for(unsigned int i=1;i<Nx;i++)for(unsigned int j=0;j<Ny;j++){
            //MAC格子座標でのサンプリング点の計算。MACグリッドの速度成分は、水平成分が(0.0,0.5)鉛直成分が(0.5,0.0)に設定してある。
            double x = i*dx;
            double y = (j + 0.5)*dx;
            //サンプリング点の移流
            double newx = x - dt*interpolation(Nx+1, Ny, x,y-0.5*dx, old_u);
            double newy = y - dt*interpolation(Nx, Ny+1, x-0.5*dx,y, old_v);
            
            u[i][j] = interpolation(Nx+1, Ny, newx, newy-0.5*dx, old_u);
        }
        for(unsigned int i=0;i<Nx;i++)for(unsigned int j=1;j<Ny;j++){
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
    double interpolation(int Nx,int Ny,double x,double y,double** q){
        double s = fmax(0.0,fmin(Nx-1-1e-6,x/dx));
        double t = fmax(0.0,fmin(Ny-1-1e-6,y/dx));
        
        unsigned int i = s;
        unsigned int j = t;
        Eigen::Vector4d f = {q[i][j],q[i][j+1],q[i+1][j],q[i+1][j+1]};
        s -= i;
        t -= j;
        Eigen::Vector4d c = { (1-s)*(1-t) , (1-s)*t ,
                                s*(1-t) , s*t};
        return f.dot(c);
    }
    
    
    
    void project(unsigned int Nx,unsigned int Ny){
        double scale = dt/(rho*dx*dx);//左辺の係数部分
        double eps = 1.0e-4;//ガウスザイデル法の精度
        double err;//ガウスザイデル法の残差
        do{
        //圧力場を求める
            err = 0.0;
            for(unsigned int j=0;j<Ny;j++)for(unsigned int i=0;i<Nx;i++){
                double D[4] = {1.0,1.0,-1.0,-1.0};//周囲4方向に向かって働く、圧力の向き
                double F[4] = {(double)(i<Nx-1),(double)(j<Ny-1),(double)(i>0),(double)(j>0)};//境界条件。壁なら0,流体なら1
                //int F[4] = {i<Nx-1,j<Ny-1,i>0,j>0};
                double P[4] = {//pn。周囲4つのセルの圧力値
                    (F[0] ? p[i+1][j] : 0.0),
                    (F[1] ? p[i][j+1] : 0.0),
                    (F[2] ? p[i-1][j] : 0.0),
                    (F[3] ? p[i][j-1] : 0.0)};
                double U[4] = {u[i+1][j],v[i][j+1],u[i][j],v[i][j]};
                //セルの圧力値を求める
                double det = 0.0;//左辺のdeterminant
                double sum_L = 0.0;
                double sum_R = 0.0;
                for(unsigned int n=0;n<4;n++){
                    det += F[n]*scale;
                    sum_L += F[n]*P[n]*scale;
                    sum_R += F[n]*D[n]*U[n]*dx;
                }
                err = fmax(err,fabs(det*p[i][j]-sum_L+sum_R));
                p[i][j] = (sum_L-sum_R)/det;
            }
        }while(eps<err);//反復回数は初期値と収束速度？に依存。定数回なら全体で計算量はO(n^2)
        //新しい流速を求める。u(t+1) = u* - (dt/rho)grad(p)
        for(unsigned int i=1; i<Nx;i++)for(unsigned int j=0;j<Ny;j++){
            u[i][j] = u[i][j] - dt/rho * (p[i][j]-p[i-1][j])/dx;
        }
        for(unsigned int i=0; i<Nx;i++)for(unsigned int j=1;j<Ny;j++){
            v[i][j] = v[i][j] - dt/rho * (p[i][j]-p[i][j-1])/dx;
        }
    }
};

#endif /* mable_h */
