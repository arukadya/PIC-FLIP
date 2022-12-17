#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <iomanip>
#include <unordered_map>
#include "Hasher.h"//Hash関数
#include <set>
#include <functional>//std::hash
#include <type_traits>//std::remove_cvref_t(C++20)
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#include "Fluid.h"
#include "vtk.h"
#include "gnuplot.h"
#include "functions.h"
#define repeatCount 100
#define mp 1 //粒子の重さ
#ifndef Flip_h
#define Flip_h

struct particle{
    Eigen::Vector2d velocity;
    Eigen::Vector2d position;
    std::pair<int,int>gridIndex = std::make_pair(-1, -1);
    particle(Eigen::Vector2d v,Eigen::Vector2d p){
        velocity = v;
        position = p;
    }
    void setGridIndex(int x,int y){
        gridIndex.first = x;
        gridIndex.second = y;
    }
};

struct P2GG2P : Fluid{
    std::vector<int> division;//division[0] = xの分割数.division[1]=y...
    std::vector<particle> particles;//入力メッシュの頂点
    std::unordered_map<std::vector<int>,std::vector<int>,ArrayHasher<2>>map;//ハッシュテーブル
    std::vector<Eigen::Vector2d> vertices;//出力メッシュの頂点
    
    void execute(){
        for(unsigned int i=0;i<repeatCount;i++){
            std::cout << i << std::endl;
            locateParticlesOnGrid(map);
            std::cout << "mapsize = " << map.size() << std::endl;
            particlesVelocityToGrid();
            std::cout << "P2G" << std::endl;
            calGridPressure();
            std::cout << "GridPressure" << std::endl;
            gridVelocityToParticles();
            std::cout << "G2P" << std::endl;
            moveParticles();
            std::cout << "moveParticles" << std::endl;
            output(vertices);
            std::string OutputVTK = "outputVTK/output"+std::to_string(i)+".vtk";
            std::ostringstream ssdat;
            ssdat << "outputDAT/output" << std::setw(3) << std::setfill('0') << i << ".dat";
            //std::string OutputDAT = "outputDAT/output"+std::to_string(i)+".dat";
            std::string OutputDAT(ssdat.str());
            std::cout << OutputVTK.c_str() << std::endl;
            outputVTK(OutputVTK.c_str(),vertices);
            outputPLT_P(Nx, Ny, dx, OutputDAT.c_str(), p);
        }
    }
    void output(std::vector<Eigen::Vector2d> &v){
        v.resize(particles.size());
        for(int i=0;i<v.size();i++)v[i] = particles[i].position;
    }
    P2GG2P(double x,double t,double density,std::vector<std::vector<double>> &horizontal_v,std::vector<std::vector<double>>&vertical_v,std::vector<std::vector<double>>&pressure,std::vector<std::vector<double>>&gridM,std::vector<std::vector<Eigen::Vector2d>>&gridF):Fluid(x,t,density,horizontal_v,vertical_v,pressure,gridM,gridF){
        division = {Nx,Ny};
        initParticles();
    }
    
    void initParticles(){
        std::cout << "initparticles" << std::endl;
        for(unsigned int i=Nx/2;i<Nx;i++)for(unsigned int j=0;j<Ny;j++){
            Eigen::Vector2d v0 = {0.0,0};
            std::vector<Eigen::Vector2d>pos(4);//1グリッドにn^2個置くのが流儀らしい
            pos[0] = Eigen::Vector2d{(i+0.25)*dx,(j+0.25)*dx};
            pos[1] = Eigen::Vector2d{(i+0.75)*dx,(j+0.25)*dx};
            pos[2] = Eigen::Vector2d{(i+0.25)*dx,(j+0.75)*dx};
            pos[3] = Eigen::Vector2d{(i+0.75)*dx,(j+0.75)*dx};
            for(unsigned int k=0;k<4;k++){
                //if(i==Nx-1 && j == Ny/2 - 1)std::cout << pos[k].x() << "," << pos[k].y() << std::endl;
                particle tmp = particle(v0,pos[k]);
                particles.push_back(tmp);
                weights.push_back(-1);
            }
        }
    }
    
    void locateParticlesOnGrid(std::unordered_map<std::vector<int>,std::vector<int>,ArrayHasher<2>>&map){
        std::unordered_map<std::vector<int>,std::vector<int>,ArrayHasher<2>>m;
        for(int i=0;i<particles.size();i++){
            std::vector<int>key = calGridOrKey(particles[i].position);
            if(m.find(key) == m.end()){
                std::vector<int>contain = {i};
                m.emplace(key,contain);
            }
            else{
                m.at(key).push_back(i);
            }
            particles[i].setGridIndex(key[0], key[1]);
        }
        map=m;
        //std::cout << "end locating" << std::endl;
    }
    Eigen::Vector2d calCellSize(){
        return Eigen::Vector2d{dx,dx};
    }
    std::vector<int> calGridOrKey(Eigen::Vector2d v){
        std::vector<int> key;
        Eigen::Vector2d cellSize = calCellSize();
        key.push_back( floor(( (v.x() ) / cellSize.x() ) ));
        key.push_back( floor(( (v.y() ) / cellSize.y() ) ));
        return key;
    }
    void particlesVelocityToGrid(){
        for(unsigned int i=0;i<Nx;i++)for(unsigned int j=0;j<Ny;j++){
            mi[i][j] = 0;
        }
        for(const auto &[key,val] : map){
            if(key[0] < 0 || key[0] >=Nx || key[1] < 0 || key[1] >=Ny){
                std::cout << "outside map!!" << std::endl;
                std::cout << key[0] << " " << key[1] << std::endl;
//                for(auto &x:val){
//                    std::cout << "particle" << x << "(" << particles[x].position.x() << "," << particles[x].position.y() << ")" << std::endl;
//                }
                continue;
            }
            double gridU=0.0;
            double gridV=0.0;
            Eigen::Vector2d gridPos = {(key[0]+0.5)*dx,(key[1]+0.5)*dx};
            //std::cout << gridPos << std::endl;
            /*for(auto &x:val){
                //ただ平均とるだけ
                gridU += particles[x].velocity.x();
                gridV += particles[x].velocity.y();
            }
             u[key[0]][key[1]] = gridU/val.size();
             v[key[0]][key[1]] = gridV/val.size();
             */
//-----------------------------------------------------------------------------------------------
            //P2Gの方法は諸説あり。
            double sum_mp = 0;
            for(auto &x:val){
                particles[x].setGridIndex(key[0], key[1]);
                weights[x] = weightFunction(particles[x].position, gridPos, dx);
                //std::cout << particles[x].position << std::endl;
                //std::cout << weights[x] << std::endl;
                gridU += weights[x]*particles[x].velocity.x()*mp;
                gridV += weights[x]*particles[x].velocity.y()*mp;
                sum_mp += weights[x]*mp;
            }
            u[key[0]][key[1]] = gridU;
            v[key[0]][key[1]] = gridV;
            mi[key[0]][key[1]] = sum_mp;
            //std::cout << u[key[0]][key[1]] << " " << v[key[0]][key[1]] << std::endl;
//-----------------------------------------------------------------------------------------------
        }
        //print_velocity();
    }
    void calGridPressure(){
        for(int i=0;i<Nx;i++)for(int j=0;j<Ny;j++){
            std::vector<int>key = {i,j};
            if(map.find(key) == map.end()){
                p[i][j] = 0;
            }
        }
        project(map);
        //std::cout << "pressure" << std::endl;
        //print_pressure();
        //std::cout << "velocity" << std::endl;
        //print_velocity();
    }
    
    void gridVelocityToParticles(){
        //print_velocity();
        for(int i=0;i<particles.size();i++){
            //付近４グリッドの線形補完
//            particles[i].velocity.x() = interpolation(Nx+1, Ny, particles[i].position.x(), particles[i].position.y(), u);
//            particles[i].velocity.y() = interpolation(Nx, Ny+1, particles[i].position.x(), particles[i].position.y(), v);
            //重みを使った逆算
            particles[i].velocity.x() = weights[i]*u[particles[i].gridIndex.first][particles[i].gridIndex.second];
            particles[i].velocity.y() = weights[i]*v[particles[i].gridIndex.first][particles[i].gridIndex.second];
        }
    }
    void moveParticles(){
        //std::cout << L << std::endl;
        for(int i=0;i<particles.size();i++){
            particles[i].position.x() += particles[i].velocity.x()*dt;
            particles[i].position.y() += particles[i].velocity.y()*dt;
            pushout(particles[i].position, L,dx);
            //std::cout <<"x1:"<< particles[i].position.x() << " " << particles[i].position.y() << std::endl;
        }
    }
};
#endif /* Flip_h */
