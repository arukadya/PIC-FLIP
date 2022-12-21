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
#define alpha 0.5
#define mp 1 //粒子の重さ
#ifndef Flip_h
#define Flip_h

struct particle{
    Eigen::Vector2d PIC_velocity;
    Eigen::Vector2d FLIP_velocity;
    Eigen::Vector2d velocity;
    Eigen::Vector2d position;
    std::pair<int,int>gridIndex = std::make_pair(-1, -1);
    particle(Eigen::Vector2d v,Eigen::Vector2d p){
        Eigen::Vector2d zero = {0,0};
        velocity = v;
        PIC_velocity = zero;
        FLIP_velocity = zero;
        position = p;
    }
    void setGridIndex(int x,int y){
        gridIndex.first = x;
        gridIndex.second = y;
    }
};

struct PIC_FLIP : Fluid{
    std::vector<int> division;//division[0] = xの分割数.division[1]=y...
    std::vector<particle> particles;//入力メッシュの頂点
    std::unordered_map<std::vector<int>,std::vector<int>,ArrayHasher<2>>map;//ハッシュテーブル
    std::vector<Eigen::Vector2d> vertices;//出力メッシュの頂点
    
    void execute(){
        for(unsigned int i=0;i<repeatCount;i++){
            std::cout << i << std::endl;
            locateParticlesOnGrid(map);
            //preprocessingParticles();
            std::cout << "mapsize = " << map.size() << std::endl;
            particlesVelocityToGrid();
            std::cout << "P2G" << std::endl;
            calGridPressure();
            std::cout << "GridPressure" << std::endl;
            PICgridVelocityToParticles();
            FLIP_gridVelocityToParticles();
            std::cout << "G2P" << std::endl;
            advectParticles();
            std::cout << "advectParticles" << std::endl;
            output(vertices);
            std::string OutputVTK = "outputVTK/output"+std::to_string(i)+".vtk";
            std::ostringstream ssPressure;
            ssPressure << "outputPressure/output" << std::setw(3) << std::setfill('0') << i << ".dat";
            std::ostringstream ssMap;
            ssMap << "outputMap/output" << std::setw(3) << std::setfill('0') << i << ".dat";
            std::string OutputPressure(ssPressure.str());
            std::string OutputMap(ssMap.str());
            std::cout << OutputVTK.c_str() << std::endl;
            
            outputVTK(OutputVTK.c_str(),vertices);
            outputPLT_P(Nx, Ny, dx, OutputPressure.c_str(), p);
            outputPLT_M(Nx, Ny,OutputMap.c_str(), map);
        }
    }
    void output(std::vector<Eigen::Vector2d> &v){
        v.resize(particles.size());
        for(int i=0;i<v.size();i++)v[i] = particles[i].position;
    }
    PIC_FLIP(double x,double t,double density,std::vector<std::vector<double>> &horizontal_v,std::vector<std::vector<double>>&vertical_v,std::vector<std::vector<double>>&pressure,std::vector<std::vector<double>>&gridUM,std::vector<std::vector<double>>&gridVM,std::vector<std::vector<double>>&gridUF,std::vector<std::vector<double>>&gridVF):Fluid(x,t,density,horizontal_v,vertical_v,pressure,gridUM,gridVM,gridUF,gridVF){
        division = {Nx,Ny};
        initParticles();
    }
    
    void initParticles(){
        //std::cout << "initparticles" << std::endl;
//        for(unsigned int i=Nx/2;i<Nx;i++)for(unsigned int j=0;j<Ny;j++){
        for(unsigned int i=Nx/3;i<Nx/3*2;i++)for(unsigned int j=0;j<Ny;j++){
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
//----------------速度も初期化した方がいい気がする----------------------------------------------------------
//グリッドの切り方はよく考えるべき．スタッカート格子に則って，u,v別々にやるのが合ってる気がする．ブランチはきれ（戒め）
//-----------------------------------------------------------------------------------------------
            //P2Gの方法は諸説あり。
//-----------------------------------------------------------------------------------------------
        //初期化
        for(unsigned int i=0;i<Nx+1;i++)for(unsigned int j=0;j<Ny;j++){
            umi[i][j] = 0;
            old_u[i][j] = u[i][j];
            u[i][j] = 0;
            //delta_u[i][j] = 0;
        }
        for(unsigned int i=0;i<Nx;i++)for(unsigned int j=0;j<Ny+1;j++){
            vmi[i][j] = 0;
            old_v[i][j] = v[i][j];
            v[i][j] = 0;
            //delta_v[i][j] = 0;
        }
        //u
        for(int i=1;i<Nx+1;i++)for(int j=0;j<Ny;j++){
            std::vector<Eigen::Vector2d>gx_list = {{(i-1)*dx,(j+0.5)*dx},{(i)*dx,(j+0.5)*dx}};
            std::vector<std::vector<int>>key_list = {{i-1,j},{i,j}};
            if(i != 1){
                gx_list.push_back({(i-2)*dx,(j+0.5)*dx});
                key_list.push_back({i-2,j});
            }
            if(i != Nx){
                gx_list.push_back({(i+1)*dx,(j+0.5)*dx});
                key_list.push_back({i+1,j});
            }
            for(int k=0;k<gx_list.size();k++){
                if(map.find(key_list[k]) != map.end()){
                    auto val = map.at(key_list[k]);
                    for(auto x:val){
                        Eigen::Vector2d px = particles[x].position;
                        Eigen::Vector2d pv = particles[x].velocity;
                        double weight = weightFunction(px, gx_list[k], dx);
                        umi[i][j] += weight*mp;
                        u[i][j] += weight*mp*pv.x();
                    }
                    u[i][j] /= umi[i][j];
                }
                //else std::cout << key_list[k][0] << "," << key_list[k][1] << ":OutsideMap" << std::endl;
                u[i][j] += ufi[i][j] * dt;
            }
        }
        std::cout << u[Nx/2][0]/dx << "," <<u[Nx/2][1]/dx<< std::endl;
        //v
        for(int i=0;i<Nx;i++)for(int j=1;j<Ny+1;j++){
            std::vector<Eigen::Vector2d>gx_list = {{(i+0.5)*dx,(j)*dx},{(i+0.5)*dx,(j)*dx}};
            std::vector<std::vector<int>>key_list = {{i,j-1},{i,j}};
            if(j != 1){
                gx_list.push_back({(i+0.5)*dx,(j-2)*dx});
                key_list.push_back({i,j-2});
            }
            if(j != Ny){
                gx_list.push_back({(i+0.5)*dx,(j+1)*dx});
                key_list.push_back({i,j+1});
            }
            for(int k=0;k<gx_list.size();k++){
                if(map.find(key_list[k]) != map.end()){
                    auto val = map.at(key_list[k]);
                    for(auto x:val){
                        Eigen::Vector2d px = particles[x].position;
                        Eigen::Vector2d pv = particles[x].velocity;
                        double weight = weightFunction(px, gx_list[k], dx);
                        vmi[i][j] += weight*mp;
                        v[i][j] += weight*mp*pv.y();
                    }
                    v[i][j] /= vmi[i][j];
                }
            }
            v[i][j] += vfi[i][j] * dt;
        }
    }
    
    void calGridPressure(){
        for(int i=0;i<Nx;i++)for(int j=0;j<Ny;j++){
            std::vector<int>key = {i,j};
            if(map.find(key) == map.end()){
                p[i][j] = 0;
            }
        }
        project(map);
    }
    
    void PICgridVelocityToParticles(){
//-----------------------これも速度を初期化した方がいい気がする-----------------------------------------------
        //for(auto &x:particles)x.velocity = {0,0};
        //PIC_velocityの初期化が必要で，後はいらない
        for(auto &x:particles)x.PIC_velocity = {0,0};
        for(int i=1;i<Nx+1;i++)for(int j=0;j<Ny;j++){
            std::vector<Eigen::Vector2d>gx_list = {{(i-1)*dx,(j+0.5)*dx},{(i)*dx,(j+0.5)*dx}};
            std::vector<std::vector<int>>key_list = {{i-1,j},{i,j}};
            if(i != 1){
                gx_list.insert(gx_list.begin(),{(i-2)*dx,(j+0.5)*dx});
                key_list.insert(key_list.begin(),{i-2,j});
            }
            if(i != Nx){
                gx_list.push_back({(i+1)*dx,(j+0.5)*dx});
                key_list.push_back({i+1,j});
            }
            for(int k=0;k<gx_list.size();k++){
                if(map.find(key_list[k]) != map.end()){
                    auto val = map.at(key_list[k]);
                    for(auto x:val){
                        Eigen::Vector2d px = particles[x].position;
                        double weight = weightFunction(px, gx_list[k], dx);
                        if(i == 1)particles[x].PIC_velocity.x()+= weight*mp*u[i-1+k][j];
                        else particles[x].PIC_velocity.x()+= weight*mp*u[i-2+k][j];
                    }
                }
            }
        }
        for(int i=0;i<Nx;i++)for(int j=1;j<Ny+1;j++){
            std::vector<Eigen::Vector2d>gx_list = {{(i+0.5)*dx,(j)*dx},{(i+0.5)*dx,(j)*dx}};
            std::vector<std::vector<int>>key_list = {{i,j-1},{i,j}};
            if(j != 1){
                gx_list.insert(gx_list.begin(),{(i+0.5)*dx,(j-2)*dx});
                key_list.insert(key_list.begin(),{i,j-2});
            }
            if(j != Ny){
                gx_list.push_back({(i+0.5)*dx,(j+1)*dx});
                key_list.push_back({i,j+1});
            }
            for(int k=0;k<gx_list.size();k++){
                if(map.find(key_list[k]) != map.end()){
                    auto val = map.at(key_list[k]);
                    for(auto x:val){
                        Eigen::Vector2d px = particles[x].position;
                        double weight = weightFunction(px, gx_list[k], dx);
                        if(j == 1)particles[x].PIC_velocity.y()+= weight*mp*v[i][j-1+k];
                        else particles[x].PIC_velocity.y()+= weight*mp*v[i][j-2+k];
                    }
                }
            }
        }
    }
    void FLIP_gridVelocityToParticles(){
//-----------------------速度の初期化は不要-----------------------------------------------
        //for(auto &x:particles)x.velocity = {0,0};
        for(auto &x:particles)x.FLIP_velocity = x.velocity;
        for(int i=1;i<Nx+1;i++)for(int j=0;j<Ny;j++){
            std::vector<Eigen::Vector2d>gx_list = {{(i-1)*dx,(j+0.5)*dx},{(i)*dx,(j+0.5)*dx}};
            std::vector<std::vector<int>>key_list = {{i-1,j},{i,j}};
            if(i != 1){
                gx_list.insert(gx_list.begin(),{(i-2)*dx,(j+0.5)*dx});
                key_list.insert(key_list.begin(),{i-2,j});
            }
            if(i != Nx){
                gx_list.push_back({(i+1)*dx,(j+0.5)*dx});
                key_list.push_back({i+1,j});
            }
            for(int k=0;k<gx_list.size();k++){
                if(map.find(key_list[k]) != map.end()){
                    auto val = map.at(key_list[k]);
                    for(auto x:val){
                        Eigen::Vector2d px = particles[x].position;
                        double weight = weightFunction(px, gx_list[k], dx);
//                        if(i == 1)particles[x].FLIP_velocity.x()+= weight*mp*delta_u[i-1+k][j];
//                        else particles[x].FLIP_velocity.x()+= weight*mp*delta_u[i-2+k][j];
                        if(i == 1)particles[x].FLIP_velocity.x()+= weight*mp*(u[i-1+k][j]-old_u[i-1+k][j]);
                        else particles[x].FLIP_velocity.x()+= weight*mp*(u[i-2+k][j]-old_u[i-2+k][j]);
                    }
                }
            }
        }
        for(int i=0;i<Nx;i++)for(int j=1;j<Ny+1;j++){
            std::vector<Eigen::Vector2d>gx_list = {{(i+0.5)*dx,(j)*dx},{(i+0.5)*dx,(j)*dx}};
            std::vector<std::vector<int>>key_list = {{i,j-1},{i,j}};
            if(j != 1){
                gx_list.insert(gx_list.begin(),{(i+0.5)*dx,(j-2)*dx});
                key_list.insert(key_list.begin(),{i,j-2});
            }
            if(j != Ny){
                gx_list.push_back({(i+0.5)*dx,(j+1)*dx});
                key_list.push_back({i,j+1});
            }
            for(int k=0;k<gx_list.size();k++){
                if(map.find(key_list[k]) != map.end()){
                    auto val = map.at(key_list[k]);
                    for(auto x:val){
                        Eigen::Vector2d px = particles[x].position;
                        double weight = weightFunction(px, gx_list[k], dx);
//                        if(j == 1)particles[x].FLIP_velocity.y()+= weight*mp*delta_v[i][j-1+k];
//                        else particles[x].FLIP_velocity.y()+= weight*mp*delta_v[i][j-2+k];
                        if(j == 1)particles[x].FLIP_velocity.y()+= weight*mp*(v[i][j-1+k]-old_v[i][j-1+k]);
                        else particles[x].FLIP_velocity.y()+= weight*mp*(v[i][j-2+k]-old_v[i][j-2+k]);
                    }
                }
            }
        }
    }
    void advectParticles(){
        //std::cout << L << std::endl;
        for(int i=0;i<particles.size();i++){
            //std::cout <<"x0:"<< particles[i].position.x() << " " << particles[i].position.y() << std::endl;
            particles[i].velocity = alpha*particles[i].FLIP_velocity +(1 - alpha)*particles[i].PIC_velocity;
            particles[i].position.x() += particles[i].velocity.x()*dt;
            particles[i].position.y() += particles[i].velocity.y()*dt;
            pushout(particles[i].position, L,dx);
            //std::cout <<"x1:"<< particles[i].position.x() << " " << particles[i].position.y() << std::endl;
        }
    }
//    void preprocessingParticles(){
//        locateParticlesOnGrid(map);
//        for(int i=0;i<particles.size();i++){
//            int keyx = particles[i].gridIndex.first;
//            int keyy = particles[i].gridIndex.second;
//            particles[i].velocity.x() += ufi[keyx][keyy]*dt;
//            particles[i].velocity.y() += vfi[keyx][keyy]*dt;
//            //pushout(particles[i].position, L,dx);
//        }
//        //locateParticlesOnGrid(map);
//    }
};
#endif /* Flip_h */
