#include "Hasher.h"//Hash関数
#include <functional>//std::hash
#include <type_traits>//std::remove_cvref_t(C++20)
//#include <Eigen/Core>
//#include <Eigen/Dense>
//#include <Eigen/LU>
#include "Fluid.h"
#include "vtk.h"
#include "gnuplot.h"
#include "functions.h"
#include "particle.h"
#define repeatCount 1000
#define alpha 1
//#define mp  //粒子の重さ
//#define radius 0.0025
#define gamma 1
#ifndef Flip_h
#define Flip_h
#define timer 2
#define extend 0

struct PIC_FLIP : Fluid{
    std::vector<int> division;//division[0] = xの分割数.division[1]=y...
    std::vector<particle>particles;//入力メッシュの頂点
    std::unordered_map<std::vector<int>,std::vector<int>,ArrayHasher<3>>map;//ハッシュテーブル
    std::vector<Eigen::Vector3d> vertices;//出力メッシュの頂点
    timeDisplayer TD;
    double radius = dx*2;
    double mp = pow(radius,3)/3*4*3.14;
    void execute(){
        int cnt = 0;
        for(unsigned int i=0;i<repeatCount;i++){
            //std::cout << i << std::endl;
            if(timer == 2)TD.startTimer("execute");
            //locateParticlesOnGrid(map);
            if(timer == 1)TD.startTimer("preprocess");
            preprocessingParticles();
            if(timer == 1)TD.endTimer();
            //if(timer)TD.startTimer("P2G");
            particlesVelocityToGrid();
            //if(timer)TD.endTimer();
            if(timer == 1)TD.startTimer("GridPressure");
            calGridPressure();
            //std::cout << "pressure" << std::endl;
            if(timer == 1)TD.endTimer();
            //if(timer)TD.startTimer("G2P_PIC");
            gridVelocityToParticles();
            //std::cout << "G2P" << std::endl;
            //if(timer)TD.endTimer();
            //if(timer)TD.startTimer("G2P_FLIP");
            //FLIPgridVelocityToParticles();
            //if(timer)TD.endTimer();
            //if(timer)TD.startTimer("advect");
            advectParticles();
            //if(timer)TD.endTimer();
            
            if(i%(repeatCount/100) == 0){
                if(timer == 2){
                    TD.endTimer();
                    TD.startTimer("output");
                }
                //std::cout << i << std::endl;
                std::cout << "mapsize = " << map.size() << std::endl;
                output(vertices);
                std::string OutputVTK = "outputVTK/output"+std::to_string(cnt)+".vtk";
                //std::ostringstream ssPressure;
                //ssPressure << "outputPressure/output" << std::setw(3) << std::setfill('0') << cnt << ".dat";
                //std::ostringstream ssMap;
                //ssMap << "outputMap/output" << std::setw(3) << std::setfill('0') << cnt << ".dat";
                //std::ostringstream ssParticles;
                //ssParticles << "outputParticles/output" << std::setw(3) << std::setfill('0') << cnt << ".dat";
                //std::string OutputPressure(ssPressure.str());
                //std::string OutputMap(ssMap.str());
                //std::string OutputParticles(ssParticles.str());
                std::cout << OutputVTK.c_str() << std::endl;
                outputVTK(OutputVTK.c_str(),vertices);
                //outputPLT_P(Nx, Ny, dx, OutputPressure.c_str(), p);
                //outputPLT_M(Nx, Ny,OutputMap.c_str(), map);
                //outputPLT_particles(OutputParticles.c_str(),particles);
                cnt++;
                if(timer == 2)TD.endTimer();
            }
        }
        std::cout << "dx:" << dx << " dt:" << dt << " Nx*Ny:" << Nx*Ny << " repeat:" << repeatCount << std::endl;
        std::cout << "alpha:" << alpha << " gamma:" << gamma << " radius:" << radius <<std::endl;
        std::cout << "extend:" << extend << std::endl;
    }
    void output(std::vector<Eigen::Vector3d> &v){
        v.resize(particles.size());
        for(int i=0;i<v.size();i++)v[i] = particles[i].position;
    }
    PIC_FLIP(double x,double t,double density):Fluid(x,t,density){
        division = {Nx,Ny,Nz};
        initParticles();
    }
    void initParticles(){
        for(unsigned int i=0;i<Nx;i++){
            for(unsigned int j=0;j<Ny;j++){
                for(unsigned int k=0;k<Nz;k++){
                    if((i>Nx/4 && i < Nx/4 * 3) && (j > Ny/4) && ((k > Nz/4 && k < Nz/4*3))){
                        Eigen::Vector3d v0 = {0.0,0.0,0.0};
                        std::vector<Eigen::Vector3d>pos(8);//1グリッドにn^2個置くのが流儀らしい
                        pos[0] = Eigen::Vector3d{(i+0.25)*dx,(j+0.25)*dx,(k+0.25)*dx};
                        pos[1] = Eigen::Vector3d{(i+0.75)*dx,(j+0.25)*dx,(k+0.25)*dx};
                        pos[2] = Eigen::Vector3d{(i+0.25)*dx,(j+0.75)*dx,(k+0.25)*dx};
                        pos[3] = Eigen::Vector3d{(i+0.75)*dx,(j+0.75)*dx,(k+0.25)*dx};
                        pos[4] = Eigen::Vector3d{(i+0.25)*dx,(j+0.25)*dx,(k+0.75)*dx};
                        pos[5] = Eigen::Vector3d{(i+0.75)*dx,(j+0.25)*dx,(k+0.75)*dx};
                        pos[6] = Eigen::Vector3d{(i+0.25)*dx,(j+0.75)*dx,(k+0.75)*dx};
                        pos[7] = Eigen::Vector3d{(i+0.75)*dx,(j+0.75)*dx,(k+0.75)*dx};
                        for(unsigned int k=0;k<8;k++){
                            //if(i==Nx-1 && j == Ny/2 - 1)std::cout << pos[k].x() << "," << pos[k].y() << std::endl;
                            particle tmp = particle(v0,pos[k]);
                            particles.push_back(tmp);
                            weights.push_back(-1);
                        }
                    }
                }
            }
        }
    }

    void locateParticlesOnGrid(std::unordered_map<std::vector<int>,std::vector<int>,ArrayHasher<3>>&map){
        std::unordered_map<std::vector<int>,std::vector<int>,ArrayHasher<3>>m;
        for(int i=0;i<particles.size();i++){
            std::vector<int>key = calGridOrKey(particles[i].position);
            if(m.find(key) == m.end()){
                std::vector<int>contain = {i};
                m.emplace(key,contain);
            }
            else{
                m.at(key).push_back(i);
            }
            particles[i].setGridIndex(key[0], key[1], key[2]);
        }
        map=m;
        //std::cout << "end locating" << std::endl;
    }
    Eigen::Vector3d calCellSize(){
        return Eigen::Vector3d{dx,dx,dx};
    }
    std::vector<int> calGridOrKey(Eigen::Vector3d v){
        std::vector<int> key;
        Eigen::Vector3d cellSize = calCellSize();
        key.push_back( floor(( (v.x() ) / cellSize.x() ) ));
        key.push_back( floor(( (v.y() ) / cellSize.y() ) ));
        key.push_back( floor(( (v.z() ) / cellSize.z() ) ));
        return key;
    }
    void particlesVelocityToGrid(){
//----------------速度も初期化した方がいい気がする----------------------------------------------------------
//グリッドの切り方はよく考えるべき．スタッカート格子に則って，u,v別々にやるのが合ってる気がする．ブランチはきれ（戒め）
//-----------------------------------------------------------------------------------------------
            //P2Gの方法は諸説あり。
//-----------------------------------------------------------------------------------------------
        //初期化
        //copyVelocity();
        for(unsigned int i=0;i<Nx+1;i++){
            for(unsigned int j=0;j<Ny;j++){
                for(unsigned int k=0;k<Nz;k++){
                    umi.value[i][j][k] = 0;
                    old_u.value[i][j][k] = 0;
                    u.value[i][j][k]= 0;
                }
            }
        }
        for(unsigned int i=0;i<Nx;i++){
            for(unsigned int j=0;j<Ny+1;j++){
                for(unsigned int k=0;k<Nz;k++){
                vmi.value[i][j][k] = 0;
                old_v.value[i][j][k] = 0;
                v.value[i][j][k] = 0;
                }
            }
        }
        for(unsigned int i=0;i<Nx;i++){
            for(unsigned int j=0;j<Ny;j++){
                for(unsigned int k=0;k<Nz+1;k++){
                wmi.value[i][j][k] = 0;
                old_w.value[i][j][k] = 0;
                w.value[i][j][k] = 0;
                }
            }
        }
        //u
        for(int i=1;i<Nx+1;i++){
            for(int j=0;j<Ny;j++){
                for(int k=0;k<Nz;k++){
                    std::vector<Eigen::Vector3d>gx_list = {{(i-1)*dx,(j+0.5)*dx,(k+0.5)*dx},{(i)*dx,(j+0.5)*dx,(k+0.5)*dx}};
                    std::vector<std::vector<int>>key_list = {{i-1,j,k},{i,j,k}};
                    bool flg = false;
                    if(extend){
                        if(i != 1){
                            gx_list.push_back({(i-2)*dx,(j+0.5)*dx,(k+0.5)*dx});
                            key_list.push_back({i-2,j,k});
                        }
                        if(i != Nx){
                            gx_list.push_back({(i+1)*dx,(j+0.5)*dx,(k+0.5)*dx});
                            key_list.push_back({i+1,j,k});
                        }
                    }
                    
                    for(int l=0;l<gx_list.size();l++){
                        if(map.find(key_list[l]) != map.end()){
                            flg = true;
                            auto val = map.at(key_list[l]);
                            for(auto x:val){
                                Eigen::Vector3d px = particles[x].position;
                                Eigen::Vector3d pv = particles[x].velocity;
                                double weight = weightFunction(px, gx_list[l], dx);
                                umi.value[i][j][k] += weight*mp;
                                u.value[i][j][k] += weight*mp*pv.x();
                            }
                        }
                    }
                    if(flg){
                        u.value[i][j][k] /= umi.value[i][j][k];
                        old_u.value[i][j][k] = u.value[i][j][k];
                        u.value[i][j][k] += ufi.value[i][j][k] * dt/umi.value[i][j][k];
                    }
                }
            }
        }
        //v
        for(int i=0;i<Nx;i++){
            for(int j=1;j<Ny+1;j++){
                for(int k=0;k<Nz;k++){
                    std::vector<Eigen::Vector3d>gx_list = {{(i+0.5)*dx,(j-1)*dx,(k+0.5)*dx},{(i+0.5)*dx,(j)*dx,(k+0.5)*dx}};
                    std::vector<std::vector<int>>key_list = {{i,j-1,k},{i,j,k}};
                    bool flg = false;
                    if(extend){
                        if(j != 1){
                            gx_list.push_back({(i+0.5)*dx,(j-2)*dx,(k+0.5)*dx});
                            key_list.push_back({i,j-2,k});
                        }
                        if(j != Ny){
                            gx_list.push_back({(i+0.5)*dx,(j+1)*dx,(k+0.5)*dx});
                            key_list.push_back({i,j+1,k});
                        }
                    }
                    for(int l=0;l<gx_list.size();l++){
                        if(map.find(key_list[l]) != map.end()){
                            flg = true;
                            auto val = map.at(key_list[l]);
                            for(auto x:val){
                                Eigen::Vector3d px = particles[x].position;
                                Eigen::Vector3d pv = particles[x].velocity;
                                double weight = weightFunction(px, gx_list[l], dx);
                                vmi.value[i][j][k] += weight*mp;
                                v.value[i][j][k] += weight*mp*pv.y();
                            }
                        }
                    }
                    if(flg){
                        v.value[i][j][k] /= vmi.value[i][j][k];
                        old_v.value[i][j][k] = v.value[i][j][k];
                        v.value[i][j][k] += vfi.value[i][j][k] * dt/vmi.value[i][j][k];
                    }
                }
            }
        }
        for(int i=0;i<Nx;i++){
            for(int j=0;j<Ny;j++){
                for(int k=1;k<Nz+1;k++){
                    std::vector<Eigen::Vector3d>gx_list = {{(i+0.5)*dx,(j+0.5)*dx,(k-1)*dx},{(i+0.5)*dx,(j+0.5)*dx,(k)*dx}};
                    std::vector<std::vector<int>>key_list = {{i,j,k-1},{i,j,k}};
                    bool flg = false;
                    if(extend){
                        if(j != 1){
                            gx_list.push_back({(i+0.5)*dx,(j+0.5)*dx,(k-2)*dx});
                            key_list.push_back({i,j,k-2});
                        }
                        if(j != Ny){
                            gx_list.push_back({(i+0.5)*dx,(j+0.5)*dx,(k+1)*dx});
                            key_list.push_back({i,j,k+1});
                        }
                    }
                    for(int l=0;l<gx_list.size();l++){
                        if(map.find(key_list[l]) != map.end()){
                            flg = true;
                            auto val = map.at(key_list[l]);
                            for(auto x:val){
                                Eigen::Vector3d px = particles[x].position;
                                Eigen::Vector3d pv = particles[x].velocity;
                                double weight = weightFunction(px, gx_list[l], dx);
                                wmi.value[i][j][k] += weight*mp;
                                w.value[i][j][k] += weight*mp*pv.z();
                            }
                        }
                    }
                    if(flg){
                        w.value[i][j][k] /= wmi.value[i][j][k];
                        old_w.value[i][j][k] = w.value[i][j][k];
                        w.value[i][j][k] += wfi.value[i][j][k] * dt/wmi.value[i][j][k];
                    }
                }
            }
        }
    }
    
    void calGridPressure(){
        //copyVelocity();
        for(int i=0;i<Nx;i++){
            for(int j=0;j<Ny;j++){
                for(int k=0;k<Nz;k++){
                    std::vector<int>key = {i,j};
                    if(map.find(key) == map.end()){
                        p.value[i][j][k] = 0;
                    }
                }
            }
        }
        project(map);
    }
    
    void gridVelocityToParticles(){
//-----------------------これも速度を初期化した方がいい気がする-----------------------------------------------
        for(auto &p:particles){
            p.PIC_velocity = {0,0,0};
            p.FLIP_velocity = p.velocity;
        }
        for(int i=1;i<Nx+1;i++){
            for(int j=0;j<Ny;j++){
                for(int k=0;k<Nz;k++){
                    std::vector<Eigen::Vector3d>gx_list = {{(i-1)*dx,(j+0.5)*dx,(k+0.5)*dx},{(i)*dx,(j+0.5)*dx,(k+0.5)*dx}};
                    std::vector<std::vector<int>>key_list = {{i-1,j,k},{i,j,k}};
                    if(extend){
                        if(i != 1){
                            gx_list.push_back({(i-2)*dx,(j+0.5)*dx,(k+0.5)*dx});
                            key_list.push_back({i-2,j,k});
                        }
                        if(i != Nx){
                            gx_list.push_back({(i+1)*dx,(j+0.5)*dx,(k+0.5)*dx});
                            key_list.push_back({i+1,j,k});
                        }
                    }
                    for(int l=0;l<gx_list.size();l++){
                        if(map.find(key_list[l]) != map.end()){
                            auto val = map.at(key_list[l]);
                            for(auto x:val){
                                Eigen::Vector3d px = particles[x].position;
                                double weight = weightFunction(px, gx_list[l], dx);
                                if(i == 1){
                                    particles[x].PIC_velocity.x()+= weight*mp*u.value[i-1+l][j][k];
                                    particles[x].FLIP_velocity.x()+= weight*mp*(u.value[i-1+l][j][k]-old_u.value[i-1+l][j][k]);
                                }
                                else{
                                    particles[x].PIC_velocity.x()+= weight*mp*u.value[i-2+l][j][k];
                                    particles[x].FLIP_velocity.x()+= weight*mp*(u.value[i-1+l][j][k]-old_u.value[i-2+l][j][k]);
                                }
                            }
                        }
                    }
                }
            }
        }
        for(int i=0;i<Nx;i++){
            for(int j=1;j<Ny+1;j++){
                for(int k=0;k<Nz;k++){
                    std::vector<Eigen::Vector3d>gx_list = {{(i+0.5)*dx,(j-1)*dx,(k+0.5)*dx},{(i+0.5)*dx,(j)*dx,(k+0.5)*dx}};
                    std::vector<std::vector<int>>key_list = {{i,j-1,k},{i,j,k}};
                    if(extend){
                        if(j != 1){
                            gx_list.push_back({(i+0.5)*dx,(j-2)*dx,(k+0.5)*dx});
                            key_list.push_back({i,j-2,k});
                        }
                        if(j != Ny){
                            gx_list.push_back({(i+0.5)*dx,(j+1)*dx,(k+0.5)*dx});
                            key_list.push_back({i,j+1,k});
                        }
                    }
                    for(int l=0;l<gx_list.size();l++){
                        if(map.find(key_list[l]) != map.end()){
                            auto val = map.at(key_list[l]);
                            for(auto x:val){
                                Eigen::Vector3d px = particles[x].position;
                                double weight = weightFunction(px, gx_list[l], dx);
                                if(j == 1){
                                    particles[x].PIC_velocity.y()+= weight*mp*v.value[i][j-1+l][k];
                                    particles[x].FLIP_velocity.y()+= weight*mp*(v.value[i][j-1+l][k]-old_v.value[i][j-1+l][k]);
                                }
                                else{
                                    particles[x].PIC_velocity.y()+= weight*mp*v.value[i][j-2+l][k];
                                    particles[x].FLIP_velocity.y()+= weight*mp*(v.value[i][j-2+l][k]-old_v.value[i][j-2+l][k]);
                                }
                            }
                        }
                    }
                }
            }
        }
        for(int i=0;i<Nx;i++){
            for(int j=0;j<Ny;j++){
                for(int k=1;k<Nz+1;k++){
                    std::vector<Eigen::Vector3d>gx_list = {{(i+0.5)*dx,(j+0.5)*dx,(k-1)*dx},{(i+0.5)*dx,(j+0.5)*dx,(k)*dx}};
                    std::vector<std::vector<int>>key_list = {{i,j,k-1},{i,j,k}};
                    if(extend){
                        if(j != 1){
                            gx_list.push_back({(i+0.5)*dx,(j+0.5)*dx,(k-2)*dx});
                            key_list.push_back({i,j,k-2});
                        }
                        if(j != Ny){
                            gx_list.push_back({(i+0.5)*dx,(j+0.5)*dx,(k+1)*dx});
                            key_list.push_back({i,j,k+1});
                        }
                    }
                    for(int l=0;l<gx_list.size();l++){
                        if(map.find(key_list[l]) != map.end()){
                            auto val = map.at(key_list[l]);
                            for(auto x:val){
                                Eigen::Vector3d px = particles[x].position;
                                double weight = weightFunction(px, gx_list[l], dx);
                                if(k == 1){
                                    particles[x].PIC_velocity.z()+= weight*mp*w.value[i][j][k-1+l];
                                    particles[x].FLIP_velocity.z()+= weight*mp*(w.value[i][j][k-1+l]-old_w.value[i][j][k-1+l]);
                                }
                                else {
                                    particles[x].PIC_velocity.z()+= weight*mp*w.value[i][j][k-2+l];
                                    particles[x].FLIP_velocity.z()+= weight*mp*(w.value[i][j][k-2+l]-old_w.value[i][j][k-2+l]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
//    void FLIPgridVelocityToParticles(){
////-----------------------速度の初期化は不要-----------------------------------------------
//        for(auto &p:particles){
//            p.FLIP_velocity = {0,0};
//            p.FLIP_velocity = p.velocity;
//        }
//        for(int i=1;i<Nx+1;i++){
//            for(int j=0;j<Ny;j++){
//                //std::cout << u[i][j]-old_u[i][j] << ",";
//                std::vector<Eigen::Vector2d>gx_list = {{(i-1)*dx,(j+0.5)*dx},{(i)*dx,(j+0.5)*dx}};
//                std::vector<std::vector<int>>key_list = {{i-1,j},{i,j}};
//                if(extend){
//                    if(i != 1){
//                        gx_list.insert(gx_list.begin(),{(i-2)*dx,(j+0.5)*dx});
//                        key_list.insert(key_list.begin(),{i-2,j});
//                    }
//                    if(i != Nx){
//                        gx_list.push_back({(i+1)*dx,(j+0.5)*dx});
//                        key_list.push_back({i+1,j});
//                    }
//                }
//                for(int k=0;k<gx_list.size();k++){
//                    if(map.find(key_list[k]) != map.end()){
//                        auto val = map.at(key_list[k]);
//                        for(auto x:val){
//                            Eigen::Vector2d px = particles[x].position;
//                            double weight = weightFunction(px, gx_list[k], dx);
//    //                        if(i == 1)particles[x].FLIP_velocity.x()+= weight*mp*delta_u[i-1+k][j];
//    //                        else particles[x].FLIP_velocity.x()+= weight*mp*delta_u[i-2+k][j];
//                            if(i == 1)particles[x].FLIP_velocity.x()+= weight*mp*(u[i-1+k][j]-old_u[i-1+k][j]);
//                            else particles[x].FLIP_velocity.x()+= weight*mp*(u[i-2+k][j]-old_u[i-2+k][j]);
//                        }
//                    }
//                }
//            }
//            //std::cout << std::endl;
//        }
//        for(int i=0;i<Nx;i++)for(int j=1;j<Ny+1;j++){
//            std::vector<Eigen::Vector2d>gx_list = {{(i+0.5)*dx,(j)*dx},{(i+0.5)*dx,(j)*dx}};
//            std::vector<std::vector<int>>key_list = {{i,j-1},{i,j}};
//            if(extend){
//                if(j != 1){
//                    gx_list.insert(gx_list.begin(),{(i+0.5)*dx,(j-2)*dx});
//                    key_list.insert(key_list.begin(),{i,j-2});
//                }
//                if(j != Ny){
//                    gx_list.push_back({(i+0.5)*dx,(j+1)*dx});
//                    key_list.push_back({i,j+1});
//                }
//            }
//            for(int k=0;k<gx_list.size();k++){
//                if(map.find(key_list[k]) != map.end()){
//                    auto val = map.at(key_list[k]);
//                    for(auto x:val){
//                        Eigen::Vector2d px = particles[x].position;
//                        double weight = weightFunction(px, gx_list[k], dx);
//                        if(j == 1)particles[x].FLIP_velocity.y()+= weight*mp*(v[i][j-1+k]-old_v[i][j-1+k]);
//                        else particles[x].FLIP_velocity.y()+= weight*mp*(v[i][j-2+k]-old_v[i][j-2+k]);
//                    }
//                }
//            }
//        }
//    }
    void advectParticles(){
        //std::cout << L << std::endl;
        for(int i=0;i<particles.size();i++){
            //std::cout <<"x0:"<< particles[i].position.x() << " " << particles[i].position.y() << std::endl;
            particles[i].velocity = alpha*particles[i].FLIP_velocity + (1 - alpha)*particles[i].PIC_velocity;
            //if(particles[i].velocity.norm() * dt > dx)std::cout << "CFLError" << std::endl;
            particles[i].position.x() += particles[i].velocity.x()*dt;
            particles[i].position.y() += particles[i].velocity.y()*dt;
            particles[i].position.z() += particles[i].velocity.z()*dt;
            pushout(particles[i].position, L,dx);
            //std::cout << particles[i].velocity.x() << " " << particles[i].velocity.y() << std::endl;
        }
    }
    void preprocessingParticles(){
        locateParticlesOnGrid(map);
        //各粒子について，４近傍の密度分布を修正するベクトルの計算
        for(int i=0;i<particles.size();i++){
            std::vector<int>key = {particles[i].gridIndex[0],particles[i].gridIndex[1],particles[i].gridIndex[2]};
            std::vector<bool>F = {key[0]<Nx-1,key[1]<Ny-1,key[0]>0,key[1]>0,key[2]>0,key[2]<Nz-1};
            std::vector<std::vector<int>>keys = {{key[0]+1,key[1],key[2]},{key[0],key[1]+1,key[2]},{key[0]-1,key[1],key[2]},{key[0],key[1]-1,key[2]},{key[0],key[1],key[2]-1},{key[0],key[1],key[2]+1}};
            auto val = map.at(key);
            particles[i].fixVector = {0,0,0};
            for(auto &j:val){
                if(j == i)continue;
                particles[i].fixVector += fixDensityVector(particles[j].position, particles[i].position, radius, gamma, dt);
            }
            for(int k=0;k<4;k++){
                if(F[k] && map.find(keys[k]) != map.end()){
                    auto val = map.at(keys[k]);
                    //particles[i].fixVector = {0,0};
                    for(auto &j:val){
                        if(j == i)continue;
                        particles[i].fixVector += fixDensityVector(particles[j].position, particles[i].position, radius, gamma, dt);
                    }
                }
            }
//            std::cout << particles[i].fixVector.x() << "," << particles[i].fixVector.y() << std::endl;
        }
        //再サンプリングして速度を修正
        for(int i=0;i<particles.size();i++){
            particles[i].fixVelocity = {0,0,0};
            std::vector<int>key = {particles[i].gridIndex[0],particles[i].gridIndex[1],particles[i].gridIndex[2]};
            std::vector<bool>F = {key[0]<Nx-1,key[1]<Ny-1,key[0]>0,key[1]>0,key[2]>0,key[2]<Nz-1};
            std::vector<std::vector<int>>keys = {{key[0]+1,key[1],key[2]},{key[0],key[1]+1,key[2]},{key[0]-1,key[1],key[2]},{key[0],key[1]-1,key[2]},{key[0],key[1],key[2]-1},{key[0],key[1],key[2]+1}};
            auto val = map.at(key);
            Eigen::Vector3d sumA = {0,0,0};
            double sumB = 0;
            //bool flg = false;
            for(auto &j:val){
                if(j == i)continue;
                //std::cout << particles[i].fixVector.x() << "," << particles[i].fixVector.y() << std::endl;
                //fixParticleVelocity(particles[i], particles[j], radius, sumA, sumB, mp);
                sumB += fixParticleVelocity(particles[i], particles[j], radius, mp);
                sumA += fixParticleVelocity(particles[i], particles[j], radius, mp)*particles[j].velocity;
                //flg = true;
            }
            for(int k=0;k<4;k++){
                if(F[k] && map.find(keys[k]) != map.end()){
                    auto val = map.at(keys[k]);
                    for(auto &j:val){
                        if(j == i)continue;
                        //fixParticleVelocity(particles[i], particles[j], radius, sumA, sumB, mp);
                        sumB += fixParticleVelocity(particles[i], particles[j], radius, mp);
                        sumA += fixParticleVelocity(particles[i], particles[j], radius, mp)*particles[j].velocity;
                        //flg = true;
                        //std::cout << particles[j].velocity.x() << "," << particles[i].velocity.y()<< std::endl;
                    }
                }
            }
            if(sumB > 1.0e-4){
                //std::cout << sumA.x() << "," << sumA.y() << "," << sumB << std::endl;
                particles[i].fixVelocity = sumA/sumB;
                particles[i].velocity = particles[i].fixVelocity;
            }
        }
        //計算したベクトルで位置を修正
        for(int i=0;i<particles.size();i++){
            particles[i].position += particles[i].fixVector;
            //particles[i].velocity = particles[i].fixVelocity;
            pushout(particles[i].position, L,dx);
        }
        locateParticlesOnGrid(map);
    }
};
#endif /* Flip_h */
