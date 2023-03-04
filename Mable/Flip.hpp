//
//  Flip.hpp
//  Mable
//
//  Created by 須之内俊樹 on 2023/02/26.
//

#ifndef Flip_hpp
#define Flip_hpp

#include "Hasher.h"//Hash関数
#include <functional>//std::hash
#include <type_traits>//std::remove_cvref_t(C++20)

#include "Fluid.hpp"
#include "vtk.hpp"
#include "functions.hpp"
#include "particle.hpp"
#include "watersurface.hpp"
#define repeatCount 2000
#define alpha 1
#ifndef Flip_h
#define Flip_h
#define timer 0
#define extend 0
#define th_d 0.1
struct Fluid;
struct PIC_FLIP : Fluid{
    std::vector<int> division;//division[0] = xの分割数.division[1]=y...
    std::vector<Particle>particles;//入力メッシュの頂点
    std::unordered_map<std::vector<int>,std::vector<int>,ArrayHasher<3>>map;//ハッシュテーブル
    std::vector<Eigen::Vector3d> vertices;//出力メッシュの頂点
    myArray3d implicit_function = myArray3d(Nx,Ny,Nz,0);
    std::vector<std::vector<double>> surfaceMesh;
    std::vector<double>volumes;
    TimeDisplayer TD;
    double radius;
    double mp;
    double gamma;
    double threshold;
    PIC_FLIP(double x,double t,double density,double r,double th,double ga):Fluid(x,t,density){
        radius = dx/2*3;
        mp = pow(radius,3)/3*4*3.14;
        gamma = 1;
        threshold = 1;
        division = {Nx,Ny,Nz};
        initParticles();
        radius = r;
        threshold = th;
        gamma = ga;
    };
    void execute(std::string foldername,std::string filename);
    void output(std::vector<Eigen::Vector3d> &v);
    void initParticles();
    void locateParticlesOnGrid(std::unordered_map<std::vector<int>,std::vector<int>,ArrayHasher<3>>&map);
    Eigen::Vector3d calCellSize();
    std::vector<int> calGridOrKey(Eigen::Vector3d v);
    void particlesVelocityToGrid();
    void calGridPressure();
    void gridVelocityToParticles();
    void advectParticles();
    void preprocessingParticles();
};
#endif /* Flip_h */

#endif /* Flip_hpp */
