//
//  particle.cpp
//  Mable
//
//  Created by 須之内俊樹 on 2023/02/26.
//

#include "particle.hpp"

Particle::Particle(Eigen::Vector3d &v,Eigen::Vector3d &p){
    gridIndex = {-1, -1, -1};
    Eigen::Vector3d zero = {0,0,0};
    velocity = v;
    PIC_velocity = zero;
    FLIP_velocity = zero;
    position = p;
}
void Particle::setGridIndex(int x,int y,int z){
    gridIndex[0] = x;
    gridIndex[1] = y;
    gridIndex[2] = z;
}

