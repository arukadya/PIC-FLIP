//
//  particle.hpp
//  Mable
//
//  Created by 須之内俊樹 on 2023/02/26.
//

#ifndef particle_hpp
#define particle_hpp
#include "Eigen/Core"
#include <stdio.h>
#include <vector>
struct Particle{
    Eigen::Vector3d PIC_velocity;
    Eigen::Vector3d FLIP_velocity;
    Eigen::Vector3d velocity;
    Eigen::Vector3d fixVector;
    Eigen::Vector3d fixVelocity;
    Eigen::Vector3d position;
    std::vector<int> gridIndex;
    Particle(Eigen::Vector3d &v,Eigen::Vector3d &p);
    void setGridIndex(int x,int y,int z);
};
#endif /* particle_hpp */
