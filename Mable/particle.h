//
//  particle.h
//  Mable
//
//  Created by 須之内俊樹 on 2022/12/24.
//

#ifndef particle_h
#define particle_h

struct Particle{
    Eigen::Vector3d PIC_velocity;
    Eigen::Vector3d FLIP_velocity;
    Eigen::Vector3d velocity;
    Eigen::Vector3d fixVector;
    Eigen::Vector3d fixVelocity;
    Eigen::Vector3d position;
    std::vector<int> gridIndex {-1, -1, -1};
    Particle(Eigen::Vector3d &v,Eigen::Vector3d &p);
    void setGridIndex(int x,int y,int z);
};

#endif /* particle_h */
