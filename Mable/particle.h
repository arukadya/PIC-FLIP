//
//  particle.h
//  Mable
//
//  Created by 須之内俊樹 on 2022/12/24.
//

#ifndef particle_h
#define particle_h

struct particle{
    Eigen::Vector2d PIC_velocity;
    Eigen::Vector2d FLIP_velocity;
    Eigen::Vector2d velocity;
    Eigen::Vector2d fixVector;
    Eigen::Vector2d fixVelocity;
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

#endif /* particle_h */
