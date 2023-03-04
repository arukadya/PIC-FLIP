//
//  functions.hpp
//  Mable
//
//  Created by 須之内俊樹 on 2023/02/26.
//

#ifndef functions_hpp
#define functions_hpp
#include <chrono>
#include "Eigen/Core"
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <math.h>
#include <sstream>
#include <iomanip>
#include <queue>
#include <iostream>
#include "particle.hpp"
struct TimeDisplayer{
    std::chrono::system_clock::time_point startTime;
    std::chrono::system_clock::time_point endTime;
    const char* str;
    void startTimer(const char* s);
    void endTimer();
};
double kernelFunction(double x);
double sharp_kernel(double x);
double weightFunction(Eigen::Vector3d &px,Eigen::Vector3d &gx,double dx);
double smoothFunction(Eigen::Vector3d r,double dx);
double sharpFunction(Eigen::Vector3d r,double dx);
Eigen::Vector3d fixDensityVector(Eigen::Vector3d &pj,Eigen::Vector3d &pi,double dx,double gamma,double dt);

double fixParticleVelocity(Particle &pi,Particle &pj,double dx,double mp);
void pushout(Eigen::Vector3d &x,double L,double dx);
#endif /* functions_hpp */
