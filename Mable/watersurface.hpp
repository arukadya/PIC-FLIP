//
//  watersurface.hpp
//  Mable
//
//  Created by 須之内俊樹 on 2023/02/26.
//

#ifndef watersurface_hpp
#define watersurface_hpp

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
#include "Array3d.hpp"
#include <unordered_map>
#include "Hasher.h"
#include "functions.hpp"
myArray3d cal_implicitFunction(std::vector<Particle> &particles,myMap &map,double radius,double dx,int nx,int ny,int nz);

std::vector<std::vector<double>> makeSurface(std::vector<Particle> &particles,myMap &map,double radius,double dx,int nx,int ny,int nz,double threshold,double th_d);

std::vector<std::vector<double>> makeSurface_imp(std::vector<Particle> &particles,myMap &map,double radius,double dx,int nx,int ny,int nz);

enum{
    inWater = 0,
    inAir = 1
};
double cal_volume(std::vector<Particle> &particles,myMap &map,double radius,double dx,int nx,int ny,int nz,double threshold,int cnt);
#endif /* watersurface_hpp */
