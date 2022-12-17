//
//  HashMapClasstaring.h
//  HashMapClasstaring
//
//  Created by 須之内俊樹 on 2022/10/29.
//
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
#define repeatCount 100
//#define Np 100 //粒子の数
#ifndef Flip_hpp
#define Flip_hpp
double kernelFunction(double x);
double weightFunction(Eigen::Vector2d px,Eigen::Vector2d gx,double dx);
#endif /* Flip_hpp */
