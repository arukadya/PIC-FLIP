//
//  ioOFF.hpp
//  Mable
//
//  Created by 須之内俊樹 on 2023/03/16.
//

#ifndef ioOFF_hpp
#define ioOFF_hpp
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <Eigen/Core>
#include <fstream>
void inputOFF(const char* InputFileName,std::vector<Eigen::Vector3d> &Vertices,std::vector<std::vector<int>> &Faces);

void outputOFF(std::string OutputFileName,std::vector<Eigen::Vector3d> &Vertices,std::vector<std::vector<int>> &Faces);

#endif /* ioOFF_hpp */
