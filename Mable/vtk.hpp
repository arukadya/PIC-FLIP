//
//  vtk.hpp
//  Mable
//
//  Created by 須之内俊樹 on 2023/02/26.
//

#ifndef vtk_hpp
#define vtk_hpp
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <Eigen/Core>
#include "Array3d.hpp"
void inputVTK(const char* InputFileName,std::vector<Eigen::Vector3d> &Vertices,std::vector<std::vector<int>> &Faces);

void outputVTK(std::string OutputFileName,myArray3d &Vertices,double dx,int nx,int ny,int nz);

void outputVTK_implicit(std::string OutputFileName,myArray3d &Vertices,double dx,double threshold);

void outputVolume(std::string OutputFileName,std::vector<double> &volumes);
#endif /* vtk_hpp */
