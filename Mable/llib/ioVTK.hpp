//
//  ioVTK.hpp
//  marchining_cubes
//
//  Created by 須之内俊樹 on 2023/03/04.
//

#ifndef ioVTK_hpp
#define ioVTK_hpp

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <math.h>
#include <sstream>
#include <Eigen/Core>
int inputVTK(int &num_cells,std::vector<int> &nums,Eigen::Vector3d &dists,Eigen::Vector3d &origin,std::vector<double> &v,std::string InputFlieName);
#endif /* ioVTK_hpp */
