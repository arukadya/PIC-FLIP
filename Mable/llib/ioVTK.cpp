//
//  ioVTK.cpp
//  marchining_cubes
//
//  Created by 須之内俊樹 on 2023/03/04.
//

#include "ioVTK.hpp"
int inputVTK(int &num_cells,std::vector<int> &nums,Eigen::Vector3d &dists,Eigen::Vector3d &origin,std::vector<double> &v,std::string InputFlieName){
    std::ifstream Inputfile(InputFlieName);
    if (!Inputfile.is_open()) {
        std::cerr << "Could not open the file - '"
             << InputFlieName << "'" << std::endl;
        return EXIT_FAILURE;
    }
    std::string line;
    std::string word;
    for(int i=0;i<4;++i)std::getline(Inputfile,line);//4行は無視
    
    std::getline(Inputfile,line);
    std::stringstream ss_nums{line};
    getline(ss_nums,word,' ');
    for(int i=0;i<3;++i){
        getline(ss_nums,word,' ');
        nums.push_back(std::stoi(word));
    }
    //std::cout << nums.size() << std::endl;
    std::getline(Inputfile,line);
    getline(ss_nums,word,' ');
    double x = std::stod(word);
    getline(ss_nums,word,' ');
    double y = std::stod(word);
    getline(ss_nums,word,' ');
    double z = std::stod(word);
    origin = {x,y,z};
    
    std::getline(Inputfile,line);
    std::stringstream ss_dists{line};
    getline(ss_dists,word,' ');
    getline(ss_dists,word,' ');
    x = std::stod(word);
    getline(ss_dists,word,' ');
    y = std::stod(word);
    getline(ss_dists,word,' ');
    z = std::stod(word);
    dists = {x,y,z};
    
    std::getline(Inputfile,line);
    std::stringstream ss_num{line};
    getline(ss_num,word,' ');
    getline(ss_num,word,' ');
    num_cells = std::stoi(word);
    for(int i=0;i<2;++i)std::getline(Inputfile,line);//2行は無視
    
    while(getline(Inputfile,line)){
        std::stringstream ss{line};
        getline(ss,word);
        v.push_back(std::stod(word));
    }
    Inputfile.close();
    return EXIT_SUCCESS;
}
