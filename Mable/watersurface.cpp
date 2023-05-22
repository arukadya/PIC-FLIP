//
//  watersurface.cpp
//  Mable
//
//  Created by 須之内俊樹 on 2023/02/26.
//

#include "watersurface.hpp"

myArray3d cal_implicitFunction(std::vector<Particle> &particles,myMap &map,double radius,double dx,int nx,int ny,int nz){
    myArray3d implicit_function = myArray3d(nz+1,ny+1,nx+1,0);
    for(int i=1;i<nz+1;i++){
        for(int j=1;j<ny+1;j++){
            for(int k=1;k<nx+1;k++){
                Eigen::Vector3d gridPos = {(k+0.5)*dx,(j+0.5)*dx,(i+0.5)*dx};
                std::vector<std::vector<int>>key_list = {{k-1,j,i},{k,j,i}};
                double func_val = 0;
                for(int l=0;l<key_list.size();l++){
                    if(map.contains(key_list[l])){
                        auto val = map.at(key_list[l]);
                        for(auto x:val){
                            Eigen::Vector3d px = particles[x].position;
                            Eigen::Vector3d pv = particles[x].velocity;
                            double weight = weightFunction(px, gridPos, dx);
                            func_val += weight;
                        }
                    }
                }
                implicit_function.value[i][j][k] = func_val;
                //std::cout << func_val << std::endl;
            }
        }
    }
    return implicit_function;
}
std::vector<std::vector<double>> makeSurface(std::vector<Particle> &particles,myMap &map,double radius,double dx,int nx,int ny,int nz,double threshold,double th_d){
    myArray3d implicit_function = cal_implicitFunction(particles, map, radius, dx, nx, ny, nz);
    std::vector<std::vector<double>> outputMesh;
    for(int i=0;i<nx;i++){
        for(int j=0;j<ny;j++){
            for(int k=0;k<nz;k++){
//                if(implicit_function.value[i][j][k] < threshold + th_d &&implicit_function.value[i][j][k] > threshold){
                outputMesh.push_back({(i+0.5)*dx,(j+0.5)*dx,(k+0.5)*dx,implicit_function.value[i][j][k]});
//                }
            }
        }
    }
    return outputMesh;
}
std::vector<std::vector<double>> makeSurface_imp(std::vector<Particle> &particles,myMap &map,double radius,double dx,int nx,int ny,int nz){
    myArray3d implicit_function = cal_implicitFunction(particles, map, radius, dx, nx, ny, nz);
    std::vector<std::vector<double>> outputMesh;
    for(int i=0;i<nx;i++){
        for(int j=0;j<ny;j++){
            for(int k=0;k<nz;k++){
                outputMesh.push_back({(i+0.5)*dx,(j+0.5)*dx,(k+0.5)*dx,implicit_function.value[i][j][k]});
            }
        }
    }
//    outputMesh.push_back({(0.5)*dx,(0.5)*dx,(0.5)*dx});
//    outputMesh.push_back({(nx+0.5)*dx,(0.5)*dx,(0.5)*dx});
//    outputMesh.push_back({(nx+0.5)*dx,(0.5)*dx,(nz+0.5)*dx});
//    outputMesh.push_back({(0.5)*dx,(0.5)*dx,(nz+0.5)*dx});
    return outputMesh;
}
double cal_volume(std::vector<Particle> &particles,myMap &map,double radius,double dx,int nx,int ny,int nz,double threshold,int cnt){
    myArray3d implicit_function = cal_implicitFunction(particles, map, radius, dx, nz, ny, nx);
    double sum = 0;
    for(int i=0;i<nx;i++){
        for(int k=0;k<nz;k++){
            double WS_under = 0;
            int state = inWater;
            for(int j=0;j<ny;j++){
                if(j == 0){
                    if(implicit_function.value[i][j][k] < threshold)state = inAir;
                    else{
                        state = inWater;
                        WS_under = (j+0.5)*dx;
                    }
                }
                else{
                    if(implicit_function.value[i][j][k] < threshold){
                        if(state == inWater){
                            state = inAir;
                            sum += (j+0.5)*dx - WS_under;
                        }
                        //if(cnt < 2)std::cout << "(" << i <<"," << j << "," << k << ") = " <<implicit_function.value[i][j][k] << std::endl;
                        //break;
                    }
                    else{
                        if(state == inAir){
                            state = inWater;
                            WS_under = (j+0.5)*dx;
                        }
                    }
                    if(j == ny-1  && implicit_function.value[i][j][k] > threshold){
                        sum += (j-0.5)*dx - WS_under;
                    }
                }
            }
            //std::cout << "(" << i << "," << k << ") = " << c << std::endl;
        }
    }
    return sum;
}
