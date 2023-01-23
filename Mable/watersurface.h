//
//  watersurface.h
//  Mable
//
//  Created by 須之内俊樹 on 2023/01/19.
//

#ifndef watersurface_h
#define watersurface_h

#include "particle.h"

myArray3d cal_implicitFunction(std::vector<particle> &particles,std::unordered_map<std::vector<int>,std::vector<int>,ArrayHasher<3>>&map,double radius,double dx,int nx,int ny,int nz){
    myArray3d implicit_function = myArray3d(nx,ny,nz,0);
    for(int i=0;i<nx;i++){
        for(int j=0;j<ny;j++){
            for(int k=0;k<nz;k++){
                Eigen::Vector3d gridPos = {(i+0.5)*dx,(j+0.5)*dx,(k+0.5)*dx};
                std::vector<std::vector<int>>key_list = {{i-1,j,k},{i,j,k}};
                double func_val = 0;
                for(int l=0;l<key_list.size();l++){
                    if(map.find(key_list[l]) != map.end()){
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
std::vector<std::vector<double>> makeSurface(std::vector<particle> &particles,std::unordered_map<std::vector<int>,std::vector<int>,ArrayHasher<3>>&map,double radius,double dx,int nx,int ny,int nz,double threshold,double th_d){
    myArray3d implicit_function = cal_implicitFunction(particles, map, radius, dx, nx, ny, nz);
    std::vector<std::vector<double>> outputMesh;
    for(int i=0;i<nx;i++){
        for(int j=0;j<ny;j++){
            for(int k=0;k<nz;k++){
                if(implicit_function.value[i][j][k] < threshold + th_d &&implicit_function.value[i][j][k] > threshold){
                    outputMesh.push_back({(i+0.5)*dx,(j+0.5)*dx,(k+0.5)*dx});
                }
            }
        }
    }
    outputMesh.push_back({(0.5)*dx,(0.5)*dx,(0.5)*dx});
    outputMesh.push_back({(nx+0.5)*dx,(0.5)*dx,(0.5)*dx});
    outputMesh.push_back({(nx+0.5)*dx,(0.5)*dx,(nz+0.5)*dx});
    outputMesh.push_back({(0.5)*dx,(0.5)*dx,(nz+0.5)*dx});
    return outputMesh;
}
std::vector<std::vector<double>> makeSurface_imp(std::vector<particle> &particles,std::unordered_map<std::vector<int>,std::vector<int>,ArrayHasher<3>>&map,double radius,double dx,int nx,int ny,int nz){
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
enum {
    inWater = 0,
    inAir = 1
};
double cal_volume(std::vector<particle> &particles,std::unordered_map<std::vector<int>,std::vector<int>,ArrayHasher<3>>&map,double radius,double dx,int nx,int ny,int nz,double threshold,int cnt){
    myArray3d implicit_function = cal_implicitFunction(particles, map, radius, dx, nx, ny, nz);
    double sum = 0;
    for(int i=0;i<nx;i++){
        for(int k=0;k<nz;k++){
            bool flg = true;
            bool cross = false;
            int c = 0;
            for(int j=0;j<ny;j++){
                if(implicit_function.value[i][j][k] < threshold){
                    //if()
                    sum += (j+0.5)*dx;
                    if(cross)c++;
                    cross = false;
                    //if(cnt < 2)std::cout << "(" << i <<"," << j << "," << k << ") = " <<implicit_function.value[i][j][k] << std::endl;
                    //break;
                }
                else{
                    if(!cross)c++;
                    cross = true;
                }
            }
            //std::cout << "(" << i << "," << k << ") = " << c << std::endl;
            if(flg){
                sum += (ny-0.5)*dx;
                //if(cnt < 2)std::cout << "(" << i <<"," << ny-1 << "," << k << ") = " <<implicit_function.value[i][ny-1][k] << std::endl;
            }
        }
    }
    return sum;
}
#endif /* watersurface_h */
