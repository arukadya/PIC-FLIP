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
            }
        }
    }
    return implicit_function;
}
void makeSurface(std::vector<particle> &particles,std::unordered_map<std::vector<int>,std::vector<int>,ArrayHasher<3>>&map,double radius,double dx,int nx,int ny,int nz,double threshold){
    myArray3d implicit_function = cal_implicitFunction(particles, map, radius, dx, nx, ny, nz);
    std::vector<myArray3d> outputMesh;
    for(int i=0;i<nx;i++){
        for(int j=0;j<ny;j++){
            for(int k=0;k<nz;k++){
                if(implicit_function.value[i][j][k] > threshold){
                    //outputMesh.push_back({(i+0.5)*dx,(j+0.5)*dx,(k+0.5)*dx});
                }
            }
        }
    }
}
#endif /* watersurface_h */
