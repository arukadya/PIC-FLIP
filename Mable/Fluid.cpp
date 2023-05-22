//
//  Fluid.cpp
//  Mable
//
//  Created by 須之内俊樹 on 2023/02/26.
//

#include "Fluid.hpp"
Fluid::Fluid(){
    
}
Fluid::Fluid(double x,double t,double density){
    dx = x;
    dt = t;
    rho = density;
    L = dx*Nx;
    vfi.reset(f0.y());
}

std::vector<int>Fluid::DirichletBoundaryCondition(int i,int j,int k,myMap &map){
    std::vector<int>ret(6,1);
    std::vector<std::vector<int>>keys = {{i+1,j,k},{i,j+1,k},{i-1,j,k},{i,j-1,k},{i,j,k-1},{i,j,k+1}};
    if(i == Nx-1)ret[0] = 0;
    else if(!map.contains(keys[0]))ret[0] = 0;
    if(j == Ny-1)ret[1] = 0;
    else if(!map.contains(keys[1]))ret[1] = 0;

    if(i == 0)ret[2] = 0;
    else if(!map.contains(keys[2]))ret[2] = 0;
    if(j == 0)ret[3] = 0;
    else if(!map.contains(keys[3]))ret[3] = 0;

    if(k == 0)ret[4] = 0;
    else if(!map.contains(keys[4]))ret[4] = 0;
    if(k == Nz-1)ret[5] = 0;
    else if(!map.contains(keys[5]))ret[5] = 0;
    return ret;
}
void Fluid::project(myMap &map){
    SparseMatrix A(Nx*Ny*Nz,Nx*Ny*Nz),B(Nx*Ny*Nz,Nx*Ny*Nz);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(Nx*Ny*Nz);
    Eigen::VectorXd px;
    std::set<int> DirichletKey;
    std::vector<std::vector<int>>keys;
    //Tripletの計算
    std::vector<Triplet> triplets;
    for(int i=0;i<Nx;i++){
        for(int j=0;j<Ny;j++){
            for(int k=0;k<Nz;k++){
                std::vector<int>key = {i,j,k};
                if(!map.contains(key)){
                    //前処理でAが変更されてしまうので，境界条件として別で無理矢理設定する．
                    triplets.emplace_back(i+j*Nx+k*Nx*Ny,i+j*Nx+k*Nx*Ny,1);
                    keys.push_back(key);
                    DirichletKey.insert(i+j*Nx+k*Nx*Ny);
                    continue;
                }
                double scale = dt/(rho*dx*dx);
                //std::cout << i << "," << j << std::endl;
                double D[6] = {1.0,1.0,-1.0,-1.0,-1.0,1.0};//周囲6方向に向かって働く、圧力の向き
                //double F[4] = {(double)(i<Nx-1),(double)(j<Ny-1),(double)(i>0),(double)(j>0)};//境界条件。壁なら0,流体なら1
                
                std::vector<int> F = {i<Nx-1,j<Ny-1,i>0,j>0,k>0,k<Nz-1};
                double U[6] = {
                    u.value[i+1][j][k],
                    v.value[i][j+1][k],
                    u.value[i][j][k],
                    v.value[i][j][k],
                    w.value[i][j][k],
                    w.value[i][j][k+1]};
                double sumP = 0;
                for(int n=0;n<6;n++){
                    sumP += -F[n]*scale;
                    //sumP += scale;
                    b(i+j*Nx+k*Nx*Ny) += D[n]*F[n]*U[n]/(dx);
                }
                F = DirichletBoundaryCondition(i,j,k,map);
//                    for(int n=0;n<6;n++){
//                        std::cout << F_pri[n] << "," << F[n] << std::endl;
//                    }
                triplets.emplace_back(i+j*Nx+k*Nx*Ny,i+j*Nx+k*Nx*Ny, sumP);
                if(F[0])triplets.emplace_back(i+j*Nx+k*Nx*Ny,i+1+j*Nx+k*Nx*Ny, F[0]*scale);
                if(F[1])triplets.emplace_back(i+j*Nx+k*Nx*Ny,i+(j+1)*Nx+k*Nx*Ny, F[1]*scale);
                if(F[2])triplets.emplace_back(i+j*Nx+k*Nx*Ny,i-1+j*Nx+k*Nx*Ny, F[2]*scale);
                if(F[3])triplets.emplace_back(i+j*Nx+k*Nx*Ny,i+(j-1)*Nx+k*Nx*Ny, F[3]*scale);
                if(F[4])triplets.emplace_back(i+j*Nx+k*Nx*Ny,i+j*Nx+(k-1)*Nx*Ny, F[4]*scale);
                if(F[5])triplets.emplace_back(i+j*Nx+k*Nx*Ny,i+j*Nx+(k+1)*Nx*Ny, F[5]*scale);
            }
        }
    }
    A.setFromTriplets(triplets.begin(), triplets.end());
    Eigen::ConjugateGradient<SparseMatrix> solver;
    //Eigen::BiCGSTAB<SparseMatrix> solver;
    
    for(int i=0;i<A.outerSize();++i){
        for(SparseMatrix::InnerIterator it(A,i);it;++it){
            if(it.row() == *DirichletKey.begin()){
                it.valueRef() = 1;
                //std::cout << "DirichletKey.size()=" << DirichletKey.size() << std::endl;
                DirichletKey.erase(DirichletKey.begin());
            }
        }
    }
    solver.compute(A);
    px = solver.solve(b);
    for(int i=0;i<Nx;i++){
        for(int j=0;j<Ny;j++){
            for(int k=0;k<Nz;k++){
                p.value[i][j][k] = px(i+j*Nx+k*Nx*Ny);
            }
        }
    }
//    for(auto x:keys){
//        if(px(x[0] + x[1]*Ny + x[2]*Nz*Nz) > 1.0e-4)std::cout << x[0] + x[1]*Ny + x[2]*Nz*Nz << "," << px(x[0] + x[1]*Ny + x[2]*Nz*Nz) <<std::endl;
        //else std::cout << "success:" << x[0] + x[1]*Ny + x[2]*Nz*Nz << std::endl;
//   }
    
    for(int i=1; i<Nx;i++){
        for(int j=0;j<Ny;j++){
            for(int k=0;k<Nz;k++)u.value[i][j][k] = u.value[i][j][k] - dt/rho * (p.value[i][j][k]-p.value[i-1][j][k])/dx;
        }
    }
    for(int i=0;i<Nx;i++){
        for(int j=1;j<Ny;j++){
            for(int k=0;k<Nz;k++)v.value[i][j][k] = v.value[i][j][k] - dt/rho * (p.value[i][j][k]-p.value[i][j-1][k])/dx;
        }
    }
    for(int i=0;i<Nx;i++){
        for(int j=0;j<Ny;j++){
            for(int k=1;k<Nz;k++)w.value[i][j][k] = w.value[i][j][k] - dt/rho * (p.value[i][j][k]-p.value[i][j][k-1])/dx;
        }
    }
}

