#include "mable.h"
#include "testdata.h"
#include "Flip.h"
#define Margin 0
int main(int argc, const char * argv[]) {
    double dx = 0.01;//セルの大きさ
    double dt = 0.1;//時間の刻み幅
    double rho = 1.0;
//    double **u;
//    double **v;
//    double **p;
    std::vector<std::vector<double>>u(Nx+1);
    std::vector<std::vector<double>>umi(Nx+1);
    std::vector<std::vector<double>>v(Nx);
    std::vector<std::vector<double>>vmi(Nx);
    std::vector<std::vector<double>>p(Nx);
    std::vector<std::vector<Eigen::Vector2d>>fi(Nx);
    
    for(int i=0;i<Nx+1;i++)u[i] = std::vector<double>(Ny);
    for(int i=0;i<Nx+1;i++)umi[i] = std::vector<double>(Ny);
    for(int i=0;i<Nx;i++)v[i] = std::vector<double>(Ny+1);
    for(int i=0;i<Nx;i++)vmi[i] = std::vector<double>(Ny+1);
    for(int i=0;i<Nx;i++)p[i] = std::vector<double>(Ny);
    for(int i=0;i<Nx;i++)fi[i] = std::vector<Eigen::Vector2d>(Ny);
    
    //二次元配列のアドレス確保
//    u = (double **)malloc(sizeof(double *)*Ny+Margin);
//    v = (double **)malloc(sizeof(double *)*Ny+1+Margin);
//    p = (double **)malloc(sizeof(double *)*Ny+Margin);
//    for(unsigned int i=0;i<Ny+Margin;i++)u[i] =(double *)malloc(sizeof(double)*Nx+1+Margin);
//    for(unsigned int i=0;i<Ny+1+Margin;i++)v[i]=(double *)malloc(sizeof(double)*Nx+Margin);
//    for(unsigned int i=0;i<Ny+Margin;i++)p[i] =(double *)malloc(sizeof(double)*Nx+Margin);
//
    PIC_FLIP simulator = PIC_FLIP(dx,dt,rho,u,v,p,umi,vmi,fi);
    //std::cout << "simulator initialized!\n";
    simulator.execute();
//    for(unsigned int i=0;i<Ny+Margin;i++)free(u[i]);
//    free(u);
//    for(unsigned int i=0;i<Ny+1+Margin;i++)free(v[i]);
//    free(v);
//
    return 0;
}
