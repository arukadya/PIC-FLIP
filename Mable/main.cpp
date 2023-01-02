#include "mable.h"
#include "testdata.h"
#include "Flip.h"
int main(int argc, const char * argv[]) {
    double dx = 0.1;//セルの大きさ
    double dt = 0.1;//時間の刻み幅
    double rho = 1.0;

    std::vector<std::vector<double>>u(Nx+1);
    std::vector<std::vector<double>>umi(Nx+1);
    std::vector<std::vector<double>>v(Nx);
    std::vector<std::vector<double>>vmi(Nx);
    std::vector<std::vector<double>>p(Nx);
    std::vector<std::vector<double>>ufi(Nx+1);
    std::vector<std::vector<double>>vfi(Nx);
    
    for(int i=0;i<Nx+1;i++)u[i] = std::vector<double>(Ny,0);
    for(int i=0;i<Nx+1;i++)umi[i] = std::vector<double>(Ny,0);
    for(int i=0;i<Nx+1;i++)ufi[i] = std::vector<double>(Ny,0);
    for(int i=0;i<Nx;i++)v[i] = std::vector<double>(Ny+1,0);
    for(int i=0;i<Nx;i++)vmi[i] = std::vector<double>(Ny+1,0);
    for(int i=0;i<Nx;i++)vfi[i] = std::vector<double>(Ny+1,0);
    for(int i=0;i<Nx;i++)p[i] = std::vector<double>(Ny,0);
    
    PIC_FLIP simulator = PIC_FLIP(dx,dt,rho,u,v,p,umi,vmi,ufi,vfi);
    simulator.execute();
    return 0;
}
