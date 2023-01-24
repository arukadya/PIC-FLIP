#include "mable.h"
#include "testdata.h"
#include "Flip.h"
int main(int argc, const char * argv[]) {
    double dx = 0.1;//セルの大きさ
    double dt = 0.01;//時間の刻み幅
    double rho = 1.0;
    //r = 2
    PIC_FLIP simulator_r2 = PIC_FLIP(dx,dt,rho,dx*2,1,1);//radius,threshold,gamma
    simulator_r2.execute("outputVTK/r2dx_th1_ga1", "volumes/volume_r2dx_th1_ga1.dat");
    //threshold
    PIC_FLIP simulator_r2_th05 = PIC_FLIP(dx,dt,rho,dx*2,1.05,1);//radius,threshold,gamma
    simulator_r2_th05.execute("outputVTK/r2dx_th05_ga1", "volumes/volume_r2dx_th05_ga1.dat");
    PIC_FLIP simulator_r2_th10 = PIC_FLIP(dx,dt,rho,dx*2,1.1,1);//radius,threshold,gamma
    simulator_r2_th10.execute("outputVTK/r2dx_th10_ga1", "volumes/volume_r2dx_th10_ga1.dat");
    PIC_FLIP simulator_r2_th95 = PIC_FLIP(dx,dt,rho,dx*2,0.95,1);//radius,threshold,gamma
    simulator_r2_th95.execute("outputVTK/r2dx_th95_ga1", "volumes/volume_r2dx_th95_ga1.dat");
    PIC_FLIP simulator_r2_th90 = PIC_FLIP(dx,dt,rho,dx*2,0.90,1);//radius,threshold,gamma
    simulator_r2_th90.execute("outputVTK/r2dx_th90_ga1", "volumes/volume_r2dx_th90_ga1.dat");
    //gamma
    PIC_FLIP simulator_r2_ga2 = PIC_FLIP(dx,dt,rho,dx*2,1,2);//radius,threshold,gamma
    simulator_r2_ga2.execute("outputVTK/r2dx_th1_ga2", "volumes/volume_r2dx_th1_ga2.dat");
    PIC_FLIP simulator_r2_ga15 = PIC_FLIP(dx,dt,rho,dx*2,1,1.5);//radius,threshold,gamma
    simulator_r2_ga15.execute("outputVTK/r2dx_th1_ga15", "volumes/volume_r2dx_th1_ga15.dat");
    PIC_FLIP simulator_r2_ga05 = PIC_FLIP(dx,dt,rho,dx*2,1,0.5);//radius,threshold,gamma
    simulator_r2_ga05.execute("outputVTK/r2dx_th1_ga05", "volumes/volume_r2dx_th1_ga05.dat");
    //r = 1.5
    PIC_FLIP simulator_r15 = PIC_FLIP(dx,dt,rho,dx*3/2,1,1);//radius,threshold,gamma
    simulator_r15.execute("outputVTK/r15dx_th1_ga1", "volumes/volume_r15dx_th1_ga1.dat");
    //threshold
    PIC_FLIP simulator_r15_th05 = PIC_FLIP(dx,dt,rho,dx*3/2,1.05,1);//radius,threshold,gamma
    simulator_r15_th05.execute("outputVTK/r15dx_th05_ga1", "volumes/volume_r15dx_th05_ga1.dat");
    PIC_FLIP simulator_r15_th10 = PIC_FLIP(dx,dt,rho,dx*3/2,1.1,1);//radius,threshold,gamma
    simulator_r15_th10.execute("outputVTK/r15dx_th10_ga1", "volumes/volume_r15dx_th10_ga1.dat");
    PIC_FLIP simulator_r15_th95 = PIC_FLIP(dx,dt,rho,dx*3/2,0.95,1);//radius,threshold,gamma
    simulator_r15_th95.execute("outputVTK/r15dx_th95_ga1", "volumes/volume_r15dx_th95_ga1.dat");
    PIC_FLIP simulator_r15_th90 = PIC_FLIP(dx,dt,rho,dx*3/2,0.90,1);//radius,threshold,gamma
    simulator_r15_th90.execute("outputVTK/r15dx_th90_ga1", "volumes/volume_r15dx_th90_ga1.dat");
    //gamma
    PIC_FLIP simulator_r15_ga2 = PIC_FLIP(dx,dt,rho,dx*3/2,1,2);//radius,threshold,gamma
    simulator_r15_ga2.execute("outputVTK/r15dx_th1_ga2", "volumes/volume_r15dx_th1_ga2.dat");
    PIC_FLIP simulator_r15_ga15 = PIC_FLIP(dx,dt,rho,dx*3/2,1,1.5);//radius,threshold,gamma
    simulator_r15_ga15.execute("outputVTK/r15dx_th1_ga15", "volumes/volume_r15dx_th1_ga15.dat");
    PIC_FLIP simulator_r15_ga05 = PIC_FLIP(dx,dt,rho,dx*3/2,1,0.5);//radius,threshold,gamma
    simulator_r15_ga05.execute("outputVTK/r15dx_th1_ga05", "volumes/volume_r15dx_th1_ga05.dat");
    
    //r=1
    PIC_FLIP simulator_r1 = PIC_FLIP(dx,dt,rho,dx*1,1,1);//radius,threshold,gamma
    simulator_r1.execute("outputVTK/r1dx_th1_ga1", "volumes/volume_r1dx_th1_ga1.dat");
    //threshold
    PIC_FLIP simulator_r1_th05 = PIC_FLIP(dx,dt,rho,dx ,1.05,1);//radius,threshold,gamma
    simulator_r1_th05.execute("outputVTK/r1dx_th05_ga1", "volumes/volume_r1dx_th05_ga1.dat");
    PIC_FLIP simulator_r1_th10 = PIC_FLIP(dx,dt,rho,dx ,1.1,1);//radius,threshold,gamma
    simulator_r1_th10.execute("outputVTK/r1dx_th10_ga1", "volumes/volume_r1dx_th10_ga1.dat");
    PIC_FLIP simulator_r1_th95 = PIC_FLIP(dx,dt,rho,dx ,0.95,1);//radius,threshold,gamma
    simulator_r1_th95.execute("outputVTK/r1dx_th95_ga1", "volumes/volume_r1dx_th95_ga1.dat");
    PIC_FLIP simulator_r1_th90 = PIC_FLIP(dx,dt,rho,dx ,0.90,1);//radius,threshold,gamma
    simulator_r1_th90.execute("outputVTK/r1dx_th90_ga1", "volumes/volume_r1dx_th90_ga1.dat");
    //gamma
    PIC_FLIP simulator_r1_ga2 = PIC_FLIP(dx,dt,rho,dx ,1,2);//radius,threshold,gamma
    simulator_r1_ga2.execute("outputVTK/r1dx_th1_ga2", "volumes/volume_r1dx_th1_ga2.dat");
    PIC_FLIP simulator_r1_ga15 = PIC_FLIP(dx,dt,rho,dx ,1,1.5);//radius,threshold,gamma
    simulator_r1_ga15.execute("outputVTK/r1dx_th1_ga15", "volumes/volume_r1dx_th1_ga15.dat");
    PIC_FLIP simulator_r1_ga05 = PIC_FLIP(dx,dt,rho,dx ,1,0.5);//radius,threshold,gamma
    simulator_r1_ga05.execute("outputVTK/r1dx_th1_ga05", "volumes/volume_r1dx_th1_ga05.dat");
    return 0;
}
